#!/usr/bin/env python3
"""
of_reader.py - read an OpenFOAM case (ASCII) with numpy only.

  * constant/polyMesh/{points,faces,owner,neighbour,boundary}, plain or .gz
  * field files (internalField uniform / nonuniform List<scalar|vector>)
  * mesh geometry for general polyhedral meshes: face centres and area vectors,
    cell centres and volumes (face-pyramid decomposition, as OpenFOAM's
    primitiveMesh does it), closedness of every cell
  * cell shapes for VTK: hexahedra (type 12, VTK node order derived from the
    faces) and face streams for general polyhedra (type 42)

Files with `format binary` in the FoamFile header are refused with a clear
message (convert with `foamFormatConvert` after setting writeFormat ascii).

Python >= 3.6, numpy. No per-cell or per-face Python loops.
"""
import gzip
import os
import re
import sys
from pathlib import Path

import numpy as np

__all__ = [
    "FoamFormatError", "find_file", "read_raw", "parse_foam_header",
    "parse_internal_field", "expand_uniform", "read_field", "read_label_list",
    "read_cell_labels", "get_time_directories", "time_name", "PolyMesh",
    "read_polymesh",
]


class FoamFormatError(RuntimeError):
    pass


# =============================================================================
# low level: files, header, number lists
# =============================================================================
def find_file(path):
    """Return the existing file `path` or `path.gz`, else None."""
    p = Path(path)
    if p.is_file():
        return p
    g = Path(str(p) + ".gz")
    if g.is_file():
        return g
    return None


def read_raw(path):
    """Whole file as bytes; `path` or `path.gz` (gzip is detected by suffix
    and by magic number)."""
    p = find_file(path)
    if p is None:
        raise FileNotFoundError("no such file: %s (also tried %s.gz)" % (path, path))
    with open(str(p), "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw


_RE_FOAMFILE = re.compile(rb"FoamFile\s*\{(.*?)\}", re.S)
_RE_ENTRY = re.compile(rb"^\s*(\w+)\s+(.*?);\s*$", re.M)


def parse_foam_header(raw, path="?", require_ascii=True):
    """Entries of the FoamFile dictionary as {str: str} plus the key '_end'
    (offset of the first byte after the header). Stops on `format binary`."""
    m = _RE_FOAMFILE.search(raw, 0, min(len(raw), 8192))
    if m is None:
        return {"_end": 0}
    hdr = {}
    for e in _RE_ENTRY.finditer(m.group(1)):
        hdr[e.group(1).decode()] = e.group(2).decode(errors="replace").strip()
    hdr["_end"] = m.end()
    if require_ascii and hdr.get("format", "ascii").strip('"') != "ascii":
        raise FoamFormatError(
            "%s is written in OpenFOAM '%s' format; this reader needs ASCII.\n"
            "  Convert the case first: set `writeFormat ascii;` in system/controlDict and run\n"
            "  `foamFormatConvert` (add `-time <range>` / `-constant` as needed), then start again."
            % (path, hdr.get("format")))
    return hdr


def _numbers(buf, dtype=float):
    """Whitespace separated numbers in bytes/str -> 1D array."""
    try:
        a = np.fromstring(buf, sep=" ", dtype=dtype)
    except Exception:
        a = None
    if a is None:
        s = buf.decode() if isinstance(buf, (bytes, bytearray)) else buf
        a = np.array(s.split(), dtype=dtype)
    return a


_PAREN_TO_SPACE = bytes.maketrans(b"()", b"  ")
_RE_LIST_START = re.compile(rb"(?:^|\n)[ \t]*(\d+)\s*([({])")


def _list_start(raw, pos, path):
    """Find `N (` or `N {` at the start of a line after `pos`.
    Returns (n, kind, index of the first byte after the bracket)."""
    m = _RE_LIST_START.search(raw, pos)
    if m is None:
        raise FoamFormatError("%s: cannot find the list size / opening bracket" % path)
    return int(m.group(1)), m.group(2), m.end()


def _flat_list(raw, pos, path, dtype):
    """Parse a flat list `N ( v v v )` or a uniform list `N{v}` starting after
    `pos`. Returns (array, index after the closing bracket)."""
    n, kind, start = _list_start(raw, pos, path)
    if kind == b"{":
        end = raw.index(b"}", start)
        v = _numbers(raw[start:end], dtype)
        if v.size != 1:
            raise FoamFormatError("%s: cannot parse uniform list N{value}" % path)
        return np.full(n, v[0], dtype=dtype), end + 1
    end = raw.index(b")", start)
    a = _numbers(raw[start:end], dtype)
    if a.size != n:
        raise FoamFormatError("%s: list announces %d entries, found %d" % (path, n, a.size))
    return a, end + 1


def read_label_list(path):
    """labelList file (owner, neighbour, cellDecomposition, ...). Returns
    (int64 array, header dict)."""
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path)
    a, _ = _flat_list(raw, hdr["_end"], path, np.int64)
    return a, hdr


def read_cell_labels(path, n_cells=None):
    """One integer per cell from a labelList (constant/cellDecomposition as
    written by `decomposePar -cellDist`) or from a volScalarField (the
    `cellDist` field of the same tool)."""
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path)
    if b"internalField" in raw:
        a = parse_internal_field(path, is_vector=False, _raw=raw)
        if n_cells is not None:
            a = expand_uniform(a, n_cells, False)
        r = np.rint(a).astype(np.int64)
        if np.abs(a - r).max() > 1e-6:
            raise FoamFormatError("%s: values are not integers" % path)
        a = r
    else:
        a, _ = _flat_list(raw, hdr["_end"], path, np.int64)
    if n_cells is not None and a.size != n_cells:
        raise FoamFormatError("%s: %d entries, mesh has %d cells" % (path, a.size, n_cells))
    return a


# =============================================================================
# field files (logic of ref/cfdem_to_ip_files_v3.py: parse_internal_field)
# =============================================================================
_RE_COUNT = re.compile(rb"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(")
_RE_UNI_S = re.compile(rb"internalField\s+uniform\s+([\d.eE+\-]+)")
_RE_UNI_V = re.compile(rb"internalField\s+uniform\s+\(([^)]+)\)")


def parse_internal_field(path, is_vector, _raw=None):
    """Parse the internalField of an OpenFOAM ASCII field file.
    Returns (n,) or (n,3). A uniform field returns a length-1 (scalar) or
    length-3 (vector) array: use expand_uniform() to tile it to n cells."""
    raw = read_raw(path) if _raw is None else _raw
    parse_foam_header(raw, path)

    at = raw.find(b"internalField")
    if at < 0:
        raise ValueError("no internalField in %s" % path)
    bf_at = raw.find(b"boundaryField", at)
    if bf_at < 0:
        bf_at = len(raw)

    # uniform?
    if raw.find(b"nonuniform", at, bf_at) < 0:
        if is_vector:
            m = _RE_UNI_V.search(raw, at, bf_at)
            if not m:
                raise ValueError("cannot parse uniform vector in %s" % path)
            v = _numbers(m.group(1))
            if v.size != 3:
                raise ValueError("uniform vector in %s has %d components" % (path, v.size))
            return v
        m = _RE_UNI_S.search(raw, at, bf_at)
        if not m:
            raise ValueError("cannot parse uniform scalar in %s" % path)
        return np.array([float(m.group(1))])

    # nonuniform: count, then the data between the outer ( ... )
    m = _RE_COUNT.search(raw, at, bf_at)
    if not m:
        raise ValueError("cannot find nonuniform count in %s" % path)
    n_cells = int(m.group(1))
    open_paren = m.end() - 1
    close_paren = raw.rfind(b")", open_paren, bf_at)   # last ")" before boundaryField
    data = raw[open_paren + 1:close_paren]

    if is_vector:
        flat = _numbers(data.translate(_PAREN_TO_SPACE))
        if flat.size != 3 * n_cells:
            raise ValueError("%s: expected %d floats, got %d (scalar field read as vector?)"
                             % (path, 3 * n_cells, flat.size))
        return flat.reshape(-1, 3)
    if b"(" in data[:256]:
        raise ValueError("%s: vector field read as scalar" % path)
    arr = _numbers(data)
    if arr.size != n_cells:
        raise ValueError("%s: expected %d floats, got %d" % (path, n_cells, arr.size))
    return arr


def expand_uniform(field, n, is_vector):
    """Tile the result of parse_internal_field for a uniform field to n cells."""
    if is_vector:
        if field.shape == (3,):
            return np.tile(field, (n, 1))
    else:
        if field.size == 1 and n != 1:
            return np.full(n, float(field[0]))
    return field


def read_field(path, n_cells, is_vector):
    """Field of exactly n_cells values, (n,) or (n,3); uniform fields are tiled."""
    a = expand_uniform(parse_internal_field(path, is_vector), n_cells, is_vector)
    if len(a) != n_cells:
        raise ValueError("%s: %d values, expected %d cells" % (path, len(a), n_cells))
    return a


# =============================================================================
# time directories
# =============================================================================
def get_time_directories(case_path, t_start, t_end):
    """Sorted names of the time directories with t_start <= t <= t_end
    (float comparison with a small tolerance)."""
    hits = []
    for item in Path(case_path).iterdir():
        if not item.is_dir():
            continue
        try:
            t = float(item.name)
        except ValueError:
            continue
        if t_start - 1e-9 <= t <= t_end + 1e-9:
            hits.append((t, item.name))
    hits.sort(key=lambda x: x[0])
    return [name for _, name in hits]


def time_name(t, digits=6):
    """Name OpenFOAM gives a time directory: 2 -> '2', 2.005 -> '2.005'."""
    s = ("%.*f" % (digits, t)).rstrip("0").rstrip(".")
    if s in ("", "-0"):
        s = "0"
    return s


# =============================================================================
# polyMesh
# =============================================================================
def _read_points(path):
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path)
    n, kind, start = _list_start(raw, hdr["_end"], path)
    if kind != b"(":
        raise FoamFormatError("%s: unexpected uniform point list" % path)
    end = raw.rfind(b")")
    a = _numbers(raw[start:end].translate(_PAREN_TO_SPACE))
    if a.size != 3 * n:
        raise FoamFormatError("%s: %d points announced, %d numbers found" % (path, n, a.size))
    return a.reshape(n, 3)


def _read_faces(path):
    """faces file -> (ptr[nF+1], nodes[...]) in CSR form."""
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path)
    cls = hdr.get("class", "faceList")

    if cls == "faceCompactList":
        # two flat lists: offsets (nF+1), then all labels
        ptr, pos = _flat_list(raw, hdr["_end"], path, np.int64)
        nodes, _ = _flat_list(raw, pos, path, np.int64)
        if ptr[-1] != nodes.size:
            raise FoamFormatError("%s: faceCompactList offsets and labels do not fit" % path)
        return ptr, nodes

    n, kind, start = _list_start(raw, hdr["_end"], path)
    if kind != b"(":
        raise FoamFormatError("%s: unexpected uniform face list" % path)
    end = raw.rfind(b")")
    body = raw[start:end]

    # fast path: every face a quad  ->  "4(a b c d)" x n  ->  5 n integers
    flat = _numbers(body.translate(_PAREN_TO_SPACE), np.int64)
    if flat.size == 5 * n and n > 0:
        flat = flat.reshape(n, 5)
        if (flat[:, 0] == 4).all():          # exact: the first non-quad would show here
            return np.arange(n + 1, dtype=np.int64) * 4, np.ascontiguousarray(flat[:, 1:]).ravel()
        flat = flat.ravel()

    # general: mark the end of every face with -1, the count is the entry after the previous end
    flat = _numbers(body.replace(b"(", b" ").replace(b")", b" -1 "), np.int64)
    ends = np.flatnonzero(flat == -1)
    if ends.size != n:
        raise FoamFormatError("%s: %d faces announced, %d found" % (path, n, ends.size))
    starts = np.empty(n, dtype=np.int64)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    counts = flat[starts]
    if not (counts == ends - starts - 1).all():
        bad = int(np.flatnonzero(counts != ends - starts - 1)[0])
        raise FoamFormatError("%s: face %d: size does not match its label list" % (path, bad))
    keep = np.ones(flat.size, dtype=bool)
    keep[starts] = False
    keep[ends] = False
    ptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=ptr[1:])
    return ptr, flat[keep]


_RE_PATCH_KEY = {
    "type": re.compile(rb"\btype\s+([^;\s]+)\s*;"),
    "nFaces": re.compile(rb"\bnFaces\s+(\d+)\s*;"),
    "startFace": re.compile(rb"\bstartFace\s+(\d+)\s*;"),
}


def _read_boundary(path):
    """polyMesh/boundary -> list of dicts {name, type, nFaces, startFace}."""
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path)
    body = re.sub(rb"//[^\n]*", b"", raw[hdr["_end"]:])
    body = re.sub(rb"/\*.*?\*/", b"", body, flags=re.S)
    m = re.search(rb"(\d+)\s*\(", body)
    if m is None:
        raise FoamFormatError("%s: cannot find the patch list" % path)
    n = int(m.group(1))
    pos = m.end()
    patches = []
    while len(patches) < n:
        mm = re.compile(rb"\s*([^\s{}();]+)\s*\{").match(body, pos)
        if mm is None:
            raise FoamFormatError("%s: cannot parse patch %d" % (path, len(patches)))
        depth, i = 1, mm.end()
        while depth > 0:
            j_open, j_close = body.find(b"{", i), body.find(b"}", i)
            if j_close < 0:
                raise FoamFormatError("%s: unbalanced braces" % path)
            if 0 <= j_open < j_close:
                depth, i = depth + 1, j_open + 1
            else:
                depth, i = depth - 1, j_close + 1
        block = body[mm.end():i - 1]
        p = {"name": mm.group(1).decode()}
        for key, rx in _RE_PATCH_KEY.items():
            k = rx.search(block)
            if k is None:
                raise FoamFormatError("%s: patch %s has no '%s'" % (path, p["name"], key))
            p[key] = k.group(1).decode() if key == "type" else int(k.group(1))
        patches.append(p)
        pos = i
    return patches


def _sum_by(index, values, n):
    """Sum rows of `values` (m,) or (m,k) into n bins."""
    if values.ndim == 1:
        return np.bincount(index, weights=values, minlength=n)
    out = np.empty((n, values.shape[1]))
    for k in range(values.shape[1]):
        out[:, k] = np.bincount(index, weights=values[:, k], minlength=n)
    return out


class PolyMesh:
    """OpenFOAM polyMesh with geometry.

    points (nP,3); face_ptr (nF+1), face_nodes: nodes of face f are
    face_nodes[face_ptr[f]:face_ptr[f+1]], right-hand normal pointing out of
    the owner cell; owner (nF); neighbour (nInternal); patches (list of dicts).
    """

    def __init__(self, points, face_ptr, face_nodes, owner, neighbour, patches, n_cells=None):
        self.points = np.ascontiguousarray(points, dtype=float)
        self.face_ptr = np.asarray(face_ptr, dtype=np.int64)
        self.face_nodes = np.asarray(face_nodes, dtype=np.int64)
        self.owner = np.asarray(owner, dtype=np.int64)
        self.neighbour = np.asarray(neighbour, dtype=np.int64)
        self.patches = patches
        self.n_points = len(self.points)
        self.n_faces = len(self.face_ptr) - 1
        self.n_internal = len(self.neighbour)
        nc = int(self.owner.max()) + 1 if self.owner.size else 0
        if self.neighbour.size:
            nc = max(nc, int(self.neighbour.max()) + 1)
        if n_cells is not None:
            if nc > n_cells:
                raise FoamFormatError("owner/neighbour refer to cell %d, but nCells is %d" % (nc - 1, n_cells))
            nc = n_cells
        self.n_cells = nc
        self._geom = None
        self._cf = None
        self._hex = None
        self._check()

    # ------------------------------------------------------------------ checks
    def _check(self):
        if len(self.owner) != self.n_faces:
            raise FoamFormatError("owner has %d entries, faces has %d" % (len(self.owner), self.n_faces))
        if self.n_internal > self.n_faces:
            raise FoamFormatError("neighbour is longer than owner")
        if self.face_nodes.size and (self.face_nodes.min() < 0 or self.face_nodes.max() >= self.n_points):
            raise FoamFormatError("faces refer to points outside 0..%d" % (self.n_points - 1))
        if self.owner.size and self.owner.min() < 0:
            raise FoamFormatError("negative owner label")
        if self.neighbour.size and self.neighbour.min() < 0:
            raise FoamFormatError("negative neighbour label (old-style neighbour file with -1 entries "
                                  "is not supported; internal faces must come first)")
        if (np.diff(self.face_ptr) < 3).any():
            raise FoamFormatError("a face has fewer than 3 points")
        for p in self.patches:
            if p["startFace"] < self.n_internal or p["startFace"] + p["nFaces"] > self.n_faces:
                raise FoamFormatError("patch %s: faces %d..%d are outside the boundary face range %d..%d"
                                      % (p["name"], p["startFace"], p["startFace"] + p["nFaces"] - 1,
                                         self.n_internal, self.n_faces - 1))

    # ---------------------------------------------------------------- geometry
    def _face_geometry(self, chunk=250000):
        nF = self.n_faces
        fc = np.empty((nF, 3))
        fa = np.empty((nF, 3))
        counts = np.diff(self.face_ptr)
        pts = self.points
        for f0 in range(0, nF, chunk):
            f1 = min(nF, f0 + chunk)
            s, e = int(self.face_ptr[f0]), int(self.face_ptr[f1])
            cnt = counts[f0:f1]
            P = pts[self.face_nodes[s:e]]                       # (m,3) points of all faces of the chunk
            if (cnt == cnt[0]).all():
                # all faces of the chunk have the same number of points (the usual all-quad mesh)
                k = int(cnt[0])
                P = P.reshape(f1 - f0, k, 3)
                est = P.mean(axis=1)
                Pn = np.roll(P, -1, axis=1)
                E = est[:, None, :]
                n2 = np.cross(Pn - P, E - P)                    # twice the area vector of every triangle
                a = np.sqrt((n2 * n2).sum(axis=2))
                sum_n = n2.sum(axis=1)
                sum_a = a.sum(axis=1)
                sum_ac = (a[:, :, None] * (P + Pn + E)).sum(axis=1)
            else:
                lp = self.face_ptr[f0:f1 + 1] - s
                est = np.add.reduceat(P, lp[:-1], axis=0) / cnt[:, None]
                nxt = np.arange(1, e - s + 1)
                nxt[lp[1:] - 1] = lp[:-1]                       # last point of a face -> its first point
                Pn = P[nxt]
                E = np.repeat(est, cnt, axis=0)
                n2 = np.cross(Pn - P, E - P)
                a = np.sqrt((n2 * n2).sum(axis=1))
                sum_n = np.add.reduceat(n2, lp[:-1], axis=0)
                sum_a = np.add.reduceat(a, lp[:-1])
                sum_ac = np.add.reduceat(a[:, None] * (P + Pn + E), lp[:-1], axis=0)
            ok = sum_a > 1e-300
            c = est.copy()
            c[ok] = sum_ac[ok] / (3.0 * sum_a[ok, None])
            fc[f0:f1] = c
            fa[f0:f1] = 0.5 * sum_n
        return fc, fa

    def _cell_geometry(self, fc, fa):
        nC, nI = self.n_cells, self.n_internal
        own, nei = self.owner, self.neighbour
        nfc = np.bincount(own, minlength=nC) + np.bincount(nei, minlength=nC)
        if (nfc < 4).any():
            raise FoamFormatError("cell %d has only %d faces" % (int(np.argmin(nfc)), int(nfc.min())))
        est = (_sum_by(own, fc, nC) + _sum_by(nei, fc[:nI], nC)) / nfc[:, None]
        # pyramids face -> estimated centre; 3 x pyramid volume = A . (fc - est), sign by side
        d = fc - est[own]
        p3 = (fa * d).sum(axis=1)
        vol3 = np.bincount(own, weights=p3, minlength=nC)
        ctr = _sum_by(own, p3[:, None] * (0.75 * fc + 0.25 * est[own]), nC)
        d = fc[:nI] - est[nei]
        p3 = -(fa[:nI] * d).sum(axis=1)
        vol3 += np.bincount(nei, weights=p3, minlength=nC)
        ctr += _sum_by(nei, p3[:, None] * (0.75 * fc[:nI] + 0.25 * est[nei]), nC)
        ok = np.abs(vol3) > 1e-300
        cc = est.copy()
        cc[ok] = ctr[ok] / vol3[ok, None]
        return cc, vol3 / 3.0, nfc

    def geometry(self):
        """(face_centres, face_areas, cell_centres, cell_volumes)."""
        if self._geom is None:
            fc, fa = self._face_geometry()
            cc, vol, nfc = self._cell_geometry(fc, fa)
            self._geom = (fc, fa, cc, vol)
            self._nfc = nfc
        return self._geom

    @property
    def face_centres(self):
        return self.geometry()[0]

    @property
    def face_areas(self):
        return self.geometry()[1]

    @property
    def cell_centres(self):
        return self.geometry()[2]

    @property
    def cell_volumes(self):
        return self.geometry()[3]

    def closedness(self):
        """Per cell |sum of outward face area vectors| / (sum of |face area|).
        All faces of the cell enter, boundary faces included."""
        fa = self.face_areas
        nC, nI = self.n_cells, self.n_internal
        s = _sum_by(self.owner, fa, nC) - _sum_by(self.neighbour, fa[:nI], nC)
        mag = np.sqrt((fa * fa).sum(axis=1))
        surf = np.bincount(self.owner, weights=mag, minlength=nC) + \
            np.bincount(self.neighbour, weights=mag[:nI], minlength=nC)
        return np.sqrt((s * s).sum(axis=1)) / np.maximum(surf, 1e-300)

    def boundary_volume(self):
        """Volume enclosed by the boundary faces, (1/3) sum fc . A (divergence
        theorem). Independent of the cells: the sum of cell volumes must match."""
        fc, fa = self.geometry()[:2]
        nI = self.n_internal
        return float((fc[nI:] * fa[nI:]).sum() / 3.0)

    def bounding_box(self):
        return self.points.min(axis=0), self.points.max(axis=0)

    # ----------------------------------------------------- cell -> faces (CSR)
    def cell_faces(self):
        """(ptr[nC+1], face[...], flip[...]): faces of every cell; flip is True
        where the cell is the neighbour of the face (face normal points INTO
        the cell). Owner faces first, each group in ascending face order."""
        if self._cf is None:
            nC, nI, nF = self.n_cells, self.n_internal, self.n_faces
            cells = np.concatenate([self.owner, self.neighbour])
            faces = np.concatenate([np.arange(nF, dtype=np.int64), np.arange(nI, dtype=np.int64)])
            flip = np.zeros(nF + nI, dtype=bool)
            flip[nF:] = True
            order = np.argsort(cells, kind="stable")
            ptr = np.zeros(nC + 1, dtype=np.int64)
            np.cumsum(np.bincount(cells, minlength=nC), out=ptr[1:])
            self._cf = (ptr, faces[order], flip[order])
        return self._cf

    # --------------------------------------------------------------- hexahedra
    def hex_nodes(self, chunk=100000):
        """(is_hex[nC], nodes[nC,8]): VTK_HEXAHEDRON node order for every cell
        that is a topological hexahedron (6 quads, 8 points, every corner
        with one edge to the opposite face); nodes are -1 for other cells.

        VTK order: points 0-3 are one face with the right-hand normal pointing
        into the cell, 4-7 the opposite face, point i+4 joined to point i."""
        if self._hex is not None:
            return self._hex
        nC = self.n_cells
        ptr, cface, cflip = self.cell_faces()
        fcount = np.diff(self.face_ptr)
        nfc = np.diff(ptr)
        quad_faces = np.bincount(np.repeat(np.arange(nC), nfc), weights=(fcount[cface] == 4), minlength=nC)
        cand = np.flatnonzero((nfc == 6) & (quad_faces == 6))
        is_hex = np.zeros(nC, dtype=bool)
        nodes = np.full((nC, 8), -1, dtype=np.int64)
        for c0 in range(0, cand.size, chunk):
            cells = cand[c0:c0 + chunk]
            m = cells.size
            idx = ptr[cells][:, None] + np.arange(6)[None, :]
            F = cface[idx]                                     # (m,6) faces
            flip0 = cflip[idx[:, 0]]
            Q = self.face_nodes[self.face_ptr[F][:, :, None] + np.arange(4)[None, None, :]]   # (m,6,4)
            B = Q[:, 0, :].copy()
            # owner side: normal points out of the cell -> reverse to make it point inwards
            rev = ~flip0
            B[rev] = B[rev][:, [0, 3, 2, 1]]
            S = Q[:, 1:, :]                                    # (m,5,4) the other faces
            u = np.concatenate([S, np.roll(S, -1, axis=2)], axis=2).reshape(m, 40)
            v = np.concatenate([np.roll(S, -1, axis=2), S], axis=2).reshape(m, 40)
            v_in_bottom = (v[:, :, None] == B[:, None, :]).any(axis=2)
            T = np.empty((m, 4), dtype=np.int64)
            good = np.ones(m, dtype=bool)
            for i in range(4):
                hit = (u == B[:, i, None]) & ~v_in_bottom
                nh = hit.sum(axis=1)
                vmax = np.where(hit, v, -1).max(axis=1)
                vmin = np.where(hit, v, np.iinfo(np.int64).max).min(axis=1)
                good &= (nh == 2) & (vmax == vmin)
                T[:, i] = vmax
            H = np.concatenate([B, T], axis=1)
            Hs = np.sort(H, axis=1)
            good &= (np.diff(Hs, axis=1) > 0).all(axis=1)       # 8 different points
            # all points of the cell are these 8
            allp = np.sort(Q.reshape(m, 24), axis=1)
            n_unique = 1 + (np.diff(allp, axis=1) > 0).sum(axis=1)
            good &= (n_unique == 8)
            is_hex[cells[good]] = True
            nodes[cells[good]] = H[good]
        self._hex = (is_hex, nodes)
        return self._hex

    # -------------------------------------------------------------- polyhedra
    def poly_stream(self, cells):
        """Face streams of the given cells for VTK_POLYHEDRON (type 42).

        Returns (n_faces[len(cells)], face_sizes[sum n_faces], nodes[sum sizes]):
        for every cell its faces (owner faces first), every face with its
        points ordered so that the right-hand normal points OUT of the cell."""
        cells = np.asarray(cells, dtype=np.int64)
        ptr, cface, cflip = self.cell_faces()
        nf = ptr[cells + 1] - ptr[cells]
        tot = int(nf.sum())
        if tot == 0:
            return nf, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
        first = np.zeros(cells.size, dtype=np.int64)
        np.cumsum(nf[:-1], out=first[1:])
        e = np.repeat(ptr[cells] - first, nf) + np.arange(tot)   # entries of the cell-face CSR
        f = cface[e]
        flip = cflip[e]
        size = self.face_ptr[f + 1] - self.face_ptr[f]
        start = np.zeros(tot, dtype=np.int64)
        np.cumsum(size[:-1], out=start[1:])
        ntot = int(size.sum())
        j = np.arange(ntot) - np.repeat(start, size)             # position inside the face
        sz = np.repeat(size, size)
        j = np.where(np.repeat(flip, size), (sz - j) % sz, j)    # reversed: 0, n-1, n-2, ..., 1
        return nf, size, self.face_nodes[np.repeat(self.face_ptr[f], size) + j]

    def summary(self):
        counts = np.diff(self.face_ptr)
        return ("%d points, %d faces (%d internal, %d boundary; %d..%d points per face), %d cells, %d patches"
                % (self.n_points, self.n_faces, self.n_internal, self.n_faces - self.n_internal,
                   int(counts.min()), int(counts.max()), self.n_cells, len(self.patches)))


_RE_NOTE = re.compile(r"nCells:\s*(\d+)")


def read_polymesh(case, region=None):
    """Read <case>/constant/[region/]polyMesh. `case` may also be the polyMesh
    directory itself."""
    case = Path(case)
    d = case
    if not (find_file(case / "points") and find_file(case / "faces")):
        d = case / "constant"
        if region:
            d = d / region
        d = d / "polyMesh"
    for name in ("points", "faces", "owner", "neighbour"):
        if find_file(d / name) is None:
            raise FileNotFoundError("%s: no '%s' (or '%s.gz') - is this an OpenFOAM case with a "
                                    "reconstructed (not decomposed) mesh?" % (d, name, name))
    points = _read_points(d / "points")
    ptr, nodes = _read_faces(d / "faces")
    owner, hdr = read_label_list(d / "owner")
    neighbour, _ = read_label_list(d / "neighbour")
    n_cells = None
    m = _RE_NOTE.search(hdr.get("note", ""))
    if m:
        n_cells = int(m.group(1))
    if neighbour.size == owner.size and neighbour.size and neighbour.min() < 0:
        raise FoamFormatError("%s/neighbour has -1 entries (very old mesh format); run "
                              "`foamFormatConvert` or `renumberMesh -overwrite` with a current OpenFOAM" % d)
    patches = _read_boundary(d / "boundary") if find_file(d / "boundary") else []
    return PolyMesh(points, ptr, nodes, owner, neighbour, patches, n_cells=n_cells)


# =============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: of_reader.py CASE     (prints a mesh summary with geometry checks)")
    msh = read_polymesh(sys.argv[1])
    print(msh.summary())
    for p in msh.patches:
        print("  patch %-20s type %-12s nFaces %8d startFace %9d" % (p["name"], p["type"], p["nFaces"], p["startFace"]))
    vol = msh.cell_volumes
    lo, hi = msh.bounding_box()
    is_hex, _ = msh.hex_nodes()
    print("hexahedra: %d of %d cells" % (int(is_hex.sum()), msh.n_cells))
    print("sum of cell volumes   %.15e" % vol.sum())
    print("volume inside boundary %.15e (divergence theorem), relative difference %.2e"
          % (msh.boundary_volume(), abs(vol.sum() - msh.boundary_volume()) / abs(vol.sum())))
    print("bounding box volume   %.15e   [%g %g %g] .. [%g %g %g]" % ((np.prod(hi - lo),) + tuple(lo) + tuple(hi)))
    print("cell volume min %.6e max %.6e" % (vol.min(), vol.max()))
    print("closedness max |sum A| / sum |A| = %.3e" % msh.closedness().max())

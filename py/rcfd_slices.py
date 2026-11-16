#!/usr/bin/env python3
"""Slices through an rCFD result, straight from the snapshots: no ParaView, no OpenFOAM.

    cd <case>                                   # the folder with data/ and out/
    python3 <rcfd_foam2>/py/rcfd_slices.py --mark 0 10

writes into ./slices/:
    slice_x0_t<time>.png, slice_y0_t<time>.png, slice_z<pos>_t<time>.png
                              all species on one plane at one time
    <species>_x0_times.png ...  one species on one plane at all times, one colour scale
    cell_centres_f32.npy      positions of all cells (written once, read again afterwards)

Planes: --x, --y, --z take one or more positions in mm (defaults: x 0, y 0, z middle of the
domain). The flow direction z is drawn from left to right. --mark z1 z2 draws two lines at
these z positions (the reacting part of the lattice case: --mark 0 10).

Fields: by default the concentration fields c0, c1, ... of the first phase, labelled A, B, C ...
(--labels to change); --extra adds every other field of the snapshot (for the lattice case
sel_CD, share_E, wall_rate). --times 1 2.5 5 picks snapshots near these times.

Many snapshots: --every 5 takes every fifth one; the <species>_times pictures show --montage
times (default 6, evenly spread); --gif also writes an animation slice_<plane>.gif (needs Pillow).

How a slice is made: a cell belongs to a plane if the plane passes through it (distance of the
cell centre below half the cell size, cell size = cube root of the volume). Each such cell is
painted onto a pixel grid (--pixel, default 0.02 mm), large cells first. Pixels that no cell
covers are solid and drawn grey.

Needs numpy and matplotlib. Reads data/topo/ (cell positions) and out/snap_*.bin; if out/ is
gone, the binary fields in rCFD_results/<time>/ are used. Memory: about 0.3 GB per million
cells is enough; the first run on 16 million cells reads the positions for a few minutes.
"""
import argparse, glob, itertools, os, re, sys, time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

SOLID = "#d0d0d0"
PALE = "#f2f2f2"
NICE = {"sel_CD": "Selectivity (C + D) / (C + D + E)", "share_E": "Share of the reacting B that becomes E",
        "wall_rate": "A converted at the wall"}
UNIT = {"sel_CD": "share", "share_E": "share", "wall_rate": "1/s"}


# ----------------------------------------------------------------------------- input
def read_centres(case, cache):
    """-> array [n_cells, 4] float32: x, y, z [m] and cell size (cube root of the volume) [m],
    in the cell numbering of the snapshots (the undecomposed mesh)."""
    if cache and os.path.isfile(cache):
        c = np.load(cache)
        print("cell positions from %s: %d cells" % (cache, len(c)), flush=True)
        return c
    meta = open(os.path.join(case, "data/topo/meta.txt")).read().split()
    n_global, n_ranks = int(meta[0]), int(meta[1])
    c = np.full((n_global, 4), np.nan, dtype=np.float32)
    t0 = time.time()
    for r in range(n_ranks):
        gid = np.fromfile(os.path.join(case, "data/topo/gid_%d.bin" % r), dtype="<i4")
        with open(os.path.join(case, "data/topo/topo_%d" % r)) as f:
            n_cells = int(f.readline().split()[0])
            block = "".join(itertools.islice(f, n_cells))
        a = np.array(block.split(), dtype=np.float64).reshape(n_cells, -1)
        own = a[:, 4].astype(int) == r
        g = gid[:n_cells][own]
        c[g, :3] = a[own, :3]
        c[g, 3] = np.cbrt(a[own, 3])
        if r % 8 == 0 or r == n_ranks - 1:
            print("  cell positions: processor %d of %d (%.0f s)" % (r + 1, n_ranks, time.time() - t0), flush=True)
    if np.isnan(c[:, 0]).any():
        sys.exit("data/topo: %d cells without a position" % int(np.isnan(c[:, 0]).sum()))
    if cache:
        np.save(cache, c)
        print("written %s (%.0f MB)" % (cache, os.path.getsize(cache) / 1e6), flush=True)
    return c


class Snapshots:
    """Fields of all cells at the snapshot times, read on demand for a list of cells."""

    def __init__(self, case):
        self.case = case
        self.snap = sorted(glob.glob(os.path.join(case, "out", "snap_*.bin")))
        self.times, self.maps, self.names = [], [], []
        if self.snap:
            f = os.path.join(case, "out", "fields.txt")
            names = [l.strip() for l in open(f)] if os.path.isfile(f) else []
            for p in self.snap:
                head = np.fromfile(p, dtype=np.uint8, count=24)
                n = int(head[:8].view("<i8")[0]); nf = int(head[8:12].view("<i4")[0])
                size = os.path.getsize(p)
                t_off, d_off = (12, 20) if size == 20 + 4 * n * nf else (16, 24)
                self.times.append(float(head[t_off:t_off + 8].view("<f8")[0]))
                self.maps.append(np.memmap(p, dtype="<f4", mode="r", offset=d_off, shape=(nf, n)))
            self.names = (names + ["field%d" % i for i in range(len(names), nf)])[:nf]
            self.source = "out/snap_*.bin"
        else:                                           # binary fields written by snap_to_foam.py --binary
            dirs = [d for d in glob.glob(os.path.join(case, "rCFD_results", "*")) if re.fullmatch(r"[0-9.eE+-]+", os.path.basename(d))]
            dirs.sort(key=lambda d: float(os.path.basename(d)))
            if not dirs:
                sys.exit("no snapshots: neither out/snap_*.bin nor rCFD_results/<time>/ in %s" % case)
            self.names = sorted(os.path.basename(f)[5:] for f in glob.glob(os.path.join(dirs[-1], "rCFD_*")))
            for d in dirs:
                self.times.append(float(os.path.basename(d)))
                self.maps.append(d)
            self.source = "rCFD_results/<time>/"

    def values(self, i_time, name, cells):
        m = self.maps[i_time]
        if not isinstance(m, str):
            return np.asarray(m[self.names.index(name)][cells], dtype=np.float64)
        path = os.path.join(m, "rCFD_" + name)
        with open(path, "rb") as f:
            head = f.read(4096)
        if b"binary" not in head:
            sys.exit("%s is a text file: run snap_to_foam.py with --binary, or keep out/" % path)
        k = head.index(b"nonuniform List<scalar>")
        mm = re.search(rb"\s*(\d+)\s*\(", head[k + 23:])
        n = int(mm.group(1)); off = k + 23 + mm.end()
        return np.asarray(np.memmap(path, dtype="<f8", mode="r", offset=off, shape=(n,))[cells], dtype=np.float64)


# ----------------------------------------------------------------------------- slice
class Plane:
    """Cells cut by the plane normal to an axis at a position, and the pixel map of the cut."""

    def __init__(self, centres, lo, hi, axis, pos, pixel):
        self.axis, self.pos = axis, pos
        u_ax, v_ax = {0: (2, 1), 1: (2, 0), 2: (0, 1)}[axis]          # horizontal, vertical
        self.u_name, self.v_name = "xyz"[u_ax], "xyz"[v_ax]
        cut = np.abs(centres[:, axis] - pos) <= 0.5 * centres[:, 3]
        self.cells = np.nonzero(cut)[0]
        c = centres[self.cells]
        self.u0, self.u1, self.v0, self.v1 = lo[u_ax], hi[u_ax], lo[v_ax], hi[v_ax]
        nu = int(np.ceil((self.u1 - self.u0) / pixel)); nv = int(np.ceil((self.v1 - self.v0) / pixel))
        owner = np.full((nv, nu), -1, dtype=np.int64)
        iu = np.floor((c[:, u_ax] - self.u0) / pixel).astype(np.int64)
        iv = np.floor((c[:, v_ax] - self.v0) / pixel).astype(np.int64)
        k = np.maximum(1, np.ceil(1.15 * c[:, 3] / pixel).astype(np.int64))      # footprint in pixels
        for kk in np.unique(k)[::-1]:                                          # large cells first
            sel = np.nonzero(k == kk)[0]
            for dv in range(kk):
                jv = iv[sel] - kk // 2 + dv
                okv = (jv >= 0) & (jv < nv)
                for du in range(kk):
                    ju = iu[sel] - kk // 2 + du
                    ok = okv & (ju >= 0) & (ju < nu)
                    owner[jv[ok], ju[ok]] = sel[ok]
        self.owner = owner
        self.solid = owner < 0

    def image(self, values):
        return np.ma.array(values[np.maximum(self.owner, 0)].astype(np.float32), mask=self.solid)

    def extent_mm(self):
        return [1e3 * self.u0, 1e3 * self.u1, 1e3 * self.v0, 1e3 * self.v1]

    def tag(self):
        return "%s%g" % ("xyz"[self.axis], round(1e3 * self.pos, 3))

    def title(self):
        return "plane %s = %g mm" % ("xyz"[self.axis], round(1e3 * self.pos, 3))


def nice_max(v):
    if not np.isfinite(v) or v <= 0:
        return 1.0
    e = 10.0 ** np.floor(np.log10(v))
    for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if v <= m * e * 1.0001:
            return m * e
    return 10 * e


def field_norm(name, stack):
    """One colour scale per field for all times of a plane."""
    data = np.concatenate([s.compressed() for s in stack])
    if name == "wall_rate":
        pos = data[data > 0]
        hi = nice_max(np.percentile(pos, 99.5)) if pos.size else 1.0
        return LogNorm(vmin=hi * 1e-3, vmax=hi), "inferno"
    if name == "sel_CD":                             # zero = no product yet: drawn pale
        pos = data[data > 0]
        if not pos.size:
            return Normalize(0.0, 1.0), "viridis"
        lo, hi = np.floor(20 * np.percentile(pos, 1)) / 20, np.ceil(20 * np.percentile(pos, 99)) / 20
        return Normalize(lo, max(hi, lo + 0.05)), "viridis"
    if name == "share_E":
        return Normalize(0.0, nice_max(np.percentile(data, 99.9)) if data.size else 1.0), "magma"
    hi = nice_max(np.percentile(data, 99.9)) if data.size else 1.0
    return Normalize(0.0, min(hi, 1.0) if hi <= 1.0001 else hi), "viridis"


def draw(ax, plane, img, norm, cmap, mark, pale_zero=False):
    cm = matplotlib.colormaps[cmap].copy() if hasattr(matplotlib, "colormaps") else plt.get_cmap(cmap).copy()
    cm.set_bad(SOLID)
    if pale_zero:                                   # zero (nothing there yet) is drawn pale, not as the lowest colour
        cm.set_under(PALE)
        img = np.ma.array(np.where(img.filled(0.0) > 0.0, img.filled(0.0), norm.vmin - 1.0 if not isinstance(norm, LogNorm) else 0.5 * norm.vmin), mask=plane.solid)
    im = ax.imshow(img, origin="lower", extent=plane.extent_mm(), cmap=cm, norm=norm, interpolation="nearest", aspect="equal")
    if mark and plane.axis != 2:
        for z in mark:
            ax.axvline(z, color="white", lw=0.9, ls=(0, (4, 3)))
            ax.axvline(z, color="black", lw=0.35, ls=(0, (4, 3)))
    ax.tick_params(labelsize=7, length=2.5, width=0.5)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    return im


def colourbar(fig, im, ax, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.012, aspect=18)
    cb.ax.tick_params(labelsize=7, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.set_label(label, fontsize=7.5)
    return cb


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default=".", help="case folder (default: here)")
    p.add_argument("--out", default=None, help="folder for the pictures (default: <case>/slices)")
    p.add_argument("--x", type=float, nargs="*", default=[0.0], help="positions of the planes normal to x [mm]")
    p.add_argument("--y", type=float, nargs="*", default=[0.0], help="positions of the planes normal to y [mm]")
    p.add_argument("--z", type=float, nargs="*", default=None, help="positions of the planes normal to z [mm]")
    p.add_argument("--times", type=float, nargs="*", default=None, help="snapshot times [s] (default: all)")
    p.add_argument("--every", type=int, default=1, help="only every N-th snapshot")
    p.add_argument("--montage", type=int, default=6, help="number of times in the <species>_times pictures (evenly spread)")
    p.add_argument("--gif", action="store_true", help="also an animation slice_<plane>.gif of the per-time pictures (needs Pillow)")
    p.add_argument("--labels", nargs="*", default=None, help="names of the concentration fields (default A B C ...)")
    p.add_argument("--extra", action="store_true", help="also the other fields of the snapshot")
    p.add_argument("--mark", type=float, nargs=2, default=None, help="two z positions [mm] marked by lines")
    p.add_argument("--pixel", type=float, default=0.02, help="pixel size [mm]")
    p.add_argument("--dpi", type=int, default=220)
    p.add_argument("--no-cache", action="store_true", help="do not write cell_centres_f32.npy")
    a = p.parse_args()

    out = a.out or os.path.join(a.case, "slices")
    os.makedirs(out, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                         "axes.labelsize": 8, "figure.dpi": 100, "savefig.dpi": a.dpi, "savefig.bbox": "tight",
                         "savefig.pad_inches": 0.06})

    snaps = Snapshots(a.case)
    centres = read_centres(a.case, None if a.no_cache else os.path.join(out, "cell_centres_f32.npy"))
    half = 0.5 * centres[:, 3]
    lo = (centres[:, :3] - half[:, None]).min(0).astype(float); hi = (centres[:, :3] + half[:, None]).max(0).astype(float)
    print("domain [mm]: x %.2f..%.2f, y %.2f..%.2f, z %.2f..%.2f; %d snapshots from %s, t = %g .. %g s"
          % (1e3 * lo[0], 1e3 * hi[0], 1e3 * lo[1], 1e3 * hi[1], 1e3 * lo[2], 1e3 * hi[2], len(snaps.times), snaps.source,
             snaps.times[0], snaps.times[-1]), flush=True)

    conc = [n for n in snaps.names if re.fullmatch(r"c\d+\..*", n)]
    phase = conc[0].split(".", 1)[1] if conc else ""
    conc = [n for n in conc if n.endswith("." + phase)]
    labels = a.labels or [chr(ord("A") + i) for i in range(len(conc))]
    if len(labels) != len(conc):
        sys.exit("%d labels for %d concentration fields (%s)" % (len(labels), len(conc), " ".join(conc)))
    fields = list(zip(conc, labels))
    if a.extra:
        fields += [(n, n) for n in snaps.names if n not in conc and not n.startswith("vof.") and not re.fullmatch(r"c\d+\..*", n)]
    if not fields:
        sys.exit("no fields to draw; the snapshots hold: %s" % " ".join(snaps.names))

    it = list(range(len(snaps.times)))[::max(1, a.every)]
    if a.times:
        it = sorted(set(int(np.argmin(np.abs(np.array(snaps.times) - t))) for t in a.times))
    # times of the montage pictures; the colour scales are taken from these states
    im_idx = sorted(set(int(round(v)) for v in np.linspace(0, len(it) - 1, max(2, min(a.montage, len(it))))))
    it_m = [it[k] for k in im_idx]
    zpos = a.z if a.z is not None else [round(0.5e3 * (lo[2] + hi[2]), 3)]
    planes = [(0, 1e-3 * v) for v in a.x] + [(1, 1e-3 * v) for v in a.y] + [(2, 1e-3 * v) for v in zpos]

    for axis, pos in planes:
        pl = Plane(centres, lo, hi, axis, pos, 1e-3 * a.pixel)
        if pl.cells.size == 0:
            print("%s: no cells" % pl.title()); continue
        print("%s: %d cells, %d x %d pixels" % (pl.title(), pl.cells.size, pl.owner.shape[1], pl.owner.shape[0]), flush=True)
        imgs = {n: [pl.image(snaps.values(i, n, pl.cells)) for i in it_m] for n, _ in fields}
        norms = {n: field_norm(n, imgs[n]) for n, _ in fields}
        frames = []
        ext = pl.extent_mm(); w_mm, h_mm = ext[1] - ext[0], ext[3] - ext[2]
        wide = w_mm > 1.5 * h_mm                       # planes along the flow: panels one below the other

        def panel_size(n_panels):
            if wide:
                w = 7.6; h = w * h_mm / w_mm
                return 1, n_panels, (w, n_panels * (h + 0.42) + 0.6)
            ncol = min(n_panels, 4); nrow = -(-n_panels // ncol)
            w = 2.5; h = w * h_mm / w_mm
            return ncol, nrow, (ncol * (w + 0.85) + 0.2, nrow * (h + 0.55) + 0.6)

        def label_axes(axs, n_used, ncol):
            for k, ax in enumerate(axs):
                if k >= n_used:
                    ax.set_axis_off(); continue
                if wide:
                    ax.set_ylabel("%s in mm" % pl.v_name)
                    if k == n_used - 1: ax.set_xlabel("%s in mm (flow from left to right)" % pl.u_name)
                    else: ax.tick_params(labelbottom=False)
                else:
                    if k + ncol >= n_used: ax.set_xlabel("%s in mm" % pl.u_name)
                    else: ax.tick_params(labelbottom=False)
                    if k % ncol == 0: ax.set_ylabel("%s in mm" % pl.v_name)
                    else: ax.tick_params(labelleft=False)

        note = "grey: solid" + ("; dashed lines: z = %g and %g mm" % tuple(a.mark) if a.mark and axis != 2 else "")
        pale = lambda n: n in ("sel_CD", "wall_rate")

        # all fields at one time
        for j, i in enumerate(it):
            ncol, nrow, size = panel_size(len(fields))
            fig, axs = plt.subplots(nrow, ncol, figsize=size, squeeze=False, constrained_layout=True)
            axs = axs.ravel()
            for ax, (n, lab) in zip(axs, fields):
                norm, cmap = norms[n]
                im = draw(ax, pl, pl.image(snaps.values(i, n, pl.cells)), norm, cmap, a.mark, pale(n))
                ax.set_title(NICE.get(lab, lab), loc="left")
                colourbar(fig, im, ax, UNIT.get(n, "mass fraction"))
            label_axes(axs, len(fields), ncol)
            fig.suptitle("%s, t = %.2f s    (%s)" % (pl.title(), snaps.times[i], note), fontsize=9.5, x=0.01, ha="left")
            name = os.path.join(out, "slice_%s_t%07.3f.png" % (pl.tag(), snaps.times[i]))
            fig.savefig(name); plt.close(fig)
            frames.append(name)
            print("  written %s" % name, flush=True)

        if a.gif and len(frames) > 1:
            try:
                from PIL import Image
                pics = [Image.open(f).convert("RGB") for f in frames]
                w, h = pics[0].size; sc = min(1.0, 1400.0 / max(w, h))
                pics = [q.resize((int(w * sc), int(h * sc))) for q in pics]
                name = os.path.join(out, "slice_%s.gif" % pl.tag())
                pics[0].save(name, save_all=True, append_images=pics[1:], duration=200, loop=0)
                print("  written %s (%d frames)" % (name, len(pics)), flush=True)
            except Exception as e:
                print("  no animation: %s" % e, flush=True)

        # one field at the montage times
        if len(it_m) > 1:
            it_keep, it = it, it_m
            for n, lab in fields:
                ncol, nrow, size = panel_size(len(it))
                fig, axs = plt.subplots(nrow, ncol, figsize=size, squeeze=False, constrained_layout=True)
                axs = axs.ravel()
                norm, cmap = norms[n]
                for ax, j in zip(axs, range(len(it))):
                    im = draw(ax, pl, imgs[n][j], norm, cmap, a.mark, pale(n))
                    ax.set_title("t = %.2f s" % snaps.times[it[j]], loc="left", fontweight="normal")
                label_axes(axs, len(it), ncol)
                colourbar(fig, im, list(axs[:len(it)]), UNIT.get(n, "mass fraction"))
                fig.suptitle("%s, %s    (%s)" % (NICE.get(lab, lab), pl.title(), note), fontsize=9.5, x=0.01, ha="left")
                name = os.path.join(out, "%s_%s_times.png" % (re.sub(r"[^A-Za-z0-9_]+", "_", lab), pl.tag()))
                fig.savefig(name); plt.close(fig)
                print("  written %s" % name, flush=True)
            it = it_keep


if __name__ == "__main__":
    main()

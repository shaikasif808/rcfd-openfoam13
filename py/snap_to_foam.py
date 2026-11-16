#!/usr/bin/env python3
"""
snap_to_foam.py - rCFD run snapshots -> OpenFOAM time directories (for ParaView's
                  OpenFOAM reader or OpenFOAM post-processing).

    python3 snap_to_foam.py --run CASE --case CASE --out CASE/rCFD_results [--copy-mesh] [--every K] [--binary] [--skip vof]

    --run    directory with out/snap_%06d.bin and out/fields.txt (written by rcfd_run)
    --case   the OpenFOAM case the mesh comes from (constant/polyMesh, ascii or binary)
    --out    results case that is created:
                 <out>/constant/polyMesh     link to CASE/constant/polyMesh (copy with --copy-mesh)
                 <out>/system/controlDict    minimal
                 <out>/<time>/rCFD_<name>    volScalarField per snapshot field: rCFD_vof.<phase>,
                                             rCFD_c<i>.<phase> (data i of the phase) and the
                                             arrays made with rcfd_user_field() in the user header
                 <out>/<name of out>.foam    empty file: open it in ParaView

boundaryField: every patch gets `type zeroGradient;`, constraint patches (empty, wedge,
symmetry, cyclic, processor ...) get their own type, as OpenFOAM requires.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from of_reader import find_file, read_raw, parse_foam_header, read_label_list, _read_boundary, time_name  # noqa: E402
from rcfd_snap import list_snapshots, read_snapshot, field_names  # noqa: E402

BANNER = (
    "/*--------------------------------*- C++ -*----------------------------------*\\\n"
    "  =========                 |\n"
    "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"
    "   \\\\    /   O peration     |\n"
    "    \\\\  /    A nd           | rCFD results written by snap_to_foam.py\n"
    "     \\\\/     M anipulation  |\n"
    "\\*---------------------------------------------------------------------------*/\n")
SEP = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
END = "\n// ************************************************************************* //\n"

CONSTRAINT_TYPES = {"empty", "wedge", "symmetry", "symmetryPlane", "cyclic", "cyclicAMI", "cyclicACMI",
                    "cyclicSlip", "processor", "processorCyclic", "nonuniformTransformCyclic"}


def header(cls, obj, location):
    return (BANNER + "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       %s;\n"
            "    location    \"%s\";\n    object      %s;\n}\n" % (cls, location, obj) + SEP + "\n")


def boundary_field_text(patches):
    s = "boundaryField\n{\n"
    for p in patches:
        t = p["type"] if p["type"] in CONSTRAINT_TYPES else "zeroGradient"
        s += "    %s\n    {\n        type            %s;\n    }\n" % (p["name"], t)
    return s + "}\n"


def write_scalar_field(path, name, tname, values, bf_text, binary=False):
    if binary:
        # OpenFOAM binary list: the values as 8-byte floats between the brackets
        head = header("volScalarField", name, tname).replace(
            "format      ascii;", "format      binary;\n    arch        \"LSB;label=32;scalar=64\";")
        with open(str(path), "wb") as f:
            f.write((head + "dimensions      [0 0 0 0 0 0 0];\n\n"
                     "internalField   nonuniform List<scalar> \n%d\n(" % len(values)).encode())
            f.write(np.asarray(values, dtype="<f8").tobytes())
            f.write((")\n;\n\n" + bf_text + END).encode())
        return
    with open(str(path), "w") as f:
        f.write(header("volScalarField", name, tname))
        f.write("dimensions      [0 0 0 0 0 0 0];\n\n")
        f.write("internalField   nonuniform List<scalar> \n%d\n(\n" % len(values))
        np.savetxt(f, values, fmt="%.8g")
        f.write(")\n;\n\n" + bf_text + END)


def read_boundary_any_format(path):
    """polyMesh/boundary -> list of {name, type}. The patch list is text also in a mesh
    written with writeFormat binary, so the format entry of the header is not checked."""
    import re
    raw = read_raw(path)
    hdr = parse_foam_header(raw, path, require_ascii=False)
    body = re.sub(rb"//[^\n]*", b"", raw[hdr["_end"]:])
    body = re.sub(rb"/\*.*?\*/", b"", body, flags=re.S)
    m = re.search(rb"(\d+)\s*\(", body)
    if m is None:
        sys.exit("%s: cannot find the patch list" % path)
    n, pos, patches = int(m.group(1)), m.end(), []
    while len(patches) < n:
        mm = re.compile(rb"\s*([^\s{}();]+)\s*\{").match(body, pos)
        if mm is None:
            sys.exit("%s: cannot parse patch %d" % (path, len(patches)))
        depth, i = 1, mm.end()
        while depth > 0:
            j_open, j_close = body.find(b"{", i), body.find(b"}", i)
            if j_close < 0:
                sys.exit("%s: unbalanced braces" % path)
            if 0 <= j_open < j_close:
                depth, i = depth + 1, j_open + 1
            else:
                depth, i = depth - 1, j_close + 1
        k = re.search(rb"\btype\s+([^;\s]+)\s*;", body[mm.end():i - 1])
        if k is None:
            sys.exit("%s: patch %s has no type" % (path, mm.group(1).decode()))
        patches.append({"name": mm.group(1).decode(), "type": k.group(1).decode()})
        pos = i
    return patches


def unique_time_names(times):
    for digits in range(6, 16):
        names = [time_name(t, digits) for t in times]
        if len(set(names)) == len(names):
            return names
    raise RuntimeError("snapshot times are not distinct")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="run directory (snapshots in RUN/out)")
    ap.add_argument("--case", required=True, help="OpenFOAM case with constant/polyMesh")
    ap.add_argument("--out", required=True, help="results case to create")
    ap.add_argument("--copy-mesh", action="store_true", help="copy constant/polyMesh instead of linking it")
    ap.add_argument("--every", type=int, default=1, help="only every K-th snapshot")
    ap.add_argument("--binary", action="store_true", help="write the fields in binary format (large meshes: smaller and faster)")
    ap.add_argument("--skip", default="", help="comma-separated beginnings of field names that are not written, e.g. vof")
    args = ap.parse_args()
    skip = tuple("rCFD_" + x.strip() for x in args.skip.split(",") if x.strip())

    src_pm = Path(args.case).resolve() / "constant" / "polyMesh"
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        if find_file(src_pm / name) is None:
            sys.exit("%s: no '%s'" % (src_pm, name))
    patches = read_boundary_any_format(src_pm / "boundary")
    n_cells = None          # taken from the first snapshot; all snapshots must agree
    bf_text = boundary_field_text(patches)

    snaps = list_snapshots(args.run)[::max(1, args.every)]
    if not snaps:
        sys.exit("no snap_*.bin in %s/out" % args.run)

    out = Path(args.out)
    (out / "constant").mkdir(parents=True, exist_ok=True)
    (out / "system").mkdir(parents=True, exist_ok=True)
    dst_pm = out / "constant" / "polyMesh"
    if dst_pm.is_symlink() or dst_pm.is_file():
        dst_pm.unlink()
    elif dst_pm.is_dir():
        shutil.rmtree(str(dst_pm))
    if args.copy_mesh:
        shutil.copytree(str(src_pm), str(dst_pm))
    else:
        os.symlink(str(src_pm), str(dst_pm), target_is_directory=True)

    data = []
    for idx, path in snaps:
        t, fields = read_snapshot(path)
        if n_cells is None:
            n_cells = fields.shape[1]
        if fields.shape[1] != n_cells:
            sys.exit("%s has %d cells, the first snapshot %d" % (path, fields.shape[1], n_cells))
        data.append((t, path))
    names = unique_time_names([t for t, _ in data])

    with open(str(out / "system" / "controlDict"), "w") as f:
        f.write(header("dictionary", "controlDict", "system"))
        f.write("application     rCFD;\n\nstartFrom       startTime;\n\nstartTime       %s;\n\nstopAt          endTime;\n\n"
                "endTime         %s;\n\ndeltaT          1;\n\nwriteControl    adjustableRunTime;\n\nwriteInterval   1;\n\n"
                "adjustTimeStep  yes;\n\nwriteFormat     ascii;\n\nwritePrecision  8;\n\nwriteCompression off;\n\n"
                "timeFormat      general;\n\ntimePrecision   10;\n" % (names[0], names[-1]))
        f.write(END)
    foam_file = out / (out.resolve().name + ".foam")
    foam_file.touch()

    for (t, path), tname in zip(data, names):
        _, fields = read_snapshot(path)
        d = out / tname
        d.mkdir(exist_ok=True)
        fnames = field_names(len(fields), path.parent)
        written = []
        for name, values in zip(fnames, fields):
            if skip and name.startswith(skip):
                continue
            write_scalar_field(d / name, name, tname, values, bf_text, args.binary)
            written.append(name)
        print("%s -> %s/{%s}   t = %.10g" % (path.name, d, ",".join(written), t), flush=True)
    print("results case %s: %d time directories (%s .. %s), %d cells, %d patches (%s); mesh %s; open %s in ParaView"
          % (out, len(names), names[0], names[-1], n_cells, len(patches),
             ", ".join("%s: %s" % (p["name"], p["type"] if p["type"] in CONSTRAINT_TYPES else "zeroGradient") for p in patches),
             "copied" if args.copy_mesh else "linked", foam_file))


if __name__ == "__main__":
    main()

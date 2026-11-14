#!/usr/bin/env python3
"""
rcfd_snap.py - read / write rCFD run snapshots CASE/out/snap_%06d.bin written by rcfd_run:

    int64 n_cells, int32 n_fields, float64 time,
    n_fields x float32[n_cells]          little endian, OpenFOAM cell order
    names of the fields: CASE/out/fields.txt

The header is read as 20 packed bytes; a header padded to 24 bytes (a C struct written
in one piece) is recognised from the file size and read as well.
"""
import os
import re
from pathlib import Path

import numpy as np

_RE_SNAP = re.compile(r"^snap_(\d+)\.bin$")


def list_snapshots(run_or_dir):
    """Sorted list of (index, path) of the snapshots in RUN/out (or in the directory itself)."""
    d = Path(run_or_dir)
    if (d / "out").is_dir() and not any(_RE_SNAP.match(f) for f in os.listdir(str(d))):
        d = d / "out"
    hits = []
    for name in os.listdir(str(d)):
        m = _RE_SNAP.match(name)
        if m:
            hits.append((int(m.group(1)), d / name))
    hits.sort()
    return hits


def read_snapshot(path):
    """-> (time, fields[n_fields][n_cells] float32)."""
    size = os.path.getsize(str(path))
    with open(str(path), "rb") as f:
        head = f.read(24)
    if len(head) < 20:
        raise ValueError("%s: too short for a snapshot" % path)
    n = int(np.frombuffer(head, dtype="<i8", count=1)[0])
    nf = int(np.frombuffer(head, dtype="<i4", count=1, offset=8)[0])
    if n <= 0 or nf <= 0 or nf > 1000:
        raise ValueError("%s: implausible header (n_cells %d, n_fields %d)" % (path, n, nf))
    body = 4 * n * nf
    if size == 20 + body:
        t_off, d_off = 12, 20
    elif size == 24 + body:
        t_off, d_off = 16, 24
    else:
        raise ValueError("%s: file size %d does not fit n_cells %d x n_fields %d (expected %d)" % (path, size, n, nf, 20 + body))
    t = float(np.frombuffer(head, dtype="<f8", count=1, offset=t_off)[0])
    with open(str(path), "rb") as f:
        f.seek(d_off)
        data = np.fromfile(f, dtype="<f4", count=n * nf).reshape(nf, n)
    return t, data


def write_snapshot(path, t, fields):
    """fields: sequence of arrays (n,) -> packed snapshot file."""
    fields = [np.ascontiguousarray(a, dtype="<f4") for a in fields]
    n = len(fields[0])
    with open(str(path), "wb") as f:
        np.array([n], dtype="<i8").tofile(f)
        np.array([len(fields)], dtype="<i4").tofile(f)
        np.array([t], dtype="<f8").tofile(f)
        for a in fields:
            a.tofile(f)


def field_names(n_fields, snap_dir=None):
    """Names of the fields of a snapshot: out/fields.txt written by rcfd_run (one name per
    line, e.g. vof.particles, c0.particles, vf_diam_0), prefixed with rCFD_."""
    names = []
    if snap_dir is not None:
        f = Path(snap_dir) / "fields.txt"
        if f.is_file():
            names = ["rCFD_" + line.strip() for line in open(str(f)) if line.strip()]
    names = names[:n_fields]
    return names + ["rCFD_field%d" % i for i in range(len(names), n_fields)]

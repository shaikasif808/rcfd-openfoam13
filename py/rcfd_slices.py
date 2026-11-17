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

Curves (no snapshots needed, from the tables the run writes: rCFD_*_monitor.out and
rCFD_*_profile_history.out), also into ./slices/:
    plot_outlet_time.png                    species at the outlet and their mass in the domain against time
    plot_stations_time.png                  species against time at positions along the flow (--stations, mm)
    plot_conversion_selectivity_time.png    converted share of the feed and selectivity at these positions
    plot_profiles.png                       species, conversion and selectivity along the flow at six times
    plot_maps_z_time.png                    the same as maps over position and time
    curves_stations.csv                     the numbers of the curves at the positions
--plots-only writes only these (seconds; works in a folder rcfd_run_<date>_<time>/ as well),
--no-plots leaves them out. Conversion = share of the feed species that has reacted; selectivity =
desired / (desired + waste) products: --feed A --desired C D --waste E are the defaults.

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


# colours of the curves: species in a fixed order, positions and times as steps of one blue
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
INK, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#efeee9"


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


# ----------------------------------------------------------------------------- curves
def table_names(path, n):
    with open(path) as f:
        m = re.search(r"mean (.+?) in the", f.readline())
    names = m.group(1).split() if m else []
    return names if len(names) == n else [chr(ord("A") + i) for i in range(n)]


def read_monitor(path):
    """time, mean of every species at the outlet, mass of every species in the domain."""
    d = np.loadtxt(path, comments="#", ndmin=2)
    n = (d.shape[1] - 1) // 2
    d = d[np.isfinite(d[:, 0]) & (d[:, 0] >= 0)][::-1]
    d = d[np.unique(d[:, 0], return_index=True)[1]]            # a time written twice: the later line
    return table_names(path, n), d[:, 0], d[:, 1:1 + n], d[:, 1 + n:1 + 2 * n]


def read_profiles(path):
    """states of the layer averages: times, positions z [m], values[time, layer, species]."""
    times, blocks = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"time\s+([-+.\deE]+)", line)
                times.append(float(m.group(1)) if m else np.nan); blocks.append([])
            elif line.strip() and blocks:
                try:
                    blocks[-1].append([float(v) for v in line.split()])
                except ValueError:
                    pass
    width = max((len(r) for b in blocks for r in b), default=0)
    state = {}
    for t, b in zip(times, blocks):
        b = np.array([r for r in b if len(r) == width])        # a line that is being written is left out
        if b.size and np.isfinite(t) and t >= 0:
            state[t] = b
    if not state:
        return [], np.zeros(0), np.zeros(0), np.zeros((0, 0, 0))
    n = width - 2
    z = np.unique(np.concatenate([b[:, 0] for b in state.values()]))
    t = np.array(sorted(state))
    val = np.full((len(t), len(z), n), np.nan)
    for k, tk in enumerate(t):
        val[k, np.searchsorted(z, state[tk][:, 0])] = state[tk][:, 1:1 + n]
    return table_names(path, n), t, z, val


def shares(val, names, feed, desired, waste):
    """converted share of the feed and selectivity; nan where there is nothing to divide."""
    ix = {n: i for i, n in enumerate(names)}
    conv = sel = None
    with np.errstate(invalid="ignore", divide="ignore"):
        if feed in ix:
            total = val.sum(-1)
            conv = np.where(total > 1e-3, (total - val[..., ix[feed]]) / total, np.nan)
        if desired and waste and all(n in ix for n in desired + waste):
            good = sum(val[..., ix[n]] for n in desired); bad = sum(val[..., ix[n]] for n in waste)
            sel = np.where(good + bad > 1e-4, good / (good + bad), np.nan)
    return conv, sel


def style(ax, xlabel=None, ylabel=None, title=None):
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b5b4ae"); ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=MUTED, length=2.5, width=0.5, labelsize=7.5)
    if xlabel: ax.set_xlabel(xlabel, color=MUTED)
    if ylabel: ax.set_ylabel(ylabel, color=MUTED)
    if title: ax.set_title(title, loc="left", color=INK)


def end_labels(ax, x, ys, labels, colours):
    """Names at the right ends of the curves, moved apart where they would overlap."""
    lo, hi = ax.get_ylim()
    gap = 0.058 * (hi - lo)
    item = sorted((y, l, c) for y, l, c in zip(ys, labels, colours) if np.isfinite(y))
    pos = []
    for y, _, _ in item:
        pos.append(max(y, pos[-1] + gap) if pos else max(y, lo + 0.3 * gap))
    over = pos[-1] - (hi - 0.3 * gap) if pos else 0.0
    if over > 0:                                           # pushed out at the top: move the stack down
        pos = [q - over for q in pos]
        for k in range(len(pos) - 2, -1, -1):
            pos[k] = min(pos[k], pos[k + 1] - gap)
    x0, x1 = ax.get_xlim()
    for (y, l, c), q in zip(item, pos):
        ax.annotate(l, xy=(x, y), xytext=(x + 0.035 * (x1 - x0), q), color=INK, fontsize=7.5, va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.7, shrinkA=0, shrinkB=1.5), annotation_clip=False)


def ramp(n):
    return [RAMP[int(round(k))] for k in np.linspace(0, len(RAMP) - 1, max(n, 2))][:n] if n > 1 else [RAMP[4]]


def curves(case, out, a):
    """Pictures of the tables the run writes. Returns the number of pictures."""
    mon = sorted(glob.glob(os.path.join(case, "rCFD_*_monitor.out")))
    mon = [f for f in mon if "mean_height" not in f]
    his = sorted(glob.glob(os.path.join(case, "rCFD_*_profile_history.out")))
    if not mon and not his:
        print("no rCFD_*_monitor.out or rCFD_*_profile_history.out in %s: no curves" % os.path.abspath(case), flush=True)
        return 0
    done = []

    def save(fig, name):
        path = os.path.join(out, name)
        fig.savefig(path); plt.close(fig); done.append(path)
        print("  written %s" % path, flush=True)

    def legend(fig, handles, labels, title=None):
        fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), frameon=False, fontsize=7.5,
                   handlelength=1.6, columnspacing=1.2, title=title, title_fontsize=7.5)

    # --- outlet and content of the domain against time
    if mon:
        names, t, outlet, mass = read_monitor(mon[0])
        if t.size > 1:
            col = [SERIES[i % len(SERIES)] for i in range(len(names))]
            fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.5), constrained_layout=True)
            for ax, v, unit, title in ((axs[0], outlet, "mass fraction", "At the outlet"),
                                       (axs[1], 1e6 * mass, "mass in mg", "Content of the domain")):
                hs = [ax.plot(t, v[:, i], color=col[i], lw=1.7)[0] for i in range(len(names))]
                ax.plot(t, v.sum(1), color=MUTED, lw=0.9, ls=(0, (4, 2)))
                style(ax, "time in s", unit, title)
                ax.set_xlim(t[0], t[-1] + 0.09 * (t[-1] - t[0])); ax.set_ylim(bottom=0)
                end_labels(ax, t[-1], list(v[-1]) + [v[-1].sum()], names + ["sum"], col + [MUTED])
            legend(fig, hs, names)
            fig.suptitle("Species against time    (dashed: sum of all species)", fontsize=9.5, x=0.01, ha="left")
            save(fig, "plot_outlet_time.png")

    if not his:
        return len(done)
    names, t, z, val = read_profiles(his[0])
    if t.size < 2:
        print("  %s holds fewer than two states: no curves along the flow" % his[0], flush=True)
        return len(done)
    zmm = 1e3 * z
    col = [SERIES[i % len(SERIES)] for i in range(len(names))]
    conv, sel = shares(val, names, a.feed or names[0], a.desired, a.waste)
    extra = [(v, lab) for v, lab in ((conv, "Converted share of %s" % (a.feed or names[0])),
                                     (sel, "Selectivity (%s) / (%s)" % (" + ".join(a.desired), " + ".join(a.desired + a.waste))))
             if v is not None]

    if a.stations:
        want = a.stations
    elif a.mark:
        want = [0.5 * (a.mark[0] + a.mark[1]), a.mark[1], 0.5 * (a.mark[1] + zmm[-1]), zmm[-1]]
    else:
        want = list(np.linspace(zmm[0], zmm[-1], 5)[1:])
    ks = sorted(set(int(np.argmin(np.abs(zmm - w))) for w in want))
    st_lab = ["z = %.4g mm" % zmm[k] for k in ks]
    span = (lambda ax: ax.axvspan(a.mark[0], a.mark[1], color=BAND, lw=0, zorder=0)) if a.mark else (lambda ax: None)
    note = "    (shaded: z = %g to %g mm)" % tuple(a.mark) if a.mark else ""

    # --- species against time at the positions
    ncol = min(len(ks), 4); nrow = -(-len(ks) // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.75 * ncol + 0.4, 2.55 * nrow + 0.55), squeeze=False, sharey=True,
                            constrained_layout=True)
    top = max(1e-12, np.nanmax(val[:, ks])) * 1.05
    for j, ax in enumerate(axs.ravel()):
        if j >= len(ks):
            ax.set_axis_off(); continue
        hs = [ax.plot(t, val[:, ks[j], i], color=col[i], lw=1.7)[0] for i in range(len(names))]
        style(ax, "time in s" if j + ncol >= len(ks) else None, "mass fraction" if j % ncol == 0 else None, st_lab[j])
        ax.set_xlim(t[0], t[-1] + 0.16 * (t[-1] - t[0])); ax.set_ylim(0, top)
        end_labels(ax, t[-1], val[-1, ks[j]], names, col)
    legend(fig, hs, names)
    fig.suptitle("Species against time at positions along the flow    (averages over the cross-section)", fontsize=9.5,
                 x=0.01, ha="left")
    save(fig, "plot_stations_time.png")

    # --- conversion and selectivity against time at the positions
    if extra:
        cs = ramp(len(ks))
        fig, axs = plt.subplots(1, len(extra), figsize=(4.8 * len(extra), 3.5), squeeze=False, constrained_layout=True)
        for ax, (v, lab) in zip(axs.ravel(), extra):
            hs = [ax.plot(t, v[:, k], color=c, lw=1.7)[0] for k, c in zip(ks, cs)]
            style(ax, "time in s", "share", lab)
            ax.set_xlim(t[0], t[-1] + 0.2 * (t[-1] - t[0])); ax.set_ylim(0, 1)
            end_labels(ax, t[-1], [v[np.isfinite(v[:, k]), k][-1] if np.isfinite(v[:, k]).any() else np.nan for k in ks], st_lab, cs)
        legend(fig, hs, st_lab)
        fig.suptitle("Conversion and selectivity against time    (a curve starts when the species arrive)", fontsize=9.5,
                     x=0.01, ha="left")
        save(fig, "plot_conversion_selectivity_time.png")

    # --- along the flow at a few times
    its = sorted(set(int(round(k)) for k in np.linspace(0, len(t) - 1, min(6, len(t)))))
    cs = ramp(len(its))
    panels = [(val[..., i], names[i], "mass fraction") for i in range(len(names))] + [(v, lab, "share") for v, lab in extra]
    ncol = min(len(panels), 4); nrow = -(-len(panels) // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.75 * ncol + 0.4, 2.3 * nrow + 0.6), squeeze=False, sharex=True,
                            constrained_layout=True)
    for j, ax in enumerate(axs.ravel()):
        if j >= len(panels):
            ax.set_axis_off(); continue
        v, lab, unit = panels[j]
        span(ax)
        hs = [ax.plot(zmm, v[i], color=c, lw=1.5)[0] for i, c in zip(its, cs)]
        style(ax, "z in mm" if j + ncol >= len(panels) else None, unit, lab)
        ax.tick_params(labelbottom=True)
        ax.set_xlim(zmm[0], zmm[-1]); ax.set_ylim(0, 1 if unit == "share" else None)
    legend(fig, hs, ["%.2f s" % t[i] for i in its], title="time (light: early, dark: late)")
    fig.suptitle("Along the flow: averages over the cross-section%s" % note, fontsize=9.5, x=0.01, ha="left")
    save(fig, "plot_profiles.png")

    # --- maps over position and time
    edges = lambda c: np.concatenate([[c[0] - 0.5 * (c[1] - c[0])], 0.5 * (c[1:] + c[:-1]), [c[-1] + 0.5 * (c[-1] - c[-2])]])
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.9 * ncol + 0.4, 2.3 * nrow + 0.6), squeeze=False, sharex=True, sharey=True,
                            constrained_layout=True)
    for j, ax in enumerate(axs.ravel()):
        if j >= len(panels):
            ax.set_axis_off(); continue
        v, lab, unit = panels[j]
        cm = matplotlib.colormaps["viridis"].copy() if hasattr(matplotlib, "colormaps") else plt.get_cmap("viridis").copy()
        cm.set_bad(PALE)
        hi = 1.0 if unit == "share" else nice_max(np.nanpercentile(v, 99.9)) if np.isfinite(v).any() else 1.0
        im = ax.pcolormesh(edges(t), edges(zmm), np.ma.masked_invalid(v).T, cmap=cm, norm=Normalize(0.0, hi), rasterized=True)
        if a.mark:
            for zz in a.mark:
                ax.axhline(zz, color="white", lw=0.9, ls=(0, (4, 3))); ax.axhline(zz, color="black", lw=0.35, ls=(0, (4, 3)))
        ax.set_title(lab, loc="left", color=INK)
        ax.tick_params(colors=MUTED, length=2.5, width=0.5, labelsize=7.5, labelbottom=True)
        for s in ax.spines.values():
            s.set_linewidth(0.5)
        if j + ncol >= len(panels): ax.set_xlabel("time in s", color=MUTED)
        if j % ncol == 0: ax.set_ylabel("z in mm (flow upwards)", color=MUTED)
        colourbar(fig, im, ax, unit)
    fig.suptitle("Position against time: averages over the cross-section%s" % (
        "    (dashed lines: z = %g and %g mm; pale: nothing there yet)" % tuple(a.mark) if a.mark else ""),
        fontsize=9.5, x=0.01, ha="left")
    save(fig, "plot_maps_z_time.png")

    # --- the numbers of the curves at the positions
    path = os.path.join(out, "curves_stations.csv")
    with open(path, "w") as f:
        f.write("time_s,z_mm," + ",".join(names) + "".join(",%s" % n for n in ("conversion", "selectivity")[:len(extra)]) + "\n")
        for k in ks:
            for i in range(len(t)):
                row = ["%.6g" % t[i], "%.6g" % zmm[k]] + ["%.6g" % x for x in val[i, k]]
                row += ["" if not np.isfinite(v[i, k]) else "%.6g" % v[i, k] for v, _ in extra]
                f.write(",".join(row) + "\n")
    print("  written %s" % path, flush=True)
    return len(done)


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
    p.add_argument("--plots-only", action="store_true", help="only the curves from the tables of the run, no slices")
    p.add_argument("--no-plots", action="store_true", help="no curves")
    p.add_argument("--stations", type=float, nargs="*", default=None, help="z positions [mm] of the curves against time")
    p.add_argument("--feed", default=None, help="feed species of the conversion (default: the first)")
    p.add_argument("--desired", nargs="*", default=["C", "D"], help="desired products of the selectivity")
    p.add_argument("--waste", nargs="*", default=["E"], help="unwanted products of the selectivity")
    a = p.parse_args()

    out = a.out or os.path.join(a.case, "slices")
    os.makedirs(out, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                         "axes.labelsize": 8, "figure.dpi": 100, "savefig.dpi": a.dpi, "savefig.bbox": "tight",
                         "savefig.pad_inches": 0.06})

    if not a.no_plots:
        curves(a.case, out, a)
    if a.plots_only:
        return

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

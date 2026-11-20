#!/usr/bin/env python3
"""Recurrence statistics of an rCFD database: is it recurrent, long enough, and how do its shifts mix?

    python3 rcfd_recurrence.py                                   # in the case folder (the one holding data/)
    python3 rcfd_recurrence.py --plane x=0 z=0.375               # planes of your choice
    python3 rcfd_recurrence.py --case A B C                      # several databases, table at the end
    python3 rcfd_recurrence.py --every 4 --no-shifts             # quick look

Reads the database only (data/rec/matrix_0_0_0, data/vof, data/c2c, data/topo, data/rcfd_foam_prep.txt)
and writes the folder recurrence/ into the case folder (a few MB). Nothing in data/ is changed.

    summary.txt             all numbers, with one line of explanation each
    recurrence_matrix.png   the library's recurrence matrix, distance against time lag
    recurrence_lag.csv      distance against time lag: library matrix, sampled cells, planes
    cells.png, cells.csv    volume fraction at random bed cells: series, autocorrelation, spectrum
    bed.png, bed.csv        bed height and phase volume per chamber over the database (drift)
    planes.png              time mean and fluctuation of the volume fraction on the planes
    shifts.png, shifts.csv  displacement per frame and numerical diffusivity of the shifts, per phase

What the numbers mean
    distance          sum over cells of |difference of the volume fraction| between two frames (the
                      library's norm); 'plateau' is its value for frames far apart in time
    decorrelation     time lag at which the distance reaches 90 % of the plateau
    return            deepest minimum of the distance after that: the dominant period, if there is one
    partner           for every frame the most similar frame outside its own neighbourhood; given
                      relative to the distance between two successive frames. Near 1: the database
                      holds states as close as successive frames. Large: no real recurrence.
    seam              distance last frame -> first frame (replay jumps there), same scale
    halves            time mean of the first half against the second half, rms over the bed cells,
                      relative to the rms fluctuation: below about 0.1 the mean field is converged
    D_C2C             numerical diffusivity of the shift, variance of the source positions of a cell
                      / (2 * frame interval), as in the library (prep_set_artificial_diffusivity) and
                      in Lumetzberger, Pirker, Lichtenegger, Chem. Eng. Sci. 311 (2025) 121624, eq. 39

Needs numpy and matplotlib.
"""
import argparse, os, sys, time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PolyCollection

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f8fe", "#86b6ef", "#2a78d6", "#0d366b"])
AXES = {"x": 0, "y": 1, "z": 2}


def style(ax, xlabel, ylabel, title):
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b5b4ae"); ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=MUTED, length=2.5, width=0.5, labelsize=7.5)
    ax.set_xlabel(xlabel, color=MUTED); ax.set_ylabel(ylabel, color=MUTED); ax.set_title(title, loc="left", color=INK)


# ----------------------------------------------------------------------------- reading the database

def read_prep(case):
    d = {}
    with open(os.path.join(case, "data", "rcfd_foam_prep.txt")) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1); d[k.strip()] = v.strip()
    return dict(dt=float(d["frame_interval"]), n_frames=int(d["n_frames"]), n_phases=int(d["n_phases"]),
                n_ranks=int(d["n_ranks"]), phases=d.get("phases", "").split())


def read_topo(case, rank):
    """centres [n, 3], volumes, host rank of every cell of the rank (internal cells first), global numbers"""
    with open(os.path.join(case, "data", "topo", "topo_%d" % rank)) as f:
        n_cells = int(f.readline().split()[0])
        a = np.fromfile(f, sep=" ", count=5 * n_cells).reshape(n_cells, 5)
    gid = np.fromfile(os.path.join(case, "data", "topo", "gid_%d.bin" % rank), dtype="<i4")
    return a[:, :3].copy(), a[:, 3].copy(), a[:, 4].astype(int), gid


def read_vof(case, phase, rank, n_int, every=1):
    """volume fraction of one phase on one rank: [frames, internal cells], float32.
    Per frame the file holds one line with the cell count of the rank and one value per internal cell."""
    path = os.path.join(case, "data", "vof", "vof_0_%d_%d" % (phase, rank))
    a = np.fromfile(path, sep=" ", dtype=np.float32)
    n = a.size // (n_int + 1)
    if n * (n_int + 1) != a.size:
        print("  note: %s holds %d values, not a multiple of %d: the last frame is incomplete and left out" % (path, a.size, n_int + 1))
    return a[:n * (n_int + 1)].reshape(n, n_int + 1)[::every, 1:]


def shift_frames(case, phase, rank, frames, n_ranks):
    """yields (frame, c0, c1, w, source rank) for the wanted frames; walks the binary file without reading the rest"""
    path = os.path.join(case, "data", "c2c", "c2c_0_%d_%d" % (phase, rank))
    ext = np.dtype([("c0", "<i4"), ("r", "u1"), ("c1", "<i4"), ("w", "u1")]) if n_ranks <= 256 \
        else np.dtype([("c0", "<i4"), ("r", "<i4"), ("c1", "<i4"), ("w", "u1")])
    inn = np.dtype([("c0", "<i4"), ("c1", "<i4"), ("w", "u1")])
    wanted = set(int(k) for k in frames); last = max(wanted); size = os.path.getsize(path)
    with open(path, "rb") as f:
        k = 0
        while k <= last and f.tell() + 8 <= size:
            n_int, n_ext = np.fromfile(f, dtype="<i4", count=2)
            if k in wanted:
                a = np.fromfile(f, dtype=inn, count=n_int); b = np.fromfile(f, dtype=ext, count=n_ext)
                yield (k, np.concatenate([a["c0"], b["c0"]]), np.concatenate([a["c1"], b["c1"]]),
                       np.concatenate([a["w"], b["w"]]).astype(float) / 255.0,
                       np.concatenate([np.full(n_int, rank), b["r"].astype(int)]))
            else:
                f.seek(int(n_int) * inn.itemsize + int(n_ext) * ext.itemsize, 1)
            k += 1


# ----------------------------------------------------------------------------- analysis

def l1_matrix(a):
    """distance matrix as the library builds it: sum over the samples of |a_i - a_j|"""
    n = a.shape[0]; m = np.zeros((n, n))
    for i in range(n):
        m[i] = np.abs(a - a[i]).sum(axis=1)
    return m


def lag_curve(m):
    n = m.shape[0]
    return np.array([np.diagonal(m, k).mean() for k in range(n)])


def matrix_numbers(m, dt):
    """key numbers of a distance matrix"""
    n = m.shape[0]; lag = lag_curve(m)
    plateau = lag[n // 4: n // 2 + 1].mean() if n >= 8 else lag[1:].mean()
    step = lag[1]
    above = np.nonzero(lag >= 0.9 * plateau)[0]
    k_dec = int(above[0]) if above.size else n - 1
    k_dec = max(k_dec, 1)
    hi = max(n // 2, k_dec + 2)
    k_ret = k_dec + int(np.argmin(lag[k_dec: hi])) if hi > k_dec else k_dec
    depth = lag[k_ret] / plateau
    # most similar frame outside the own neighbourhood
    idx = np.arange(n)
    far = np.abs(idx[:, None] - idx[None, :]) > k_dec
    best = np.where(far, m, np.inf).min(axis=1)
    best = best[np.isfinite(best)]
    partner = np.median(best) / step if best.size and step > 0 else float("nan")
    # best loop: frame j follows frame i smoothly when frame j - 1 is close to frame i
    a0, a1 = max(1, n // 10), max(n - n // 10, n // 2)
    loop = None
    if a1 < n and a0 >= 1 and n >= 20:
        sub = m[a1:n, 0:a0]                      # [i - a1, j - 1]
        p = np.unravel_index(np.argmin(sub), sub.shape)
        loop = (int(p[1]) + 1, int(p[0]) + a1, sub[p] / step)
    return dict(lag=lag, plateau=plateau, step=step, t_dec=k_dec * dt, t_ret=k_ret * dt, depth=depth,
                partner=partner, seam=m[n - 1, 0] / step, seam_plateau=m[n - 1, 0] / plateau, loop=loop)


def autocorrelation(s):
    """mean normalised autocorrelation of the columns of s"""
    f = s - s.mean(axis=0); n = f.shape[0]
    keep = f.std(axis=0) > 1e-8
    if not keep.any():
        return np.zeros(n), np.zeros(n // 2 + 1)
    f = f[:, keep]
    spec = np.abs(np.fft.rfft(f, n=2 * n, axis=0)) ** 2
    ac = np.fft.irfft(spec, axis=0)[:n] / (np.arange(n, 0, -1)[:, None])
    ac = ac / ac[0]
    power = (np.abs(np.fft.rfft(f * np.hanning(n)[:, None], axis=0)) ** 2).mean(axis=1)
    return ac.mean(axis=1), power


def analyse(case, a):
    t0 = time.time()
    info = read_prep(case)
    dt = info["dt"] * a.every; n_ranks = info["n_ranks"]
    phase = info["n_phases"] - 1
    if a.phase is not None:
        phase = info["phases"].index(a.phase) if a.phase in info["phases"] else int(a.phase)
    pname = info["phases"][phase] if phase < len(info["phases"]) else "phase %d" % phase
    out = os.path.join(case, "recurrence"); os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    print("%s: %d frames of %g s, %d ranks, phase '%s'" % (case, info["n_frames"], info["dt"], n_ranks, pname)); sys.stdout.flush()

    topo = [read_topo(case, r) for r in range(n_ranks)]
    lo = np.min([t[0].min(axis=0) for t in topo], axis=0); hi = np.max([t[0].max(axis=0) for t in topo], axis=0)
    up = AXES[a.up]
    horiz = [d for d in range(3) if d != up]
    along = AXES[a.along] if a.along else max(horiz, key=lambda d: hi[d] - lo[d])
    third = [d for d in horiz if d != along][0]
    planes = []
    for p in (a.plane or ["%s=%g" % ("xyz"[third], 0.5 * (lo[third] + hi[third])), "%s=%g" % ("xyz"[along], 0.5 * (lo[along] + hi[along]))]):
        ax_, v = p.split("="); planes.append((AXES[ax_], float(v), "%s=%.4g" % (ax_, float(v) if abs(float(v)) > 1e-12 else 0.0)))

    n_glob = sum(int((t[2] == r).sum()) for r, t in enumerate(topo))
    g_mean = np.zeros(n_glob); g_std = np.zeros(n_glob); g_h1 = np.zeros(n_glob); g_h2 = np.zeros(n_glob)
    g_xyz = np.zeros((n_glob, 3)); g_vol = np.zeros(n_glob)
    vol_t = None; mom_t = None; chamber_t = None
    cells, cell_xyz, plane_series, plane_ids = [], [], [[] for _ in planes], [[] for _ in planes]
    edges = np.linspace(lo[along], hi[along], a.chambers + 1)

    for r in range(n_ranks):
        xyz, vol, host, gid = topo[r]
        ni = int((host == r).sum())
        A = read_vof(case, phase, r, ni, a.every)
        F = A.shape[0]
        own = host[:ni] == r
        x, v, g = xyz[:ni][own], vol[:ni][own], gid[:ni][own]
        A = A[:, own]
        if vol_t is None:
            vol_t = np.zeros(F); mom_t = np.zeros(F); chamber_t = np.zeros((F, a.chambers))
        F = min(F, vol_t.size); A = A[:F]
        m = A.mean(axis=0)
        g_mean[g] = m; g_std[g] = A.std(axis=0); g_xyz[g] = x; g_vol[g] = v
        g_h1[g] = A[:F // 2].mean(axis=0); g_h2[g] = A[F // 2:].mean(axis=0)
        vol_t[:F] += A @ v; mom_t[:F] += A @ (v * x[:, up])
        ch = np.clip(np.searchsorted(edges, x[:, along], side="right") - 1, 0, a.chambers - 1)
        for c in range(a.chambers):
            chamber_t[:F, c] += A[:, ch == c] @ v[ch == c]
        bed = np.nonzero(m > a.bed_min)[0]
        take = rng.choice(bed, size=min(bed.size, int(np.ceil(3.0 * a.cells / n_ranks)) + 2), replace=False) if bed.size else bed
        cells.append(A[:, take]); cell_xyz.append(x[take])
        for k, (axn, val, _) in enumerate(planes):
            sel = np.nonzero(np.abs(x[:, axn] - val) <= 1.0 * np.cbrt(v))[0]      # cells within one cell width of the plane
            plane_series[k].append(A[:, sel]); plane_ids[k].append(g[sel])
        del A

    F = vol_t.size; t = np.arange(F) * dt
    cells = np.hstack(cells); cell_xyz = np.vstack(cell_xyz)
    if cells.shape[1] > a.cells:
        pick = rng.choice(cells.shape[1], size=a.cells, replace=False); cells = cells[:, pick]; cell_xyz = cell_xyz[pick]
    S = ["database              %s" % os.path.abspath(case),
         "frames                %d of %g s (%g s), every %d used; %d cells, %d ranks, phase '%s'"
         % (info["n_frames"], info["dt"], info["n_frames"] * info["dt"], a.every, n_glob, n_ranks, pname)]
    res = dict(case=case, dt=info["dt"], T=info["n_frames"] * info["dt"])

    # ---- recurrence matrix of the library, and from samples
    mats = {}
    mpath = os.path.join(case, "data", "rec", "matrix_0_0_0")
    if os.path.isfile(mpath):
        m = np.fromfile(mpath, sep=" "); n = int(round(np.sqrt(m.size)))
        if n * n == m.size:
            mats["library (all sampled cells)"] = m.reshape(n, n)[::a.every, ::a.every]
    if cells.shape[1] >= 2:
        mats["%d random bed cells" % cells.shape[1]] = l1_matrix(cells.astype(float))
    for k, (_, _, label) in enumerate(planes):
        ps = np.hstack(plane_series[k])
        if ps.shape[1] > 2000:
            ps = ps[:, rng.choice(ps.shape[1], size=2000, replace=False)]
        if ps.shape[1] >= 2:
            mats["plane %s" % label] = l1_matrix(ps.astype(float))
    nums = {name: matrix_numbers(m, dt) for name, m in mats.items()}
    ref = next(iter(mats))
    S.append("")
    S.append("recurrence (distance = sum over cells of |difference of volume fraction|; 'step' = two successive frames)")
    for name, q in nums.items():
        S.append("  %-28s decorrelation %.3g s, return %.3g s (minimum at %.2f of the plateau%s), partner %.2f steps, seam %.2f steps (%.2f of the plateau)"
                 % (name, q["t_dec"], q["t_ret"], q["depth"], ", no distinct period" if q["depth"] > 0.93 else "",
                    q["partner"], q["seam"], q["seam_plateau"]))
    q = nums[ref]
    res.update(t_dec=q["t_dec"], t_ret=q["t_ret"], depth=q["depth"], partner=q["partner"], seam=q["seam"])
    S.append("  database length             %.0f decorrelation times" % (res["T"] / q["t_dec"])
             + (", %.1f return times" % (res["T"] / q["t_ret"]) if q["depth"] <= 0.93 else ""))
    if q["loop"]:
        S.append("  smoothest loop              frames %d to %d (%d frames): jump of %.2f steps instead of %.2f at the seam"
                 % (q["loop"][0] * a.every, q["loop"][1] * a.every, (q["loop"][1] - q["loop"][0] + 1) * a.every, q["loop"][2], q["seam"]))
    iu = np.triu_indices(mats[ref].shape[0], 1)
    for name, m in mats.items():
        if name != ref and m.shape == mats[ref].shape:
            S.append("  %-28s correlation with the library matrix %.3f" % (name, np.corrcoef(m[iu], mats[ref][iu])[0, 1]))

    lagt = np.arange(mats[ref].shape[0]) * dt
    with open(os.path.join(out, "recurrence_lag.csv"), "w") as f:
        f.write("# time lag [s], mean distance / plateau: " + "; ".join(mats) + "\n")
        for i in range(lagt.size):
            f.write("%.6g," % lagt[i] + ",".join("%.5f" % (nums[k]["lag"][i] / nums[k]["plateau"]) for k in mats) + "\n")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                         "savefig.dpi": 180, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.4), constrained_layout=True)
    im = ax[0].imshow(mats[ref] / q["plateau"], origin="lower", cmap=SEQ, extent=[0, lagt[-1] + dt, 0, lagt[-1] + dt], vmin=0, vmax=1.3)
    style(ax[0], "time of frame i in s", "time of frame j in s", "Recurrence matrix, %s" % ref); ax[0].grid(False)
    cb = fig.colorbar(im, ax=ax[0], shrink=0.85); cb.set_label("distance / plateau", color=MUTED); cb.ax.tick_params(labelsize=7.5, colors=MUTED); cb.outline.set_visible(False)
    for i, (name, qq) in enumerate(nums.items()):
        ax[1].plot(lagt[: lagt.size * 3 // 4], (qq["lag"] / qq["plateau"])[: lagt.size * 3 // 4], color=SERIES[i % len(SERIES)], lw=1.7 if i == 0 else 1.0, label=name)
    ax[1].axhline(1.0, color="#b5b4ae", lw=0.7)
    ax[1].axvline(q["t_dec"], color=MUTED, lw=0.7, ls=(0, (2, 2)))
    if q["depth"] <= 0.93:
        ax[1].axvline(q["t_ret"], color=MUTED, lw=0.7, ls=(0, (5, 2)))
    ax[1].annotate("decorrelation %.3g s" % q["t_dec"], (q["t_dec"], 0.08), xytext=(4, 0), textcoords="offset points", fontsize=7.5, color=MUTED)
    if q["depth"] <= 0.93:
        ax[1].annotate("return %.3g s" % q["t_ret"], (q["t_ret"], 0.2), xytext=(-4, 0), textcoords="offset points", ha="right", fontsize=7.5, color=MUTED)
    style(ax[1], "time lag in s", "mean distance / plateau", "Distance against time lag"); ax[1].set_ylim(bottom=0); ax[1].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.savefig(os.path.join(out, "recurrence_matrix.png")); plt.close(fig)

    # ---- random cells
    ac, power = autocorrelation(cells.astype(float)) if cells.shape[1] else (np.zeros(F), np.zeros(F // 2 + 1))
    freq = np.fft.rfftfreq(F, dt)
    zero = np.nonzero(ac <= 0)[0]; kz = int(zero[0]) if zero.size else F
    t_int = float(ac[:kz].sum() * dt)
    f_peak = float(freq[1 + np.argmax(power[1:])]) if power.size > 2 and power[1:].max() > 0 else 0.0
    f_mean = float((freq[1:] * power[1:]).sum() / power[1:].sum()) if power[1:].sum() > 0 else 0.0
    S.append("")
    S.append("volume fraction at %d random cells with a time mean above %g" % (cells.shape[1], a.bed_min))
    S.append("  integral time scale         %.3g s  (the database holds about %.0f independent samples)" % (t_int, res["T"] / (2 * t_int) if t_int > 0 else 0))
    S.append("  spectrum                    peak at %.3g Hz, mean frequency %.3g Hz, lowest resolved %.3g Hz" % (f_peak, f_mean, freq[1] if freq.size > 1 else 0))
    if freq.size > 2 and f_peak <= freq[1] * 1.01:
        S.append("  note                        the peak is the lowest frequency the database resolves: its slowest motion is not covered by one period")
    res.update(t_int=t_int, f_peak=f_peak)
    np.savetxt(os.path.join(out, "cells.csv"), np.column_stack([t, cells]), delimiter=",", fmt="%.5g",
               header="time [s], volume fraction at the cells; centres: " + " ".join("(%.4g %.4g %.4g)" % tuple(c) for c in cell_xyz))
    fig, ax = plt.subplots(1, 3, figsize=(12.6, 3.4), constrained_layout=True)
    for i in range(min(5, cells.shape[1])):
        ax[0].plot(t, cells[:, i], color=SERIES[i], lw=0.8)
    style(ax[0], "time in the database in s", "volume fraction", "Five of the cells")
    ax[1].plot(t[: F // 2], ac[: F // 2], color=SERIES[0], lw=1.6); ax[1].axhline(0, color="#b5b4ae", lw=0.7)
    style(ax[1], "time lag in s", "autocorrelation", "Mean autocorrelation, integral time %.3g s" % t_int)
    if power.size > 2:
        ax[2].loglog(freq[1:], power[1:] / power[1:].max(), color=SERIES[0], lw=1.2)
    style(ax[2], "frequency in Hz", "power / peak", "Mean spectrum, peak at %.3g Hz" % f_peak)
    fig.savefig(os.path.join(out, "cells.png")); plt.close(fig)

    # ---- bed: height, volume per chamber, drift
    height = mom_t / np.maximum(vol_t, 1e-300)
    bedc = g_mean > a.bed_min
    fluct = np.sqrt(np.mean(g_std[bedc] ** 2)) if bedc.any() else 0.0
    halves = np.sqrt(np.mean((g_h1[bedc] - g_h2[bedc]) ** 2)) / fluct if fluct > 0 else float("nan")
    tr = lambda y: np.polyfit(t, y, 1)[0] * (t[-1] - t[0]) if F > 2 else 0.0
    S.append("")
    S.append("bed (phase '%s', %d cells with a time mean above %g)" % (pname, int(bedc.sum()), a.bed_min))
    S.append("  phase volume                %.6g m3, std %.3f %%, trend over the database %+.3f %%" % (vol_t.mean(), 100 * vol_t.std() / vol_t.mean(), 100 * tr(vol_t) / vol_t.mean()))
    S.append("  height of the centre        %.4g m along %s, std %.2f %%, trend %+.2f %%" % (height.mean(), a.up, 100 * height.std() / height.mean(), 100 * tr(height) / height.mean()))
    for c in range(a.chambers):
        y = chamber_t[:, c]
        if y.mean() > 0:
            S.append("  chamber %d (%s %.3g..%.3g)   %5.2f %% of the phase, std %.2f %% of its mean, trend %+.2f %%"
                     % (c, "xyz"[along], edges[c], edges[c + 1], 100 * y.mean() / vol_t.mean(), 100 * y.std() / y.mean(), 100 * tr(y) / y.mean()))
    S.append("  halves                      %.3f  (rms difference of the time mean of the two halves / rms fluctuation %.3f)" % (halves, fluct))
    res.update(halves=halves, vol_std=100 * vol_t.std() / vol_t.mean())
    np.savetxt(os.path.join(out, "bed.csv"), np.column_stack([t, vol_t, height, chamber_t]), delimiter=",", fmt="%.7g",
               header="time [s], phase volume [m3], height of the centre [m], phase volume per chamber [m3]")
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.3), constrained_layout=True)
    ax[0].plot(t, height, color=SERIES[0], lw=1.2); style(ax[0], "time in the database in s", "height of the centre in m", "Bed height")
    for c in range(a.chambers):
        if chamber_t[:, c].mean() > 0:
            ax[1].plot(t, 100 * chamber_t[:, c] / vol_t.mean(), color=SERIES[c % len(SERIES)], lw=1.1, label="chamber %d" % c)
    style(ax[1], "time in the database in s", "share of the phase in %", "Phase per chamber along %s" % "xyz"[along]); ax[1].legend(frameon=False, fontsize=7.5, ncol=2)
    fig.savefig(os.path.join(out, "bed.png")); plt.close(fig)

    # ---- planes
    fig, ax = plt.subplots(len(planes), 2, figsize=(9.6, 3.6 * len(planes)), constrained_layout=True, squeeze=False)
    for k, (axn, val, label) in enumerate(planes):
        ids = np.concatenate(plane_ids[k]) if plane_ids[k] else np.zeros(0, int)
        ids = ids[np.argsort(-np.abs(g_xyz[ids, axn] - val))]           # the cells nearest to the plane are drawn last
        h = [d for d in range(3) if d not in (axn, up)][0] if axn != up else [d for d in range(3) if d != up][0]
        v_ = up if axn != up else [d for d in range(3) if d != up][1]
        for j, (field, name) in enumerate(((g_mean, "time mean"), (g_std, "standard deviation"))):
            flat = ids.size < 3 or np.ptp(g_xyz[ids, h]) < 1e-9 or np.ptp(g_xyz[ids, v_]) < 1e-9
            if not flat:
                # one square of the cell width per cell (cells are drawn as cubes of their volume)
                half = 0.5 * np.cbrt(g_vol[ids]); cx, cy = g_xyz[ids, h], g_xyz[ids, v_]
                verts = np.stack([np.stack([cx - half, cy - half], 1), np.stack([cx + half, cy - half], 1),
                                  np.stack([cx + half, cy + half], 1), np.stack([cx - half, cy + half], 1)], axis=1)
                sc = PolyCollection(verts, array=field[ids], cmap=SEQ, linewidths=0)
                ax[k, j].add_collection(sc); ax[k, j].autoscale_view()
                cb = fig.colorbar(sc, ax=ax[k, j], shrink=0.85); cb.ax.tick_params(labelsize=7.5, colors=MUTED); cb.outline.set_visible(False)
                ax[k, j].set_aspect("equal")
            else:
                ax[k, j].text(0.5, 0.5, "the plane holds one line of cells or none", transform=ax[k, j].transAxes, ha="center", color=MUTED)
            style(ax[k, j], "%s in m" % "xyz"[h], "%s in m" % "xyz"[v_], "Plane %s: %s of the volume fraction" % (label, name)); ax[k, j].grid(False)
    fig.savefig(os.path.join(out, "planes.png")); plt.close(fig)

    # ---- shifts
    if not a.no_shifts:
        frames = np.unique(np.linspace(0, info["n_frames"] - 1, min(a.shift_frames, info["n_frames"])).astype(int))
        rows = []; hist = {}
        S.append("")
        S.append("shifts (%d frames spread over the database; cell width h = volume^(1/3))" % frames.size)
        for p in range(info["n_phases"]):
            acc = {int(k): [] for k in frames}
            for r in range(n_ranks):
                xyz, vol, host, gid = topo[r]; n = xyz.shape[0]
                for k, c0, c1, w, src in shift_frames(case, p, r, frames, n_ranks):
                    if c0.size == 0:
                        continue
                    x0 = np.empty((c0.size, 3))
                    for s in np.unique(src):
                        sel = src == s; x0[sel] = topo[int(s)][0][c0[sel]]
                    sw = np.bincount(c1, weights=w, minlength=n).astype(float); cnt = np.bincount(c1, minlength=n)
                    mean = np.stack([np.bincount(c1, weights=w * x0[:, d], minlength=n) for d in range(3)], axis=1).astype(float)
                    sq = np.bincount(c1, weights=w * (x0 ** 2).sum(axis=1), minlength=n)
                    has = sw > 0
                    mean[has] /= sw[has][:, None]
                    var = np.zeros(n); var[has] = np.maximum(sq[has] / sw[has] - (mean[has] ** 2).sum(axis=1), 0.0)
                    many = cnt > 1
                    var[many] *= cnt[many] / (cnt[many] - 1.0); var[~many] = 0.0
                    hcell = np.cbrt(vol)
                    disp = np.linalg.norm(xyz - mean, axis=1) / hcell
                    self_w = np.bincount(c1[(src == r) & (c0 == c1)], weights=w[(src == r) & (c0 == c1)], minlength=n)
                    acc[k].append((disp[has], var[has] / (2.0 * info["dt"]), self_w[has] / sw[has], w[src != r].sum(), w.sum(), np.sqrt(var[has]) / hcell[has]))
            dall, Dall = [], []
            for k in frames:
                if not acc[int(k)]:
                    continue
                d = np.concatenate([x[0] for x in acc[int(k)]]); D = np.concatenate([x[1] for x in acc[int(k)]])
                sf = np.concatenate([x[2] for x in acc[int(k)]]); sg = np.concatenate([x[5] for x in acc[int(k)]])
                other = sum(x[3] for x in acc[int(k)]) / max(sum(x[4] for x in acc[int(k)]), 1e-300)
                rows.append((p, k, np.median(d), np.percentile(d, 99), d.max(), 100 * sf.mean(), 100 * other, D.mean(), np.median(D), np.percentile(D, 99), sg.mean()))
                dall.append(d); Dall.append(D)
            if not dall:
                continue
            R = np.array([x[2:] for x in rows if x[0] == p]).mean(axis=0)
            name = info["phases"][p] if p < len(info["phases"]) else "phase %d" % p
            S.append("  %-10s displacement per frame: median %.2f h, 99 %% below %.2f h, largest %.1f h; from the cell itself %.1f %% of the weight, from other ranks %.1f %%"
                     % (name, R[0], R[1], R[2], R[3], R[4]))
            S.append("  %-10s D_C2C: mean %.3g m2/s, median %.3g, 99 %% below %.3g; spread of the sources %.2f h per frame; mixing length in 1 s: %.1f mm"
                     % ("", R[5], R[6], R[7], R[8], 1000 * np.sqrt(2 * R[5])))
            hist[name] = (np.concatenate(dall), np.concatenate(Dall))
            if p == phase:
                res.update(disp=R[0], disp99=R[1], D=R[5], self=R[3])
        np.savetxt(os.path.join(out, "shifts.csv"), np.array(rows), delimiter=",", fmt="%.6g",
                   header="phase, frame, displacement median [h], 99 % [h], max [h], weight from the cell itself [%], from other ranks [%], D_C2C mean [m2/s], median, 99 %, spread of the sources [h]")
        fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.3), constrained_layout=True)
        for i, (name, (d, D)) in enumerate(hist.items()):
            ax[0].hist(d, bins=np.linspace(0, max(np.percentile(d, 99.5), 1e-6), 60), histtype="step", color=SERIES[i], lw=1.3, density=True, label=name)
            ax[1].hist(D[D > 0], bins=np.logspace(np.log10(max(D[D > 0].min(), 1e-9)), np.log10(D.max()), 60) if (D > 0).any() else 10, histtype="step", color=SERIES[i], lw=1.3, label=name)
        style(ax[0], "displacement per frame in cell widths", "density", "How far the content of a cell comes from"); ax[0].legend(frameon=False, fontsize=7.5)
        ax[1].set_xscale("log"); style(ax[1], "D_C2C in m2/s", "cells", "Numerical diffusivity of the shift (frame interval %g s)" % info["dt"]); ax[1].legend(frameon=False, fontsize=7.5)
        fig.savefig(os.path.join(out, "shifts.png")); plt.close(fig)

    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("\n".join(S) + "\n")
    print("\n".join(S)); print("written %s (%.0f s)\n" % (out, time.time() - t0)); sys.stdout.flush()
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default=["."], nargs="+", help="case folder(s) holding data/ (default: here)")
    p.add_argument("--phase", default=None, help="phase name or number (default: the last phase)")
    p.add_argument("--cells", type=int, default=200, help="number of random bed cells")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--bed-min", type=float, default=0.05, help="a bed cell has a time-mean volume fraction above this")
    p.add_argument("--plane", nargs="+", default=None, help="planes, e.g. x=0 z=0.375 (default: two mid planes)")
    p.add_argument("--up", default="y", choices="xyz", help="axis pointing up")
    p.add_argument("--along", default=None, choices="xyz", help="axis of the chambers (default: the longer horizontal one)")
    p.add_argument("--chambers", type=int, default=4)
    p.add_argument("--every", type=int, default=1, help="use every n-th frame for the volume-fraction statistics")
    p.add_argument("--shift-frames", type=int, default=10, help="number of frames for the shift statistics")
    p.add_argument("--no-shifts", action="store_true")
    a = p.parse_args()

    rows = [analyse(c, a) for c in a.case]
    if len(rows) > 1:
        print("%-44s %7s %7s %8s %8s %8s %7s %7s %7s %9s %10s" % ("database", "dt ms", "T s", "decorr s", "return s", "partner", "seam", "halves", "disp h", "self %", "D_C2C m2/s"))
        for r in rows:
            print("%-44s %7.3g %7.3g %8.3g %8.3g %8.2f %7.2f %7.3f %7.2f %9.1f %10.3g"
                  % (os.path.basename(os.path.abspath(r["case"]))[-44:], 1000 * r["dt"], r["T"], r["t_dec"], r["t_ret"], r["partner"], r["seam"], r["halves"],
                     r.get("disp", float("nan")), r.get("self", float("nan")), r.get("D", float("nan"))))


if __name__ == "__main__":
    main()

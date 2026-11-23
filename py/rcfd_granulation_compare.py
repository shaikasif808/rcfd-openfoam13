#!/usr/bin/env python3
"""Several granulation runs in one picture, e.g. the same batch at different frame intervals.

    python3 rcfd_granulation_compare.py dt05 dt10 dt15 dt20 dt25
    python3 rcfd_granulation_compare.py dt05 dt25 --measured samples.csv
    python3 rcfd_granulation_compare.py "5 ms, shifting=batch_05" "25 ms, shifting=batch_25" --out compare.png

A run is a folder with rCFD_granulation_monitor.out, or a batch folder with stage_1, stage_2, ...
(the stages are joined). Give "label=folder" to name a run; without a label the frame interval of
its header is used ("5 ms"), else the folder name. Runs are drawn in the order given.

Reads per run: rCFD_granulation_monitor.out, rCFD_granulation_chambers.out, log.rcfd_run (steps and
run time) and rCFD/rCFD_user_run.h (frame interval, growth acceleration). Writes the picture
(default rCFD_granulation_compare.png) and a table with the same name and the ending .csv into the
current folder. Nothing in the run folders is changed.

Panels: mass median, span, variance, difference between the chambers; then against the frame
interval (or run number): variance added per minute, and run time per minute of process time.

--measured: text file with  time in min, x50 in mm, span  (as for rcfd_granulation_plot.py).
A run made with _growth_acceleration F > 1 has its added variance divided by F.

Needs numpy and matplotlib.
"""
import argparse, os, re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex

INK, MUTED, GRID, ACCENT = "#0b0b0b", "#52514e", "#e4e3df", "#eb6834"
RAMP = LinearSegmentedColormap.from_list("ramp", ["#86b6ef", "#2a78d6", "#0d366b"])
MONITOR, CHAMBERS = "rCFD_granulation_monitor.out", "rCFD_granulation_chambers.out"


def table(path):
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(v) for v in line.replace(",", " ").split()])
            except ValueError:
                pass
    if not rows:
        return np.zeros((0, 0))
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def stage_folders(root):
    if os.path.isfile(os.path.join(root, MONITOR)):
        return [root]
    st = [d for d in os.listdir(root) if d.lower().startswith("stage") and os.path.isfile(os.path.join(root, d, MONITOR))]
    key = lambda s: (int(re.search(r"(\d+)\s*$", s).group(1)) if re.search(r"(\d+)\s*$", s) else 10 ** 9, s)
    return [os.path.join(root, d) for d in sorted(st, key=key)]


def define(text, name):
    m = re.search(r"^[ \t]*#define[ \t]+%s[ \t]+([-+0-9.eE]+)" % re.escape(name), text, re.M)
    return float(m.group(1)) if m else None


def load(root):
    stages = stage_folders(root)
    if not stages:
        raise SystemExit("%s: no %s here or in stage_* folders" % (root, MONITOR))
    r = dict(root=root, stages=stages, dts=[], acc=1.0, steps=0, wall=0.0, switches=[])
    mons, chs = [], []
    for s in stages:
        m = table(os.path.join(s, MONITOR))
        if m.size == 0 or m.shape[1] < 7:
            continue
        if mons:
            r["switches"].append(m[0, 0] / 60.0)
        mons.append(m)
        if os.path.isfile(os.path.join(s, CHAMBERS)):
            c = table(os.path.join(s, CHAMBERS))
            if c.size and (c.shape[1] - 1) % 3 == 0:
                chs.append(c)
        for h in (os.path.join(s, "rCFD", "rCFD_user_run.h"), os.path.join(s, "rCFD_user_run.h")):
            if os.path.isfile(h):
                text = open(h, errors="replace").read()
                dt, acc = define(text, "_global_time_step"), define(text, "_growth_acceleration")
                if dt:
                    r["dts"].append(dt)
                if acc:
                    r["acc"] = acc
                break
        log = os.path.join(s, "log.rcfd_run")
        if os.path.isfile(log):
            for line in open(log, errors="replace"):
                m_ = re.search(r"run finished: (\d+) steps.*?run time ([0-9.eE+-]+) s", line)
                if m_:
                    r["steps"] += int(m_.group(1)); r["wall"] += float(m_.group(2))
    if not mons:
        raise SystemExit("%s: the monitor files are empty" % root)
    w = min(m.shape[1] for m in mons)
    m = np.vstack([q[:, :w] for q in mons])
    F = max(r["acc"], 1.0)
    t = m[:, 0] / 60.0
    mean, sigma = m[:, 4], m[:, 6]
    var = sigma[0] ** 2 + np.maximum(sigma ** 2 - sigma[0] ** 2, 0.0) / F
    x50 = mean + 3.0 * var / mean
    r.update(t=t, mass=m[:, 2], mean=mean, var=var, x50=x50, span=2.563 * np.sqrt(var) / x50,
             limited=m[:, 18] if w > 18 else np.zeros(len(m)), deposited=m[:, 9] if w > 9 else np.zeros(len(m)),
             lost=m[:, 10] if w > 10 else np.zeros(len(m)), restored=m[:, 14] if w > 14 else np.zeros(len(m)))
    r["var_rate"] = float(np.polyfit(t, var, 1)[0]) if t[-1] > t[0] else float("nan")
    if chs and len({c.shape[1] for c in chs}) == 1:
        c = np.vstack(chs); n = (c.shape[1] - 1) // 3
        d = np.stack([c[:, 2 + 3 * k] for k in range(n)], axis=1)
        r["tc"], r["spread"] = c[:, 0] / 60.0, 1000.0 * (d.max(axis=1) - d.min(axis=1))
    dts = sorted(set(round(x, 9) for x in r["dts"]))
    r["dt"] = dts[0] if len(dts) == 1 else None
    r["dt_text"] = " / ".join("%g" % (1000 * x) for x in dts) + " ms" if dts else ""
    return r


def style(ax, xlabel, ylabel, title):
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#b5b4ae"); ax.spines[s].set_linewidth(0.6)
    ax.tick_params(colors=MUTED, length=2.5, width=0.5, labelsize=7.5)
    ax.set_xlabel(xlabel, color=MUTED); ax.set_ylabel(ylabel, color=MUTED); ax.set_title(title, loc="left", color=INK)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="+", help="run or batch folders, optionally as label=folder")
    p.add_argument("--measured", default=None, help="plant samples: time [min], x50 [mm], span")
    p.add_argument("--out", default="rCFD_granulation_compare.png")
    a = p.parse_args()

    runs = []
    for item in a.runs:
        label, folder = item.split("=", 1) if "=" in item and not os.path.isdir(item) else (None, item)
        r = load(folder)
        r["label"] = label or r["dt_text"] or os.path.basename(os.path.normpath(folder))
        runs.append(r)
    if len({r["label"] for r in runs}) < len(runs):                 # same interval twice: add the folder
        for r in runs:
            r["label"] += " (%s)" % os.path.basename(os.path.normpath(r["root"]))
    n = len(runs)
    cols = [to_hex(RAMP(k / max(n - 1, 1))) for k in range(n)]
    meas = table(a.measured) if a.measured else None
    rate_m = None
    if meas is not None and meas.shape[1] > 2 and len(meas) > 2:
        rate_m = float(np.polyfit(meas[:, 0], (meas[:, 2] * meas[:, 1] / 2.563) ** 2, 1)[0])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                         "savefig.dpi": 200, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(2, 3, figsize=(13.2, 6.6), constrained_layout=True)
    for r, c in zip(runs, cols):
        ax[0, 0].plot(r["t"], r["x50"], color=c, lw=1.6, label=r["label"])
        ax[0, 1].plot(r["t"], r["span"], color=c, lw=1.6, label=r["label"])
        ax[0, 2].plot(r["t"], 1e3 * r["var"], color=c, lw=1.6, label=r["label"])
        if "spread" in r:
            k = max(1, len(r["spread"]) // 30)                       # running mean: the replay loop shows as a sawtooth
            sm = np.convolve(r["spread"], np.ones(k) / k, mode="valid")
            ax[1, 0].plot(r["tc"][k - 1:], sm, color=c, lw=1.4, label=r["label"])
    if meas is not None:
        mk = dict(ls="none", marker="o", ms=4.5, mfc="white", mec=INK, mew=0.9, label="measured")
        ax[0, 0].plot(meas[:, 0], meas[:, 1], **mk)
        if meas.shape[1] > 2:
            ax[0, 1].plot(meas[:, 0], meas[:, 2], **mk)
            ax[0, 2].plot(meas[:, 0], 1e3 * (meas[:, 2] * meas[:, 1] / 2.563) ** 2, **mk)
    for axis in (ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0]):
        t_end = max(r["t"][-1] for r in runs)
        axis.set_xlim(-0.015 * t_end, t_end * 1.02)
    style(ax[0, 0], "process time in min", "mass median x50 in mm", "Size"); ax[0, 0].legend(frameon=False, fontsize=7.5)
    style(ax[0, 1], "process time in min", "span (x90 - x10) / x50", "Width")
    ax[0, 1].set_ylim(0, 1.15 * max([r["span"].max() for r in runs] + ([meas[:, 2].max()] if meas is not None and meas.shape[1] > 2 else [])))
    style(ax[0, 2], "process time in min", "variance in 0.001 mm2", "Variance of the size")
    style(ax[1, 0], "process time in min", "largest - smallest chamber mean in µm", "Difference between the chambers (running mean)")
    ax[1, 0].set_ylim(bottom=0)

    # against the frame interval where every run has one, else against the run number
    by_dt = all(r["dt"] for r in runs) and len({r["dt"] for r in runs}) == n
    x = np.array([1000 * r["dt"] for r in runs]) if by_dt else np.arange(1, n + 1)
    xl = "frame interval in ms" if by_dt else "run"
    rate = np.array([r["var_rate"] for r in runs])
    ax[1, 1].plot(x, 1e5 * rate, color=MUTED, lw=0.8, zorder=1)
    ax[1, 1].scatter(x, 1e5 * rate, c=cols, s=42, zorder=2, edgecolors="white", linewidths=0.8)
    if rate_m is not None:
        ax[1, 1].axhline(1e5 * rate_m, color=ACCENT, lw=1.0, ls=(0, (4, 2)))
        ax[1, 1].annotate("samples", (x.max(), 1e5 * rate_m), xytext=(0, 3), textcoords="offset points", ha="right", fontsize=7.5, color=ACCENT)
    style(ax[1, 1], xl, "variance added in 1e-5 mm2 per min", "Widening of the distribution")
    top = np.nanmax(np.append(rate, rate_m if rate_m is not None else np.nan))
    ax[1, 1].set_ylim(0, 1e5 * top * 1.25 if np.isfinite(top) and top > 0 else 1)
    cost = np.array([r["wall"] / (r["t"][-1] - r["t"][0]) if r["wall"] > 0 and r["t"][-1] > r["t"][0] else np.nan for r in runs])
    if np.isfinite(cost).any():
        ax[1, 2].plot(x, cost, color=MUTED, lw=0.8, zorder=1)
        ax[1, 2].scatter(x, cost, c=cols, s=42, zorder=2, edgecolors="white", linewidths=0.8)
        for xi, ci, r in zip(x, cost, runs):
            if np.isfinite(ci):
                ax[1, 2].annotate("%.0f ms/step" % (1000 * r["wall"] / max(r["steps"], 1)), (xi, ci), xytext=(0, 6),
                                  textcoords="offset points", ha="center", fontsize=7, color=MUTED)
    else:
        ax[1, 2].text(0.5, 0.5, "no log.rcfd_run in the run folders", transform=ax[1, 2].transAxes, ha="center", color=MUTED)
    style(ax[1, 2], xl, "run time in s per min of process time", "Cost")
    ax[1, 2].set_ylim(0, np.nanmax(cost) * 1.25 if np.isfinite(cost).any() else 1)
    if not by_dt:
        for axis in (ax[1, 1], ax[1, 2]):
            axis.set_xticks(x); axis.set_xticklabels([r["label"] for r in runs], rotation=20, ha="right")
    fig.savefig(a.out); plt.close(fig)

    head = ["run", "folder", "frame interval [ms]", "stages", "acceleration", "process time [min]", "bed mass end [kg]",
            "x50 end [mm]", "span start", "span end", "variance added [mm2/min]", "chamber difference end [um]",
            "chamber difference mean [um]", "restored by correction [% of growth]", "limited [% of deposit]",
            "not deposited [% of deposit]", "steps", "run time [s]", "run time per step [ms]", "run time per process min [s]"]
    rows = []
    for r, c in zip(runs, cost):
        g = r["mean"][-1] - r["mean"][0]
        dep = r["deposited"][-1]
        rows.append([r["label"], r["root"], r["dt_text"].replace(" ms", ""), len(r["stages"]), "%g" % r["acc"], "%.2f" % r["t"][-1],
                     "%.3f" % r["mass"][-1], "%.4f" % r["x50"][-1], "%.4f" % r["span"][0], "%.4f" % r["span"][-1],
                     "%.4g" % r["var_rate"], "%.2f" % r["spread"][-1] if "spread" in r else "",
                     "%.2f" % r["spread"].mean() if "spread" in r else "", "%.1f" % (100 * r["restored"][-1] / g) if g > 0 else "",
                     "%.3g" % (100 * r["limited"][-1] / dep) if dep > 0 else "", "%.3g" % (100 * r["lost"][-1] / dep) if dep > 0 else "",
                     r["steps"] or "", "%.0f" % r["wall"] if r["wall"] else "",
                     "%.1f" % (1000 * r["wall"] / r["steps"]) if r["steps"] else "", "%.1f" % c if np.isfinite(c) else ""])
    csv = os.path.splitext(a.out)[0] + ".csv"
    with open(csv, "w") as f:
        f.write(",".join(head) + "\n")
        for row in rows:
            f.write(",".join('"%s"' % v if "," in str(v) else str(v) for v in row) + "\n")

    print("%-22s %7s %9s %9s %8s %8s %12s %10s %9s %9s %10s" % ("run", "t min", "mass kg", "x50 mm", "span 0", "span", "var/min mm2", "chamb. um", "restor. %", "ms/step", "s/proc.min"))
    for row in rows:
        print("%-22s %7s %9s %9s %8s %8s %12s %10s %9s %9s %10s" % (str(row[0])[:22], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[13], row[18], row[19]))
    if rate_m is not None:
        print("samples: variance added %.4g mm2/min" % rate_m)
    print("written %s and %s" % (a.out, csv))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Pictures of a granulation run with the three moment fields (rCFD_user_run_granulation.h).

    python3 rcfd_granulation_plot.py                         # in the case folder
    python3 rcfd_granulation_plot.py --measured samples.csv  # with plant samples
    python3 rcfd_granulation_plot.py --acceleration 10       # run made with _growth_acceleration 10

Reads rCFD_granulation_monitor.out and rCFD_granulation_chambers.out, writes rCFD_granulation.png.

--measured: text file with the columns  time in min, x50 in mm, span (x90 - x10) / x50  (comma or
blank separated, lines starting with # are skipped). Drawn as points.

--acceleration F: the spread that growth adds is too large by about F (variance) in an accelerated
run. With this option the added variance is divided by F before spread and span are drawn; the
uncorrected span is drawn as a thin line.

Printed as well: the variance that growth adds per minute (slope of a straight line through
sigma^2 against time), for the run and for the samples. This is the number to compare when the
runs are shorter than the batch.

Needs numpy and matplotlib. Nothing in the case is changed.
"""
import argparse, os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def table(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(v) for v in line.replace(",", " ").split()])
            except ValueError:
                pass
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


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
    p.add_argument("--case", default=".", help="case folder (default: here)")
    p.add_argument("--measured", default=None, help="plant samples: time [min], x50 [mm], span")
    p.add_argument("--acceleration", type=float, default=1.0, help="_growth_acceleration of the run")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    m = table(os.path.join(a.case, "rCFD_granulation_monitor.out"))
    if len(m) < 2:
        raise SystemExit("rCFD_granulation_monitor.out holds fewer than two lines")
    t = m[:, 0] / 60.0
    mass_f, mass_b, mean, d30, sigma = m[:, 2], m[:, 3], m[:, 4], m[:, 5], m[:, 6]
    restored = m[:, 14] if m.shape[1] > 14 else np.zeros(len(m))

    # spread without the surplus of an accelerated run
    var_added = np.maximum(sigma ** 2 - sigma[0] ** 2, 0.0)
    sigma_c = np.sqrt(sigma[0] ** 2 + var_added / max(a.acceleration, 1.0))
    x50 = mean + 3.0 * sigma_c ** 2 / mean                  # mass median of a narrow, near-normal distribution
    span = 2.563 * sigma_c / x50
    span_raw = 2.563 * sigma / (mean + 3.0 * sigma ** 2 / mean)

    meas = table(a.measured) if a.measured else None

    # variance added per minute: run (without the surplus of the acceleration) and samples
    rate = np.polyfit(t, sigma_c ** 2, 1)[0] if t[-1] > t[0] else 0.0
    text = "added variance: run %.3g mm2/min" % rate
    if meas is not None and meas.shape[1] > 2 and len(meas) > 2:
        s_meas = meas[:, 2] * meas[:, 1] / 2.563
        rate_m = np.polyfit(meas[:, 0], s_meas ** 2, 1)[0]
        text += ", samples %.3g mm2/min (run / samples = %.2f)" % (rate_m, rate / rate_m if rate_m != 0 else float("nan"))
    print(text)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
                         "savefig.dpi": 200, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(2, 2, figsize=(9.6, 6.4), constrained_layout=True)

    ax[0, 0].plot(t, x50, color=SERIES[0], lw=1.8, label="rCFD, mass median")
    ax[0, 0].plot(t, d30, color=MUTED, lw=0.9, ls=(0, (4, 2)), label="volume mean (mass balance)")
    if meas is not None:
        ax[0, 0].plot(meas[:, 0], meas[:, 1], ls="none", marker="o", ms=4.5, mfc="white", mec=INK, mew=0.9, label="measured x50")
    style(ax[0, 0], "process time in min", "diameter in mm", "Size")
    ax[0, 0].legend(frameon=False, fontsize=7.5)

    ax[0, 1].plot(t, span, color=SERIES[0], lw=1.8, label="rCFD" + (", added variance / %g" % a.acceleration if a.acceleration > 1 else ""))
    if a.acceleration > 1:
        ax[0, 1].plot(t, span_raw, color=SERIES[0], lw=0.8, ls=(0, (2, 2)), label="rCFD as computed (accelerated)")
    if meas is not None and meas.shape[1] > 2:
        ax[0, 1].plot(meas[:, 0], meas[:, 2], ls="none", marker="o", ms=4.5, mfc="white", mec=INK, mew=0.9, label="measured")
    style(ax[0, 1], "process time in min", "span (x90 - x10) / x50", "Width")
    ax[0, 1].set_ylim(bottom=0); ax[0, 1].legend(frameon=False, fontsize=7.5)

    ax[1, 0].plot(t, mass_f, color=SERIES[0], lw=1.8, label="from the fields")
    ax[1, 0].plot(t, mass_b, color=MUTED, lw=0.9, ls=(0, (4, 2)), label="start + sprayed solid on the particles")
    style(ax[1, 0], "process time in min", "bed mass in kg", "Bed mass")
    ax[1, 0].legend(frameon=False, fontsize=7.5)

    cpath = os.path.join(a.case, "rCFD_granulation_chambers.out")
    if os.path.isfile(cpath):
        c = table(cpath); n = (c.shape[1] - 1) // 3
        for k in range(n):
            ax[1, 1].plot(c[:, 0] / 60.0, c[:, 2 + 3 * k], color=SERIES[k % len(SERIES)], lw=1.4, label="chamber %d" % k)
        ax[1, 1].legend(frameon=False, fontsize=7.5, ncol=2)
    style(ax[1, 1], "process time in min", "number mean diameter in mm", "Chambers along the apparatus")

    growth = mean[-1] - mean[0]
    fig.suptitle("%s\nbed %.2f -> %.2f kg, x50 %.3f -> %.3f mm, span %.3f -> %.3f; balance correction restored %.0f %% of the growth; %s"
                 % (os.path.abspath(a.case), mass_f[0], mass_f[-1], x50[0], x50[-1], span[0], span[-1],
                    100.0 * restored[-1] / growth if growth > 0 else 0.0, text), fontsize=8, x=0.01, ha="left")
    out = a.out or os.path.join(a.case, "rCFD_granulation.png")
    fig.savefig(out); plt.close(fig)
    print("written %s" % out)


if __name__ == "__main__":
    main()

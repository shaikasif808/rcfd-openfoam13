#!/usr/bin/env python3
"""Pictures of a running or finished rCFD run.

    python3 rcfd_monitor.py                      # in the case folder: writes rCFD_monitor.png once
    python3 rcfd_monitor.py --every 60           # writes it again every 60 s until stopped (Ctrl-C)
    python3 rcfd_monitor.py --shifts rCFD_c2c_step1_node1.txt     # picture of one shift file

rCFD_monitor.png, from the files the run writes:
    rCFD_mean_height_monitor.out   mean height and mass of every size class (example header)
    log.rcfd_run                   the [c2c] lines of the first steps: where the shifts take
                                   their content from, per processor

Needs numpy and matplotlib. Nothing is changed in the case; only the png files are written.
"""
import argparse, os, re, sys, time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BINS = ["0-5 mm", "5-10 mm", "10-20 mm", "20-50 mm", "> 50 mm"]
BIN_COLOURS = ["#2a78d6", "#1baf86", "#eda100", "#e25c4b", "#7a52c7"]
CLASS_COLOURS = ["#2a78d6", "#e25c4b", "#1baf86", "#eda100", "#7a52c7", "#555555"]


def read_table(path):
    rows = []
    if not os.path.isfile(path):
        return np.zeros((0, 0))
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            try:
                rows.append([float(v) for v in line.replace(",", " ").split()])
            except ValueError:
                continue
    if not rows:
        return np.zeros((0, 0))
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def read_c2c(path, step):
    """[c2c] lines of one step: per node the solids mass and, per source class, the share of the
    weight and its split over the distance bins."""
    head = re.compile(r"\[c2c\] step\s+(\d+) frame\s+(\d+) node\s+(\d+): solids\s+([-+.\deE]+) kg")
    cls = re.compile(r"\[c2c\]\s+node\s+(\d+) from (cells of this node|corona cells|other nodes)\s*:\s*(\d+) shifts,"
                     r"\s*([-+.\d]+) % of the weight.*?mm:\s*([-+.\d]+)\s+([-+.\d]+)\s+([-+.\d]+)\s+([-+.\d]+)\s+([-+.\d]+) %")
    nodes, current = {}, None
    if not os.path.isfile(path):
        return nodes
    with open(path, errors="replace") as f:
        for line in f:
            if "[c2c]" not in line:
                continue
            m = head.search(line)
            if m:
                current = int(m.group(1))
                if current == step:
                    nodes[int(m.group(3))] = {"mass": float(m.group(4)), "share": {}, "bins": np.zeros(5)}
                continue
            m = cls.search(line)
            if m and current == step and int(m.group(1)) in nodes:
                d = nodes[int(m.group(1))]
                d["share"][m.group(2)] = float(m.group(4))
                d["bins"] += np.array([float(m.group(i)) for i in range(5, 10)])
    return nodes


def last_step(path):
    step = None
    if not os.path.isfile(path):
        return None
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("[stage] step"):
                try:
                    step = int(line.split()[2])
                except (IndexError, ValueError):
                    pass
    return step


def monitor(run, log, out, c2c_step):
    h = read_table(os.path.join(run, "rCFD_mean_height_monitor.out"))
    nodes = read_c2c(os.path.join(run, log), c2c_step)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    if h.size and h.shape[1] >= 3:
        n = (h.shape[1] - 1) // 2
        t, height, mass = h[:, 0], h[:, 1:1 + n], h[:, 1 + n:1 + 2 * n]
        for i in range(n):
            c = CLASS_COLOURS[i % len(CLASS_COLOURS)]
            ax[0, 0].plot(t, 1e3 * height[:, i], color=c, label="class %d" % i)
            if i > 0:
                ax[0, 1].plot(t, 1e3 * (height[:, i] - height[:, 0]), color=c, label="class %d minus class 0" % i)
            if mass[0, i] > 0:
                ax[1, 0].plot(t, 100 * (mass[:, i] / mass[0, i] - 1), color=c, label="class %d (%.3g kg)" % (i, mass[0, i]))
        ax[0, 0].set_title("Mean height of every size class")
        ax[0, 0].set_ylabel("mean height [mm]")
        ax[0, 1].set_title("Segregation: difference of the mean heights")
        ax[0, 1].set_ylabel("difference [mm]")
        ax[0, 1].axhline(0, color="#999999", lw=0.8)
        ax[1, 0].set_title("Mass of every size class against the start")
        ax[1, 0].set_ylabel("change of the mass [%]")
        ax[1, 0].axhline(0, color="#999999", lw=0.8)
        for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
            a.set_xlabel("time [s]")
            a.grid(alpha=0.3)
            if a.get_legend_handles_labels()[0]:
                a.legend(fontsize=8)
        state = "t = %.3f s" % t[-1]
    else:
        for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
            a.text(0.5, 0.5, "rCFD_mean_height_monitor.out\nnot found or empty", ha="center", va="center", transform=a.transAxes)
            a.set_axis_off()
        state = "no monitor data"

    a = ax[1, 1]
    if nodes:
        ids = sorted(nodes)
        bottom = np.zeros(len(ids))
        for b in range(5):
            v = np.array([nodes[i]["bins"][b] for i in ids])
            a.bar(ids, v, bottom=bottom, color=BIN_COLOURS[b], label=BINS[b], width=0.8)
            bottom += v
        other = np.array([nodes[i]["share"].get("corona cells", 0) + nodes[i]["share"].get("other nodes", 0) for i in ids])
        a.plot(ids, other, "k.", ms=7, label="from other processors")
        a.set_ylim(0, 135)
        a.set_yticks(range(0, 101, 20))
        a.set_xticks(ids if len(ids) <= 16 else ids[::max(1, len(ids) // 16)])
        a.set_xlabel("processor")
        a.set_ylabel("share of the content of the cells [%]")
        a.set_title("Step %d: distance between a cell and where its content comes from" % c2c_step)
        a.legend(fontsize=8, ncol=3, loc="upper center")
        a.grid(alpha=0.3, axis="y")
    else:
        a.text(0.5, 0.5, "no [c2c] lines of step %d in %s" % (c2c_step, log), ha="center", va="center", transform=a.transAxes)
        a.set_axis_off()

    s = last_step(os.path.join(run, log))
    fig.suptitle("%s    %s%s    written %s" % (os.path.abspath(run), state, "" if s is None else ", step %d" % s,
                                              time.strftime("%Y-%m-%d %H:%M:%S")), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    tmp = out + ".tmp.png"
    fig.savefig(tmp, dpi=110)
    plt.close(fig)
    os.replace(tmp, out)


def shifts(path, out):
    """One shift file of the library diagnostics (rCFD_c2c_step<S>_node<N>.txt): one line per shift,
    destination cell c1 and source cell c0 with positions, class of the source and weight w."""
    d = read_table(path)
    if d.size == 0 or d.shape[1] < 17:
        print("%s: no shift lines" % path)
        return
    with open(path) as f:
        title = f.readline().lstrip("# ").strip()
    c1, x1, cls, x0, w, wsum = d[:, 0].astype(int), d[:, 1:4], d[:, 7].astype(int), d[:, 10:13], d[:, 15], d[:, 16]
    dist = np.linalg.norm(x1 - x0, axis=1)
    ws = np.maximum(wsum, 1e-30)
    share = w / ws                                      # share of the content of the destination cell

    cells, inv = np.unique(c1, return_inverse=True)
    pos = np.zeros((len(cells), 3)); pos[inv] = x1
    mean_dist = np.bincount(inv, share * dist, len(cells))
    outside = np.bincount(inv, share * (cls > 0), len(cells))

    ext = pos.max(0) - pos.min(0)
    a0, a1 = sorted(np.argsort(ext)[-2:])
    name = "xyz"

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    for a, v, lab in ((ax[0], 1e3 * mean_dist, "mean distance to the source cells [mm]"),
                      (ax[1], 100 * outside, "share of the content from other processors [%]")):
        sc = a.scatter(pos[:, a0], pos[:, a1], c=v, s=6, cmap="viridis")
        fig.colorbar(sc, ax=a, label=lab)
        a.set_xlabel("%s [m]" % name[a0]); a.set_ylabel("%s [m]" % name[a1]); a.set_aspect("equal", adjustable="datalim")
    ax[0].set_title("Cells of this processor: how far their content travels")
    ax[1].set_title("Cells that take content across a processor boundary")

    edges = [0, 0.005, 0.010, 0.020, 0.050, np.inf]
    total = share.sum()
    names = ["cells of this processor", "corona cells", "other processors", "other processors"]
    bottom = np.zeros(5)
    for k, col in ((0, "#2a78d6"), (1, "#eda100"), (2, "#e25c4b")):
        m = (cls == k) if k < 2 else (cls >= 2)
        v = np.array([share[m & (dist >= edges[b]) & (dist < edges[b + 1])].sum() for b in range(5)]) / total * 100
        ax[2].bar(range(5), v, bottom=bottom, color=col, label=names[k])
        bottom += v
    ax[2].set_xticks(range(5)); ax[2].set_xticklabels(BINS)
    ax[2].set_ylabel("share of the content [%]")
    ax[2].set_title("Where the content comes from, by distance")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3, axis="y")

    fig.suptitle("%s    (%d shifts into %d cells)" % (title, len(d), len(cells)), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("written %s" % out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", default=".", help="case folder of the rCFD run (default: here)")
    p.add_argument("--log", default="log.rcfd_run", help="log file of rcfd_run in the case folder")
    p.add_argument("--out", default=None, help="name of the picture (default rCFD_monitor.png in the case folder)")
    p.add_argument("--every", type=float, default=0, help="write the picture again every so many seconds")
    p.add_argument("--c2c-step", type=int, default=0, help="step of the [c2c] lines shown (default 0)")
    p.add_argument("--shifts", default=None, help="picture of one file rCFD_c2c_step<S>_node<N>.txt instead")
    a = p.parse_args()

    if a.shifts:
        shifts(a.shifts, a.out or os.path.splitext(a.shifts)[0] + ".png")
        return

    out = a.out or os.path.join(a.run, "rCFD_monitor.png")
    while True:
        try:
            monitor(a.run, a.log, out, a.c2c_step)
        except Exception as e:                          # a file that is being written at this moment
            print("rcfd_monitor: %s (trying again)" % e, file=sys.stderr)
        if a.every <= 0:
            print("written %s" % out)
            break
        time.sleep(a.every)


if __name__ == "__main__":
    main()

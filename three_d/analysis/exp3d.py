"""Tables, figures, and numbers of the experiments in R^3 (appendix Experiments in Three Dimensions of the paper;
Section 5.3 of the arXiv version).
Run from the repository root:  python3 three_d/analysis/exp3d.py          (PARTIAL=1: allow missing runs, for tests while a batch runs)
One 3D analogue of every plot of the planar experiments, on the same kind of data. Instances: seed 1 of the six
generators and one sample of each of the four scans, for every n (one seed; the seeds 2 and 3 of the
first run stay in three_d/data/ but are not used). Inputs (three_d/run3d.py, three_d/data/, no header lines; the folder in DATA3D if set):
  exp3d_sep.csv  sweep, CK box, CK exact, one process each (times); mem3d.csv: inst,algo,peakGB,wallS of these runs
                 and of lb_swap3; exp3d.csv: the same algorithms in one process (first run; checked to agree)
  lb3d_sep.csv   inst,gen,n,seed,check,<output line of lb_swap3> (LB_0, LB_1, LB and their times); lb3d.csv (first run)
  opt3d.csv      exact OPT (n <= 1000);  alpha3d.csv: the largest sets of pairwise incompatible pairs (n <= 500)
  ls3d.csv       reshape from CK exact and from the sweep (traces <start>_<inst>.csv[.gz]); lns3d.csv: LNS on the samples
                 (traces lns_<start>_<inst>.csv[.gz])
  part3d.csv     partition WSPDs: CK exact, the sweep, and the exact octree, converted (fewest/min/merge), and reshape from
                 both starts, converted (replays of the first run's reshape runs, or converted in the run itself)
  run3d.log      peak memory of the reshape runs of the first run (the later ones: peakBeforeKB of their converted row)
Reference of an instance: OPT where CPLEX proved it, else LB.
Outputs: three_d/results/tables/exp3d.tex (tab:exp3d), three_d/results/figs/{inputs,lb,summary,ratio,per_point,time,ls,lstotal,
partconv}3d.pdf, DATA3D/summary.md (every number of the text); the folder in RESULTS3D if set."""
import os
import re
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "analysis"))   # figutil.py
from figutil import ordered_handles, top_to_bottom, METHOD, INPUT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # three_d/
D3 = os.environ.get("DATA3D", os.path.join(ROOT, "data"))
PAPER = os.environ.get("RESULTS3D", os.path.join(ROOT, "results"))
for _sub in ("tables", "figs"):
    os.makedirs(os.path.join(PAPER, _sub), exist_ok=True)
CPP3D = os.path.join(ROOT, "src")
PARTIAL = os.environ.get("PARTIAL") == "1"
COLS = ["gen", "n", "seed", "eps", "algo", "size", "valid", "time", "info"]
SYN = ["uniform", "gauss", "cluster", "grid", "spiral", "shell"]
SCANS = ["bunny", "armadillo", "dragon", "buddha"]
LABEL = {"uniform": "uniform", "gauss": "Gaussian", "cluster": "clustered", "grid": "perturbed grid",
         "spiral": "conical spiral", "shell": "spherical shell", "bunny": "bunny", "armadillo": "armadillo",
         "dragon": "dragon", "buddha": "Buddha", "scans": "scans"}
ALGOS = {"ck": "CK box", "ckexact": "CK", "sweep": "Sweep", "lsck": "CK + reshape", "lssweep": "Sweep + reshape"}
CONS = ["ck", "ckexact", "sweep"]          # the constructions (bold: smallest of the row)
STYLE = {"sweep": dict(marker="o", ls="-", color=METHOD["blue"]),
         "ck": dict(marker="s", ls="--", color=METHOD["gold"]),
         "ckexact": dict(marker="D", ls=":", color=METHOD["green"]),
         "opt": dict(marker="^", ls="-.", color=METHOD["gray"]),
         "lssweep": dict(marker="x", ls="-", color=METHOD["magenta"]),
         "lsck": dict(marker="+", ls="--", color=METHOD["pink"])}
LB_STYLE = {"steps": dict(marker="o", ls=":", color=METHOD["blue"]), "lb1": dict(marker="s", ls="--", color=METHOD["gold"]),
            "lb": dict(marker="D", ls="-", color=METHOD["green"]), "alpha": dict(marker="^", ls="-.", color=METHOD["gray"])}
LB_LABEL = {"steps": r"$\mathrm{LB}_0$ (steps)", "lb1": r"$\mathrm{LB}_1$ (swaps)", "lb": r"$\mathrm{LB}$ (plateau moves)",
            "alpha": r"$\alpha$ (exact)"}
COLOR = {g: INPUT[g] for g in ("uniform", "gauss", "cluster", "grid", "spiral")}
COLOR.update(shell=INPUT["annulus"], scans=INPUT["CG:SHOP"])   # the colors of the planar analogues
INPUTS = SYN + ["scans"]
TABLE_NS = [1000, 10000]
BUDGET = {200: 100, 1000: 100, 10000: 1000}                     # local search: budget of the samples (s)
TMIN = {200: 1e-3, 1000: 1e-2, 10000: 1}                         # left end of the time axis of fig:ls3d
TOTAL_NS = (1000, 10000)
NEAR = 1.001
FIGW = 6.2
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "lines.markersize": 2.2, "lines.linewidth": 0.9, "axes.linewidth": 0.6,
                     "grid.linewidth": 0.4, "legend.frameon": False})
lines_out = []


def say(s=""):
    lines_out.append(s)


def read(name, cols=COLS):
    p = os.path.join(D3, name)
    return pd.read_csv(p, names=cols) if os.path.exists(p) and os.path.getsize(p) else pd.DataFrame(columns=cols)


def iname(gen, n, seed):
    return f"{gen}_{n}_{seed}"


def used(df):
    """seed 1 of the generators and the scans"""
    return df[(df.seed == 1) | ~df.gen.isin(SYN)].copy()


def num(info, key):
    m = re.search(key + r"=([0-9.eE+-]+)", str(info))
    return float(m.group(1)) if m else np.nan


def load():
    sep, first = used(read("exp3d_sep.csv")), used(read("exp3d.csv"))
    for x in (sep, first):
        x["inst"] = [iname(g, n, s) for g, n, s in zip(x.gen, x.n, x.seed)]
    # the runs in one process (first run) and in separate processes must agree
    m = sep.merge(first, on=["inst", "algo"], suffixes=("", "_1"))
    assert (m["size"] == m["size_1"]).all(), "sep and first run differ"
    ms = m[m.algo == "sweep"]
    assert (ms["info"].str.extract(r"steps=(\d+)")[0].values == ms["info_1"].str.extract(r"steps=(\d+)")[0].values).all()
    E = pd.concat([sep, first[~first.inst.isin(set(sep.inst))]], ignore_index=True)
    if not PARTIAL:
        assert set(first.inst) <= set(sep.inst), "missing sep runs"
    # lower bounds
    L = used(read("lb3d_sep.csv", ["inst", "gen", "n", "seed", "check", "n2", "anchors", "free", "swaps", "swapFree", "lb1",
                                   "passes", "converged", "capHits", "tSwap", "lb", "passes0", "swaps0", "swapFree0", "tSwap0"]))
    L0 = used(read("lb3d.csv", ["inst", "gen", "n", "seed", "steps", "lb1", "lb", "check", "tSwap"]))
    assert (L.check == "ok").all() and (L0.check == "ok").all()
    L = L.assign(steps=L.anchors)
    j = L.merge(L0, on="inst", suffixes=("", "_0"))
    assert (j.steps == j.steps_0).all() and (j.lb1 == j.lb1_0).all() and (j.lb == j.lb_0).all(), "lower bounds differ between runs"
    L = pd.concat([L, L0[~L0.inst.isin(set(L.inst))]], ignore_index=True).set_index("inst")
    o = used(read("opt3d.csv"))
    o["inst"] = [iname(g, n, s) for g, n, s in zip(o.gen, o.n, o.seed)]
    o = o[o["info"].str.contains("optimal=1")].drop_duplicates("inst").set_index("inst")
    ls = used(read("ls3d.csv"))
    ls["algo"] = ls.algo.map({"ls-ck-reshape": "lsck", "ls-sweep-reshape": "lssweep"})
    ls["inst"] = [iname(g, n, s) for g, n, s in zip(ls.gen, ls.n, ls.seed)]
    d = pd.concat([E[E.algo.isin(CONS)], ls], ignore_index=True)
    d = d.join(L[["steps", "lb1", "lb"]], on="inst").join(o["size"].rename("opt"), on="inst")
    d = d.dropna(subset=["lb"])
    assert (d.opt.isna() | (d.lb <= d.opt)).all() and (d["size"] >= d.lb).all()
    d["reflb"] = d.opt.isna()
    d["ref"] = d.opt.fillna(d.lb)
    d["ratio"] = d["size"] / d.ref
    d["pairs"] = d.n * (d.n - 1) / 2 / d.ref
    d["inp"] = d.gen.where(d.gen.isin(SYN), "scans")
    return d, L, o


def kn(n):
    return f"{n // 1000}k" if n >= 1000 and n % 1000 == 0 else str(n)


def kfmt(n):
    return rf"{n // 1000}\mathrm{{k}}" if n >= 1000 and n % 1000 == 0 else str(n)


def xticks(ax, xt=(200, 500, 1000, 2000, 5000, 10000)):
    ax.set_xscale("log"); ax.set_xticks(list(xt)); ax.set_xticklabels([kn(n) for n in xt]); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def logy(ax):
    if any(len(l.get_xdata()) for l in ax.get_lines()):   # no data yet (PARTIAL): keep the linear scale
        ax.set_yscale("log")


def refmarks(ax, g, colr, marker="o", **kw):
    """open markers where the reference is LB for every instance of that n, half-filled where for some"""
    ax.plot(g.index[g.lb], g.r[g.lb], ls="none", marker=marker, mfc="white", mec=colr, mew=0.7, **kw)
    m = g.mix & ~g.lb
    if m.any():
        ax.plot(g.index[m], g.r[m], ls="none", marker=marker, fillstyle="left", mfc=colr, mfcalt="white", mec=colr, mew=0.7, **kw)


def panel_title(ax, t, bottom=False):
    ax.text(0.03, 0.05 if bottom else 0.95, t, transform=ax.transAxes, ha="left", va="bottom" if bottom else "top",
            fontsize=plt.rcParams["axes.titlesize"], bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))


def save(fig, name):
    fig.savefig(os.path.join(PAPER, "figs", name)); plt.close(fig)


# ---------------------------------------------------------------- table
def table(d):
    g = d.groupby(["gen", "n", "algo"])["ratio"].mean()
    pairs = d[d.algo == "sweep"].groupby(["gen", "n"])["pairs"].mean()
    lbf = d.groupby(["gen", "n"])["reflb"].agg(list)
    out = [r"\begin{table}[htbp]", r"\centering\footnotesize\setlength{\tabcolsep}{4pt}",
           r"\caption{Size divided by the reference in $\R^3$, per input (one instance each). The column ``Pairs'' gives "
           r"$\binom n2$ divided by the reference. $\dagger$: the reference is the lower bound. Bold: smallest construction "
           r"of the row.}\label{tab:exp3d}",
           r"\begin{tabular}{@{}lrr" + "c" * len(ALGOS) + "@{}}", r"\toprule",
           r"Input & $n$ & Pairs & " + " & ".join(r"\shortstack{\strut " + v.replace(" + ", r"\\\strut + ") + "}"
                                                  for v in ALGOS.values()) + r" \\", r"\midrule"]
    for block in (SYN, SCANS):
        rows = [r for r in block if any((r, n) in pairs.index for n in TABLE_NS)]
        for i, r in enumerate(rows):
            for k, n in enumerate(n for n in TABLE_NS if (r, n) in pairs.index):
                vals = {a: f"{g.loc[(r, n, a)]:.2f}" for a in ALGOS if (r, n, a) in g.index}
                best = min((vals[a] for a in CONS if a in vals), key=float)
                cells = [(rf"\textbf{{{vals[a]}}}" if a in CONS and vals[a] == best else vals[a]) if a in vals else "--" for a in ALGOS]
                fl = np.array(lbf.loc[(r, n)])
                mark = r"$^\dagger$" if fl.all() else r"$^\ddagger$" if fl.any() else r"\phantom{$^\dagger$}"
                out.append((LABEL[r] if k == 0 else "") + f" & {kn(n)}{mark} & {pairs.loc[(r, n)]:.1f} & " + " & ".join(cells) + r" \\")
            if i < len(rows) - 1:
                out.append(r"\addlinespace[2pt]")
        if block is SYN:
            out.append(r"\midrule")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(PAPER, "tables", "exp3d.tex"), "w").write("\n".join(out) + "\n")


# ---------------------------------------------------------------- figures
def inputs_fig():
    """one instance with n = 1000 per generator (seed 1) and one sample of each scan; points dumped by three_d/src/exp3d"""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    items = []
    with tempfile.TemporaryDirectory() as tmp:
        for g in SYN + SCANS:
            gen = g if g in SYN else f"file:{os.path.join(D3, 'scans', g + '.xyz')}"
            pre = os.path.join(tmp, g)
            subprocess.run([os.path.join(CPP3D, "exp3d"), gen, "1000", "1", "1", "20", "ck"], cwd=CPP3D, check=True,
                           stdout=subprocess.DEVNULL, env=dict(os.environ, WSPD_DUMP=pre))
            items.append((g, np.loadtxt(pre + ".pts")))
    fig = plt.figure(figsize=(FIGW, 2.9))
    axes = []
    for k, (g, P) in enumerate(items):
        assert len(P) == 1000, (g, len(P))
        ax = fig.add_subplot(2, 5, k + 1, projection="3d")
        Q = P if g in SYN else P[:, [0, 2, 1]]   # scans: their y axis is the vertical one
        lo, hi = Q.min(axis=0), Q.max(axis=0); c, h = (lo + hi) / 2, 0.5 * (hi - lo).max()
        ax.scatter(Q[:, 0], Q[:, 1], Q[:, 2], s=0.4, c="black", linewidths=0, depthshade=False)
        ax.set_xlim(c[0] - h, c[0] + h); ax.set_ylim(c[1] - h, c[1] + h); ax.set_zlim(c[2] - h, c[2] + h)
        ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=20, azim=-60)
        ax.set_axis_off(); ax.set_title(LABEL[g], pad=-2)
        axes.append((g, ax))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01, wspace=0.08, hspace=0.18)
    for g, ax in axes:   # frame in the color of the input, as fig:inputs
        b = ax.get_position()
        fig.add_artist(plt.Rectangle((b.x0 + 0.005, b.y0 + 0.005), b.width - 0.01, b.height - 0.01, transform=fig.transFigure,
                                     fill=False, edgecolor=COLOR[g if g in SYN else "scans"], linewidth=1.0))
    save(fig, "inputs3d.pdf")


def summary_fig(d):
    """mean size / reference over all instances of each n, one line per method"""
    fig, ax = plt.subplots(figsize=(FIGW, 2.0))
    for a in ["ck", "ckexact", "sweep", "lsck", "lssweep"]:
        g = d[d.algo == a].groupby("n").agg(r=("ratio", "mean"), lb=("reflb", "all"), mix=("reflb", "any"))
        st = STYLE[a]
        if a == "lsck":   # coincides with sweep + reshape: dotted, on top of it, without markers (as fig:summary)
            ax.plot(g.index, g.r, label=ALGOS[a], color=st["color"], ls=":", lw=1.6, zorder=3)
            continue
        ax.plot(g.index, g.r, label=ALGOS[a], **st)
        refmarks(ax, g, st["color"], st["marker"])
    xticks(ax); ax.set_xlabel("$n$"); ax.set_ylabel("size / reference"); ax.grid(alpha=.3)
    ax.legend(*ordered_handles(ax), loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout(); save(fig, "summary3d.pdf")


def ratio_fig(d):
    """one panel per algorithm (CK, sweep, sweep + reshape), one line per input"""
    fa = {"ckexact": "CK", "sweep": "Sweep", "lssweep": "Sweep + reshape"}
    size = (FIGW, 2.2)
    fig, axes = plt.subplots(1, 3, figsize=size, sharey=True)
    for ax, a in zip(axes, fa):
        s = d[d.algo == a]
        for inp in INPUTS:
            g = s[s.inp == inp].groupby("n").agg(r=("ratio", "mean"), lb=("reflb", "all"), mix=("reflb", "any"))
            if not len(g): continue
            ax.plot(g.index, g.r, color=COLOR[inp], label=LABEL[inp], marker="o")
            refmarks(ax, g, COLOR[inp])
        xticks(ax, (200, 1000, 10000)); ax.grid(alpha=.3); ax.set_xlabel("$n$")
        panel_title(ax, fa[a])
    axes[0].set_ylabel("size / reference")
    h, l = ordered_handles(list(axes))
    fig.legend(h, l, loc="lower center", ncol=4)   # the 3D labels are longer than the planar ones: two rows
    fig.tight_layout(rect=(0, 0.36 / size[1], 1, 1)); save(fig, "ratio3d.pdf")


def per_point_fig(d):
    """(a) dumbbells per point of sweep + reshape, (b) pairs per dumbbell of the reference, per input"""
    size = (FIGW, 2.1)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=size)
    s = d[d.algo == "lssweep"]
    for inp in INPUTS:
        x = s[s.inp == inp]
        if not len(x): continue
        g = (x["size"] / x.n).groupby(x.n).mean()
        ax.plot(g.index, g.values, color=COLOR[inp], marker="o", label=LABEL[inp])
        g = x.pairs.groupby(x.n).mean()
        bx.plot(g.index, g.values, color=COLOR[inp], marker="o", label=LABEL[inp])
    ax.set_ylabel("dumbbells / point"); bx.set_ylabel(r"$\binom{n}{2}$ / reference")
    bx.set_yscale("log")
    ticks = [1, 2, 5, 10, 20, 50, 100, 200]
    bx.set_yticks(ticks); bx.set_yticklabels([str(t) for t in ticks]); bx.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for i, a in enumerate((ax, bx)):
        xticks(a); a.grid(alpha=.3); a.set_xlabel("$n$"); panel_title(a, "(%s)" % "ab"[i])
    h, l = ordered_handles([ax, bx])
    fig.legend(h, l, loc="lower center", ncol=4)
    fig.tight_layout(rect=(0, 0.36 / size[1], 1, 1)); save(fig, "per_point3d.pdf")


def trace(prefix, start, inst, row_size=None, ls_time=None):
    p = os.path.join(D3, "traces", f"{prefix}{start}_{inst}.csv")
    if not os.path.exists(p) and os.path.exists(p + ".gz"):   # the shipped traces are gzipped
        p += ".gz"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    t = pd.read_csv(p, header=None, names=["t", "size"])
    if row_size is not None and row_size < t["size"].iloc[-1]:   # the final redundancy pass is not in the trace
        tb = float(round(ls_time))
        t = pd.concat([t[t.t < tb], pd.DataFrame({"t": [tb], "size": [row_size]})], ignore_index=True)
    return t


def size_at(t, x):
    k = np.searchsorted(t.t.values, x, side="right") - 1
    return t["size"].values[max(k, 0)]


def time_to_near(d):
    """mean time for sweep + reshape to reach NEAR * OPT (instances with proven OPT): the time of the run before the local
    search (sweep, setup, final pass) plus the first time of the trace within NEAR * OPT; the sweep alone if it already is"""
    out = []
    for r in d[(d.algo == "lssweep") & d.opt.notna()].itertuples():
        t = trace("", "sweep", r.inst)
        if t is None:
            if PARTIAL: continue
            raise SystemExit(f"no trace for {r.inst}")
        hit = t[t["size"] <= NEAR * r.opt]
        if not len(hit):
            raise SystemExit(f"{r.inst}: sweep + reshape never within {NEAR} OPT")
        if hit.t.iloc[0] == 0:
            out.append((r.n, num(r.info, "tBuild"))); continue
        out.append((r.n, r.time - num(r.info, "lsTime") + hit.t.iloc[0]))
    return pd.DataFrame(out, columns=["n", "t"]).groupby("n")["t"].mean()


def peak_memory():
    """peak memory (GB) per (algo, n): the separate runs of mem3d.csv for the constructions; reshape from the log of the
    first run (no conversion in those runs) and from peakBeforeKB of the converted rows of the later runs"""
    rows = []
    M = read("mem3d.csv", ["inst", "algo", "peak", "wall"])
    for r in M.itertuples():
        g, n, s = r.inst.rsplit("_", 2)
        if int(s) == 1 or g not in SYN:
            rows.append((r.algo, int(n), r.peak))
    for line in open(os.path.join(D3, "run3d.log")):
        m = re.search(r" OK ls:(ck|sweep):(\S+)_(\d+)_(\d) \d+s peak=([\d.]+)GB", line)
        if m and line[:16] < "2026-10-28 22:11" and (m.group(4) == "1" or m.group(2) not in SYN):
            rows.append(("ls" + m.group(1), int(m.group(3)), float(m.group(5))))
    P = used(read("part3d.csv"))
    for r in P[P.algo.isin(["ls-ck-reshape+part", "ls-sweep-reshape+part"])].itertuples():
        if not np.isnan(num(r.info, "peakBeforeKB")):
            rows.append(("ls" + r.algo.split("-")[1], r.n, num(r.info, "peakBeforeKB") / 2 ** 20))
    return pd.DataFrame(rows, columns=["algo", "n", "rss"])


def time_fig(d, o):
    size = (FIGW, 2.0)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=size, gridspec_kw=dict(width_ratios=[1.25, 1]))
    for a in ["sweep", "lssweep", "ck", "ckexact", "opt"]:
        if a == "lssweep":
            s = time_to_near(d)
            ax.plot(s.index, s.values, label=f"Sweep + reshape\nto {NEAR:g}$\\,\\mathrm{{OPT}}$", **STYLE[a])
            continue
        if a == "opt":
            s = o.groupby("n")["time"].mean()
            ax.plot(s.index, s.values, label="Exact (CPLEX)", **STYLE[a])
            continue
        s = d[d.algo == a].groupby("n")["time"].mean()
        ax.plot(s.index, s.values, label=ALGOS[a], **STYLE[a])
    xticks(ax); logy(ax); ax.set_xlabel("$n$"); ax.set_ylabel("time (s)"); ax.grid(alpha=.3)
    ax.text(0.03, 0.95, "(a) time", transform=ax.transAxes, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"])
    mem = peak_memory()
    mem = mem[mem.n >= 5000]
    names = {"sweep": "Sweep", "ckexact": "CK", "lsck": "Reshape from CK", "lssweep": "Reshape from the sweep"}
    for a in ["sweep", "lssweep", "lsck", "ckexact"]:
        s = mem[mem.algo == a].groupby("n")["rss"].max()
        if len(s): bx.plot(s.index, s.values, label=names[a], **STYLE[a])
    xticks(bx, (5000, 7000, 10000))
    logy(bx); bx.set_xlabel("$n$"); bx.set_ylabel("peak memory (GB)"); bx.grid(alpha=.3)
    bx.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([1, 2, 3, 4, 6, 8, 10]))
    bx.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g")); bx.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    bx.text(0.03, 0.95, "(b) memory", transform=bx.transAxes, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"])
    h, l = ordered_handles([ax, bx], top_to_bottom(ax) + top_to_bottom(bx))
    fig.legend(h, l, loc="lower center", ncol=4)
    fig.tight_layout(rect=(0, 0.36 / size[1], 1, 1)); save(fig, "time3d.pdf")
    return mem


def lb_fig(d, L, alpha):
    size = (FIGW, 2.1)
    fig, axes = plt.subplots(1, 3, figsize=size)
    ax, bx, cx = axes
    sw = d[d.algo == "sweep"].set_index("inst")
    o = L.join(sw["opt"]).join(alpha.rename("alpha")).dropna(subset=["opt"])
    o = o[o.n <= 1000]
    for k in ("steps", "lb1", "lb", "alpha"):
        g = (o[k] / o.opt).groupby(o.n).mean().dropna()
        ax.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    ax.set_xticks(range(200, 1001, 200)); ax.set_xticklabels([kn(n) for n in range(200, 1001, 200)])
    ax.set_ylabel(r"bound / $\mathrm{OPT}$")
    u = L.join(d[d.algo == "lssweep"].set_index("inst")["size"].rename("ub")).dropna(subset=["ub"])
    for k in ("steps", "lb1", "lb"):
        g = (u.ub / u[k]).groupby(u.n).mean()
        bx.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    xticks(bx, (200, 1000, 10000)); bx.set_ylabel("sweep + reshape / bound")
    t = L.join(sw["info"].str.extract(r"tSweep=([0-9.]+)")[0].astype(float).rename("tSweep")).dropna(subset=["tSwap0", "tSweep"])
    for k, tt in (("steps", t.tSweep), ("lb1", t.tSweep + t.tSwap0), ("lb", t.tSweep + t.tSwap)):
        g = tt.groupby(t.n).mean()
        cx.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    xticks(cx, (200, 1000, 10000)); logy(cx); cx.set_ylabel("time (s)")
    for i, a in enumerate(axes):
        a.grid(alpha=.3); a.set_xlabel("$n$")
        panel_title(a, "(%s)" % "abc"[i], bottom=(i == 0))
    h, l = ordered_handles(list(axes), [LB_LABEL[k] for k in ("alpha", "lb", "lb1", "steps")])
    fig.legend(h, l, loc="lower center", ncol=len(h))
    fig.tight_layout(rect=(0, 0.22 / size[1], 1, 1)); save(fig, "lb3d.pdf")
    return o, u, t


def ls_runs(d):
    """reshape (ls3d.csv) and LNS (lns3d.csv) runs on the samples of fig:ls3d, with their traces"""
    R = d[d.algo.isin(["lsck", "lssweep"])].copy()
    R["start"] = R.algo.str[2:]; R["moves"] = "reshape"
    N = used(read("lns3d.csv"))
    N["inst"] = [iname(g, n, s) for g, n, s in zip(N.gen, N.n, N.seed)]
    N["start"] = N.algo.str.extract(r"^ls-(ck|sweep)-lns$")[0]; N["moves"] = "lns"
    N = N.join(d.drop_duplicates("inst").set_index("inst")[["ref", "reflb", "opt"]], on="inst")
    R = R[R.inst.isin(set(N.inst))]
    A = pd.concat([R, N], ignore_index=True)
    A["trace"] = [trace("" if r.moves == "reshape" else "lns_", r.start, r.inst, r.size, num(r.info, "lsTime")) for r in A.itertuples()]
    return A


def ls_fig(A):
    ns = [n for n in BUDGET if (A.n == n).any()]
    colors = {"reshape": METHOD["blue"], "lns": METHOD["gold"]}
    styles = {"ck": "--", "sweep": "-"}
    fig, axes = plt.subplots(1, max(len(ns), 1), figsize=(FIGW, 1.9), squeeze=False)
    last = {}
    for j, n in enumerate(ns):
        ax = axes[0][j]
        grid = np.logspace(np.log10(TMIN[n]), np.log10(BUDGET[n]), 80)
        for s in ("ck", "sweep"):
            for mv in ("reshape", "lns"):
                x = A[(A.start == s) & (A.moves == mv) & (A.n == n) & A.trace.notna()]
                if not len(x): continue
                y = np.mean([np.array([size_at(r.trace, t) for t in grid]) / r.ref for r in x.itertuples()], axis=0)
                ax.plot(grid, y, color=colors[mv], ls=styles[s])
                last.setdefault((s, mv), []).append(y[-1])
        ax.set_xscale("log"); ax.grid(alpha=.3)
        ax.set_title(rf"$n={kfmt(n)}$", loc="left", pad=2); ax.set_xlabel("time (s)")
        if j == 0: ax.set_ylabel("size / reference")
    end = lambda i, k: -np.mean([v for key, vs in last.items() if key[i] == k for v in vs])
    h = [Line2D([], [], color=colors[mv], label={"reshape": "Reshape", "lns": "LNS"}[mv]) for mv in sorted({k[1] for k in last}, key=lambda k: end(1, k))] + \
        [Line2D([], [], color="k", ls=styles[s], label={"ck": "from CK", "sweep": "from the sweep"}[s]) for s in sorted({k[0] for k in last}, key=lambda k: end(0, k))]
    if h: fig.legend(handles=h, loc="lower center", ncol=len(h), fontsize=6)
    fig.tight_layout(rect=(0, 0.1, 1, 1)); save(fig, "ls3d.pdf")


def lstotal_fig(d):
    """reshape from both starts against the total time (build of the start, setup, local search), n = 1k and 10k"""
    A = d[d.algo.isin(["lsck", "lssweep"]) & d.n.isin(TOTAL_NS)].copy()
    A["start"] = A.algo.str[2:]
    A["tBuild"] = [num(i, "tBuild") for i in A["info"]]
    A["tStart"] = A.tBuild + np.array([num(i, "lsSetup") for i in A["info"]])
    A["trace"] = [trace("", r.start, r.inst, r.size, num(r.info, "lsTime")) for r in A.itertuples()]
    if not PARTIAL:
        assert A.trace.notna().all()
    A = A[A.trace.notna()]
    styles = {"ck": "--", "sweep": "-"}
    size = (FIGW, 3.067)
    fig, axes = plt.subplots(1, len(TOTAL_NS), figsize=size, squeeze=False)
    ck0 = {}
    for j, n in enumerate(TOTAL_NS):
        ax = axes[0][j]
        for s in ("ck", "sweep"):
            for inp in INPUTS:
                x = A[(A.start == s) & (A.n == n) & (A.inp == inp)]
                if not len(x): continue
                t0 = x.tStart.mean()
                grid = np.logspace(np.log10(x.tBuild.mean()), np.log10(t0 + BUDGET[n]), 200)
                y = np.mean([np.array([size_at(r.trace, t) for t in grid - t0]) / r.ref for r in x.itertuples()], axis=0)
                c = COLOR[inp]
                ax.plot(grid, y, color=c, ls=styles[s], lw=0.7)
                ax.plot(grid[0], y[0], ls="none", marker="o", ms=2.2, mec=c, mew=0.6, mfc="white" if s == "ck" else c, zorder=3)
                if s == "ck": ck0.setdefault(inp, y[0])
        ax.set_xscale("log"); ax.grid(alpha=.3)
        ax.set_title(rf"$n={kfmt(n)}$", loc="left", pad=2); ax.set_xlabel("total time (s)")
        if j == 0: ax.set_ylabel("size / reference")
    h = [Line2D([], [], color="k", ls=styles[s], marker="o", ms=2.2, mew=0.6, mfc="white" if s == "ck" else "k",
                label={"ck": "from CK", "sweep": "from the sweep"}[s]) for s in ("ck", "sweep")]
    leg = axes[0][-1].legend(handles=h, loc="upper right", fontsize=6, numpoints=2, handlelength=3.6, handletextpad=0.5, borderaxespad=0.3)
    for ln in leg.get_lines(): ln.set_markevery([0])
    h = [Line2D([], [], color=COLOR[i], label=LABEL[i]) for i in sorted(ck0, key=lambda k: -ck0[k])]
    fig.legend(handles=h, loc="lower center", ncol=len(h), fontsize=6)
    fig.tight_layout(rect=(0, 0.207 / size[1], 1, 1)); save(fig, "lstotal3d.pdf")
    return A


def part_wide(d):
    """one row per instance: the cover and partition sizes of fig:partconv3d"""
    P = used(read("part3d.csv"))
    if not len(P):
        return None
    assert (P.valid == "ok").all()
    P["inst"] = [iname(g, n, s) for g, n, s in zip(P.gen, P.n, P.seed)]
    P = P.drop_duplicates(["inst", "algo"], keep="last")
    piv = P.pivot(index="inst", columns="algo", values="size")
    bm = P[P.algo.str.endswith("+part")].assign(bm=lambda x: [num(i, "beforeMerge") for i in x["info"]]).pivot(index="inst", columns="algo", values="bm")
    W = d.drop_duplicates("inst").set_index("inst")[["ref", "reflb", "n", "inp", "gen"]]
    W = W.join(d[d.algo.isin(["lsck", "lssweep"])].pivot(index="inst", columns="algo", values="size"))
    col = lambda df, c: df[c].reindex(W.index) if c in df else np.nan
    W["rs"] = W.get("lssweep"); W["rc"] = W.get("lsck")
    W["rs_p"] = col(piv, "ls-sweep-reshape+part"); W["rs_pn"] = col(bm, "ls-sweep-reshape+part")
    W["rc_p"] = col(piv, "ls-ck-reshape+part"); W["rc_pn"] = col(bm, "ls-ck-reshape+part")
    W["sw_p"] = col(piv, "sweep+part"); W["sw_pn"] = col(bm, "sweep+part")
    W["cke"] = col(piv, "ckexact"); W["cke_m"] = col(piv, "ckexact+part")
    W["qte"] = col(piv, "qtexact"); W["qte_m"] = col(piv, "qtexact+part")
    for a, c in (("ls-sweep-reshape", "rs"), ("ls-ck-reshape", "rc")):   # the replays reproduce the runs of ls3d.csv
        if a in piv:
            x = piv[a].dropna()
            assert (W.loc[x.index, c] == x).all(), f"replay of {a} differs"
    return W


def partconv_fig(W):
    size = (FIGW, 5.0)
    fig = plt.figure(figsize=size)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.45], hspace=0.5, wspace=0.06)
    ax = fig.add_subplot(gs[0, 0]); bx = fig.add_subplot(gs[0, 1], sharey=ax); cx = fig.add_subplot(gs[1, :])
    for p_, c, t in ((ax, "rs", "(a) Sweep + reshape (cover)"), (bx, "rs_p", "(b) Converted + merge (partition)")):
        for inp in INPUTS:
            x = W[W.inp == inp]
            g = pd.DataFrame({"r": x[c] / x.ref, "lb": x.reflb.astype(bool), "n": x.n}).dropna(subset=["r"]).groupby("n").agg(
                r=("r", "mean"), lb=("lb", "all"), mix=("lb", "any"))
            if not len(g): continue
            p_.plot(g.index, g.r, color=COLOR[inp], marker="o", label=LABEL[inp])
            refmarks(p_, g, COLOR[inp])
        xticks(p_); p_.grid(alpha=.3); p_.set_xlabel("$n$")
        panel_title(p_, t)
    ax.set_ylabel("size / reference"); plt.setp(bx.get_yticklabels(), visible=False)
    h, l = ordered_handles([ax, bx])
    fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.08)
    y = (ax.get_position().y0 - 0.06 + cx.get_position().y1) / 2
    fig.legend(h, l, loc="center", ncol=len(h), bbox_to_anchor=(0.5, y), fontsize=6, columnspacing=1.0, handletextpad=0.4)
    lines = [("rs", "Sweep + reshape (cover)", dict(color=METHOD["magenta"], ls=":", marker="x")),
             ("sw_pn", "Sweep, converted", dict(color=METHOD["blue"], ls="--", marker="^")),
             ("sw_p", "Sweep, converted + merge", dict(color=METHOD["blue"], ls="-", marker="^")),
             ("rc_pn", "CK + reshape, converted", dict(color=METHOD["pink"], ls="--", marker="v", ms=3.8)),
             ("rc_p", "CK + reshape, converted + merge", dict(color=METHOD["pink"], ls="-", marker="v", ms=3.8)),
             ("rs_pn", "Sweep + reshape, converted", dict(color=METHOD["magenta"], ls="--", marker="o")),
             ("rs_p", "Sweep + reshape, converted + merge", dict(color=METHOD["magenta"], ls="-", marker="o")),
             ("cke", "CK", dict(color=METHOD["green"], ls="--", marker="D")),
             ("cke_m", "CK + merge", dict(color=METHOD["green"], ls="-", marker="D")),
             ("qte_m", "Octree + merge", dict(color=METHOD["gold"], ls="-", marker="s"))]
    for c, lab, st in lines:
        g = pd.DataFrame({"r": W[c] / W.ref, "lb": W.reflb.astype(bool), "n": W.n}).dropna(subset=["r"]).groupby("n").agg(
            r=("r", "mean"), lb=("lb", "all"), mix=("lb", "any"))
        if not len(g): continue
        cx.plot(g.index, g.r, label=lab, **st)
        if st["marker"] != "x": refmarks(cx, g, st["color"], st["marker"], **({"ms": st["ms"]} if "ms" in st else {}))
    xticks(cx); cx.grid(alpha=.3); cx.set_xlabel("$n$"); cx.set_ylabel("size / reference")
    cx.legend(*ordered_handles(cx), loc="upper left", bbox_to_anchor=(0.0, 0.88))
    panel_title(cx, "(c) Mean over all inputs")
    save(fig, "partconv3d.pdf")


# ---------------------------------------------------------------- numbers of the text
def summary(d, L, o, A, W, mem, lbo, lbu, lbt):
    say("# Experiments in R^3: numbers for the text (analysis/exp3d.py; seed 1 of the generators and the scans)")
    say("")
    say("instances per n: " + ", ".join(f"{n}: {c}" for n, c in d[d.algo == "sweep"].groupby("n").size().items()))
    say("")
    say("## size / reference (mean, max) over all instances")
    for n, g in d.groupby("n"):
        say(f"n={n} (ref LB on {int(g[g.algo == 'sweep'].reflb.sum())}): " + "; ".join(
            f"{ALGOS[a]} {g[g.algo == a].ratio.mean():.4f} / {g[g.algo == a].ratio.max():.4f}" for a in ALGOS if (g.algo == a).any()))
    say("")
    say("## per input at n = 10000 (size / reference)")
    for inp, g in d[d.n == 10000].groupby("gen"):
        say(f"{inp}: " + "; ".join(f"{ALGOS[a]} {g[g.algo == a].ratio.mean():.3f}" for a in ALGOS if (g.algo == a).any())
            + f"; pairs/dumbbell {g[g.algo == 'sweep'].pairs.mean():.2f}; dumbbells/point (sweep+reshape) {(g[g.algo == 'lssweep']['size'] / 10000).mean():.1f}")
    say("")
    say("## smallest construction")
    c = d[d.algo.isin(CONS)].pivot_table(index="inst", columns="algo", values="size").join(d.drop_duplicates("inst").set_index("inst")[["n", "gen"]])
    c = c.dropna(subset=CONS)
    for n, g in c.groupby("n"):
        say(f"n={n}: sweep <= min(CK, CK box) on {(g.sweep <= g[['ck', 'ckexact']].min(axis=1)).sum()} of {len(g)}, < on {(g.sweep < g[['ck', 'ckexact']].min(axis=1)).sum()}; "
            f"CK < sweep on {(g.ckexact < g.sweep).sum()} ({', '.join(g[g.ckexact < g.sweep].gen)})")
    say("")
    say("## reshape: sweep start / CK start (final sizes); against OPT")
    r = d[d.algo.isin(["lsck", "lssweep"])].pivot_table(index="inst", columns="algo", values="size").dropna()
    if len(r):
        q = r.lssweep / r.lsck
        say(f"all: max |1 - sweep/ck| {np.abs(1 - q).max():.5f}")
    for a in ("lsck", "lssweep"):
        g = d[(d.algo == a) & d.opt.notna()]
        if len(g):
            say(f"{a} on proven OPT ({len(g)}): mean {(g['size'] / g.opt).mean():.5f} max {(g['size'] / g.opt).max():.5f}, = OPT on {(g['size'] == g.opt).sum()}")
    say("")
    say("## lower bounds")
    for n, g in lbo.groupby("n"):
        say(f"n={n} proven={len(g)}: LB0/OPT mean {(g.steps / g.opt).mean():.4f} min {(g.steps / g.opt).min():.4f}; LB1/OPT {(g.lb1 / g.opt).mean():.4f}; "
            f"LB/OPT mean {(g.lb / g.opt).mean():.5f} min {(g.lb / g.opt).min():.5f}, = OPT on {(g.lb == g.opt).sum()}; alpha/OPT mean {(g.alpha / g.opt).mean():.5f} "
            f"min {(g.alpha / g.opt).min():.5f} (N={g.alpha.notna().sum()}), LB = alpha on {(g.lb == g.alpha).sum()}")
    x = lbo.dropna(subset=["alpha"])
    if len(x):
        gap = x.alpha - x.steps
        say(f"alpha instances {len(x)}: of the total gap alpha - LB0, swaps close {(x.lb1 - x.steps).sum() / gap.sum():.3f} and plateau moves {(x.lb - x.lb1).sum() / gap.sum():.3f}")
    for n, g in lbu.groupby("n"):
        say(f"n={n}: LB/LB0 {(g.lb / g.steps).mean():.4f} (max {(g.lb / g.steps).max():.4f}); sweep+reshape / LB {(g.ub / g.lb).mean():.4f} (max {(g.ub / g.lb).max():.4f})")
    for n, g in lbt.groupby("n"):
        say(f"n={n}: time sweep {g.tSweep.mean():.2f} s, LB1 after the sweep {g.tSwap0.mean():.2f} s, LB after the sweep {g.tSwap.mean():.2f} s")
    say("")
    say("## OPT and alpha (CPLEX)")
    for n, g in o.groupby("n"):
        say(f"n={n}: proven {len(g)}, time median {g.time.median():.2f} s mean {g.time.mean():.1f} s max {g.time.max():.0f} s")
    a = used(read("alpha3d.csv"))
    if len(a):
        say(f"alpha: {len(a)} runs, optimal {a['info'].str.contains('optimal=1').sum()}, time max {a.time.max():.0f} s")
    say("")
    say("## running times (mean s) and peak memory (max GB)")
    for n, g in d.groupby("n"):
        say(f"n={n}: " + "; ".join(f"{ALGOS[a]} {g[g.algo == a].time.mean():.2f}" for a in CONS if (g.algo == a).any()))
    for (a, n), g in mem.groupby(["algo", "n"]):
        if n == 10000: say(f"memory n=10000 {a}: max {g.rss.max():.2f} GB")
    say("")
    if A is not None and len(A):
        say("## local search on the samples (fig:ls3d): mean size / reference at the end of the budget")
        for (n, s, mv), g in A.groupby(["n", "start", "moves"]):
            say(f"n={n} {mv} from {s} ({len(g)}): {(g['size'] / g.ref).mean():.5f}")
        for (n, s), g in A.groupby(["n", "start"]):
            p = g.pivot_table(index="inst", columns="moves", values="size").dropna()
            if len(p) and "lns" in p:
                say(f"n={n} from {s}: reshape <= LNS on {(p.reshape <= p.lns).sum()} of {len(p)}, equal on {(p.reshape == p.lns).sum()}")
        say("")
    if W is not None:
        say("## partition WSPDs (size / reference, mean)")
        cols = ["rs", "rs_pn", "rs_p", "rc_p", "sw_pn", "sw_p", "cke", "cke_m", "qte", "qte_m"]
        for n, g in W.groupby("n"):
            say(f"n={n}: " + "; ".join(f"{k} {(g[k] / g.ref).mean():.4f}" for k in cols if g[k].notna().any())
                + f"; N={g.rs_p.notna().sum()}")
        x = W.dropna(subset=["rs_p", "cke_m"])
        say(f"converted sweep + reshape (+ merge) larger than CK + merge on {(x.rs_p > x.cke_m).sum()} of {len(x)}; than CK on {(x.rs_p > x.cke).sum()}")
        x = W.dropna(subset=["qte_m", "cke_m"])
        say(f"octree + merge < CK + merge on {(x.qte_m < x.cke_m).sum()} of {len(x)}")
    open(os.path.join(D3, "summary.md"), "w").write("\n".join(lines_out) + "\n")
    print("\n".join(lines_out))


def main():
    d, L, o = load()
    alpha = used(read("alpha3d.csv"))
    alpha["inst"] = [iname(g, n, s) for g, n, s in zip(alpha.gen, alpha.n, alpha.seed)]
    alpha = alpha[alpha["info"].str.contains("optimal=1")].set_index("inst")["size"]
    assert (L.reindex(alpha.index).lb <= alpha).all(), "alpha below LB"
    table(d)
    summary_fig(d); ratio_fig(d); per_point_fig(d)
    mem = time_fig(d, o)
    lbo, lbu, lbt = lb_fig(d, L, alpha)
    A = ls_runs(d)
    ls_fig(A)
    lstotal_fig(d)
    W = part_wide(d)
    if W is not None: partconv_fig(W)
    if os.environ.get("NO_INPUTS") != "1": inputs_fig()
    summary(d, L, o, A, W, mem, lbo, lbu, lbt)


if __name__ == "__main__":
    main()

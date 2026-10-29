"""Local search over time: reshape (the heuristic of the paper) and LNS (the comparison), eps = 1/20.
Run from the repository root:  python3 analysis/ls_compare.py   (DATA=runs ... to analyze your own runs)
Inputs (DATA, default data/): ls_cmp_rs.csv (n = 1000, 10000; traces traces/R_*); for n = 50000, the LNS runs of
ls_cmp_big.csv and the reshape runs of the same instances in ls_all_big_5000.csv, both with a budget of 5000 s (traces
traces/RL_*); for the figure at n = 200 (not in the table), the runs of ls_all_rs.csv on all instances (traces
traces/RA_*); for the total-time figure, the reshape runs of ls_all_rs.csv on all instances with n = 1000 and 10000
(traces traces/RA_*) and the setup times of the local search (ls_setup.csv); the references (reference.py).
Outputs (RESULTS, default results/): tables/ls_time.tex (Table lstime), figs/ls_convergence.pdf (Fig. ls),
figs/ls_total_time_inputs.pdf (Fig. lstotal), ls_cmp_summary.md."""
import os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import DATA, RESULTS, KEY, MAX_N, refs, trace_files
from figutil import METHOD, INPUT

CSVS = [os.path.join(DATA, "ls_cmp_rs.csv"), os.path.join(DATA, "ls_cmp_big.csv")]
CSV_5000 = os.path.join(DATA, "ls_all_big_5000.csv")   # the reshape runs with n = BIG (traces RL_*)
CSV_ALL = os.path.join(DATA, "ls_all_rs.csv")   # convergence figure, n = FIG_SMALL, and total-time figure (all instances, traces RA_*)
CSV_SETUP = os.path.join(DATA, "ls_setup.csv")   # setup time of the local search per instance and start
FIG_SMALL = 200
BIG = MAX_N   # the largest n: LNS (ls_cmp_big.csv) and reshape (ls_all_big_5000.csv) with a budget of 5000 s, traces RL_*
BUDGET = {1000: 100, 10000: 1000, BIG: 5000}
FIG_BUDGET = {FIG_SMALL: 100, **BUDGET}   # time axis of the convergence figure
CHECK = {1000: [1, 10, 100], 10000: [10, 100, 1000], BIG: [50, 500, 5000]}   # checkpoints (s) of the time table
TMIN = {FIG_SMALL: 1e-3, 1000: 1e-2, 10000: 1, BIG: 1}   # left end of the time axis of the convergence figure (s)
EPSS = [0.05]
NAME = {"reshape": "Reshape", "lns": "LNS"}   # rows of the tables in this order
ORDER = list(NAME)
START = {"ck": "CK", "sweep": "Sweep"}
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "legend.fontsize": 6.5, "lines.linewidth": 0.9,
                     "axes.linewidth": 0.6, "legend.frameon": False})


def load(csvs=CSVS, prefix="R"):
    L = pd.concat([pd.read_csv(f) for f in csvs if os.path.exists(f)], ignore_index=True)
    bad = L[~L.valid.astype(str).str.startswith("ok")]
    if len(bad): print("WARNING: invalid outputs\n", bad[KEY + ["algo"]])
    L = L[L.valid.astype(str).str.startswith("ok")].copy()
    m = L.algo.str.extract(r"^ls-(ck|sweep)-(reshape|lns)$")
    L["start"], L["moves"] = m[0], m[1]
    L = L[L.moves.notna()]
    L = L.drop_duplicates(KEY + ["start", "moves"], keep="last")
    L["s0"] = L["info"].str.extract(r"start=(\d+)")[0].astype(float)
    L["tBuild"] = L["info"].str.extract(r"tBuild=([0-9.]+)")[0].astype(float)
    L["lsTime"] = L["info"].str.extract(r"lsTime=([0-9.]+)")[0].astype(float)   # the time the local search ran (its budget)
    opt, lb = refs()
    L = L.join(opt, on=KEY).join(lb, on=KEY)
    L["exact"] = L["opt"].notna()
    L["ref"] = L["opt"].fillna(L["lb"])
    L["r"] = L["size"] / L["ref"]
    L["r0"] = L["s0"] / L["ref"]
    L["tp"] = prefix   # trace file prefix
    return L


def load_samples():
    """The runs on the local-search samples: reshape and LNS for n = 1000 and 10000 (ls_cmp_rs.csv, traces R_*); for
    n = BIG, LNS (ls_cmp_big.csv) and the reshape runs of the same instances (ls_all_big_5000.csv), both with a budget
    of 5000 s (traces RL_*)."""
    L = load([CSVS[0]], "R")
    if os.path.exists(CSVS[1]):
        B = load([CSVS[1]], "RL")
        B = B[B.moves == "lns"]
        if os.path.exists(CSV_5000):
            R = load([CSV_5000], "RL")
            B = pd.concat([B, R[(R.moves == "reshape")].merge(B[KEY].drop_duplicates(), on=KEY)], ignore_index=True)
        L = pd.concat([L, B], ignore_index=True)
    return L


def fig_data(L):
    """Runs of the convergence figure: those of L, plus the runs on all instances with n = FIG_SMALL."""
    L = L.drop(columns="trace", errors="ignore")
    if os.path.exists(CSV_ALL):
        A = load([CSV_ALL], "RA")
        L = pd.concat([L, A[A.n == FIG_SMALL]], ignore_index=True)
    return attach_traces(L)


def attach_traces(L):
    """Attach the trace (time since the start of the local search, size) of every run, ending with the output size."""
    curves = []
    for row in L.itertuples():
        cg = row.gen.endswith(".instance")
        inst = "images_" + row.gen.replace(".instance", "") + "_instance_0_0_1_20" if cg else f"{row.gen}_{row.n}_{row.seed}_1_20"
        cand = trace_files(f"{row.tp}_{row.start}_{row.moves}_{inst}")
        if len(cand) != 1 or os.path.getsize(cand[0]) == 0: curves.append(None); continue
        t = pd.read_csv(cand[0], header=None, names=["t", "size"])
        # The trace stops before the final pass that deletes redundant dumbbells (it removes some after LNS, none after
        # reshape): the output is the size at the end of the budget, so the text, table, and figure agree.
        # The point is placed at the nominal budget (lsTime rounded), so that the last checkpoint sees it.
        if row.size < t["size"].iloc[-1]:
            tb = float(round(row.lsTime))
            t = pd.concat([t[t.t < tb], pd.DataFrame({"t": [tb], "size": [row.size]})], ignore_index=True)
        curves.append(t)
    L["trace"] = curves
    return L


def kn(n):
    """n in the tables: 1000 -> 1k, 50000 -> 50k (text mode)."""
    return f"{n // 1000}k" if n >= 1000 and n % 1000 == 0 else str(n)


def bfmt(x): return f"{x:.0f}" if x >= 10 else f"{x:.2g}"


def fmt(x, nd=3): return "--" if pd.isna(x) else f"{x:.{nd}f}"


def ftime(x):
    if pd.isna(x): return "--"
    if np.isinf(x): return r"$>$\,budget"
    return f"{x:.2g}" if x < 10 else f"{x:.0f}"


def size_at(t, x):
    """Size of the WSPD at time x of a trace (the last recorded size at or before x)."""
    k = np.searchsorted(t.t.values, x, side="right") - 1
    return t["size"].values[max(k, 0)]


def time_table(L):
    """Appendix table: mean size/reference after each checkpoint time of the local search. Columns: n x checkpoint; row
    blocks: start (build time, start, reshape, LNS). A run has no value at a checkpoint after the end of its budget:
    the cell is then "--"."""
    ns = [n for n in BUDGET if (L.n == n).any()]
    cols = [(n, e, x) for n in ns for e in EPSS for x in CHECK[n]]
    red = {}
    for row in L[L.trace.notna()].itertuples():
        for n, e, x in cols:
            if row.n == n and row.eps == e and x <= 1.001 * row.lsTime:
                red.setdefault((row.start, row.moves, n, e, x), []).append(size_at(row.trace, x) / row.ref)
    mean = {k: np.mean(v) for k, v in red.items() if not np.isnan(v).any()}
    k = len(CHECK[ns[0]])
    spec = "@{}l" + "".join(("@{\\hspace{8pt}}" if i % k == 0 else "") + "r" for i in range(len(cols))) + "@{}"
    out = [r"\begin{table}[htbp]", r"\centering\footnotesize\setlength{\tabcolsep}{3pt}",
           r"\caption{Local search over time: mean size divided by the reference. \emph{Build time}: median time to build"
           r" the start, not included in the times. Bold: best value for each $n$ and time over both starts.NOTE}",
           r"\label{tab:lstime}", r"\begin{tabular}{" + spec + "}", r"\toprule",
           "$n$ & " + " & ".join(rf"\multicolumn{{{k}}}{{c}}{{{kn(n)}}}" for n in ns) + r" \\",
           " ".join(rf"\cmidrule(lr){{{2 + k * i}-{1 + k * (i + 1)}}}" for i in range(len(ns))),
           r"Time & " + " & ".join(f"{x:g}\\,s" for _, _, x in cols) + r" \\"]
    bt = {s: L[L.start == s].drop_duplicates(KEY).groupby(["n", "eps"]).tBuild.median() for s in START}
    g0 = L.drop_duplicates(KEY + ["start"]).groupby(["start", "n", "eps"]).r0.mean()
    best = {c: min((v for (s_, m_, *c_), v in mean.items() if tuple(c_) == c), default=None) for c in cols}
    for s in START:
        out += [r"\midrule", rf"\multicolumn{{{1 + len(cols)}}}{{@{{}}l}}{{\emph{{Start: {START[s]}}}}} \\"]
        out.append(r"Build time & " + " & ".join(
            rf"\multicolumn{{{k}}}{{c}}{{{bfmt(bt[s][(n, e)])}\,s}}" if (n, e) in bt[s].index else rf"\multicolumn{{{k}}}{{c}}{{--}}"
            for n in ns for e in EPSS) + r" \\")
        out.append(r"Start & " + " & ".join(fmt(g0.get((s, n, e))) for n, e, _ in cols) + r" \\")
        for mv in ORDER:
            if not (L.moves == mv).any(): continue
            cells = []
            for c in cols:
                v = mean.get((s, mv) + c)
                if v is None: cells.append("--"); continue
                cells.append(rf"\textbf{{{fmt(v)}}}" if fmt(v) == fmt(best[c]) else fmt(v))   # ties on the displayed digits
            out.append(rf"{NAME[mv]} & " + " & ".join(cells) + r" \\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    out[2] = out[2].replace("NOTE", " --: not run." if any(" -- " in l + " " for l in out[9:]) else "")
    return out


def kfmt(n):
    """Panel label of n: 1000 -> 1k, 50000 -> 50k."""
    return rf"{n // 1000}\mathrm{{k}}" if n >= 1000 and n % 1000 == 0 else str(n)   # inside $...$ of the title


def convergence_fig(L, fname, size=(6.2, 1.9), title=r"$n={n}$", legend_inside=False, budget=BUDGET):
    """Mean size/reference over time (log scale), one panel per n with both starts: color = moves (reshape, LNS), line
    style = start (CK exact dashed, sweep solid). fname is relative to RESULTS/figs/ or absolute."""
    ns = [n for n in budget if (L.n == n).any() and L[L.n == n].trace.notna().any()]
    if not ns: return
    e = EPSS[0]
    colors = {"reshape": METHOD["blue"], "lns": METHOD["gold"]}
    styles = {"ck": "--", "sweep": "-"}
    fig, axes = plt.subplots(1, len(ns), figsize=size, squeeze=False)
    last = {}   # (start, moves) -> values at the end of the time axis, for the legend order
    for j, n in enumerate(ns):
        ax = axes[0][j]
        grid = np.logspace(np.log10(TMIN[n]), np.log10(budget[n]), 80)
        for s in START:
            for mv in ORDER:
                d = L[(L.start == s) & (L.moves == mv) & (L.n == n) & (L.eps == e) & L.trace.notna()]
                if not len(d): continue
                y = np.mean([np.array([size_at(r.trace, x) for x in grid]) / r.ref for r in d.itertuples()], axis=0)
                ax.plot(grid, y, color=colors[mv], ls=styles[s])
                last.setdefault((s, mv), []).append(y[-1])
        ax.set_xscale("log"); ax.grid(alpha=.3)
        lo, hi = ax.get_xlim()   # at least 3 labeled times: else every s-th decade from the first, largest s giving 3
        if sum(lo <= t <= hi for t in ax.get_xticks()) < 3:
            dec = np.arange(np.ceil(np.log10(lo)), np.floor(np.log10(hi)) + 1)
            s = max([s for s in range(1, len(dec) + 1) if len(dec[::s]) >= 3], default=1)
            ax.set_xticks(10.0 ** dec[::s])
        ax.set_title(title.format(n=kfmt(n)), loc="left", pad=2)
        ax.set_xlabel("time (s)")
        if j == 0: ax.set_ylabel("size / reference")
    from matplotlib.lines import Line2D
    # legend top to bottom: each entry by the mean of its curves at the end of the time axis
    end = lambda i, k: -np.mean([v for key, vs in last.items() if key[i] == k for v in vs])
    h = [Line2D([], [], color=colors[mv], label=NAME[mv]) for mv in sorted({k[1] for k in last}, key=lambda k: end(1, k))] + \
        [Line2D([], [], color="k", ls=styles[s], label={"ck": "from CK", "sweep": "from the sweep"}[s])
         for s in sorted({k[0] for k in last}, key=lambda k: end(0, k))]
    if legend_inside:
        axes[0][-1].legend(handles=h, loc="upper right", handlelength=1.6, handletextpad=0.4, borderaxespad=0.3)
        fig.tight_layout(pad=0.3)
    else:
        fig.legend(handles=h, loc="lower center", ncol=len(h), fontsize=6)
        fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(fname if os.path.isabs(fname) else os.path.join(RESULTS, "figs", fname)); plt.close(fig)


TOTAL_NS = (1000, 10000)   # the n of the total-time figure
INPUTS = {g: (l, INPUT[g]) for g, l in (("uniform", "uniform"), ("gauss", "Gaussian"), ("cluster", "clustered"),
          ("grid", "perturbed grid"), ("spiral", "spiral"), ("annulus", "annulus"), ("CG:SHOP", "CG:SHOP"))}   # label and color


def total_data():
    """Runs of the total-time figure: reshape from both starts on every instance with n in TOTAL_NS (ls_all_rs.csv,
    traces RA_*), with the columns inp (the generator, or CG:SHOP for the image-based instances) and tStart: the time
    until the first round of the local search, that is, the time to build the start (tBuild) plus the setup of the local
    search (lsSetup in ls_setup.csv, measured in separate runs). None if there are none."""
    if not os.path.exists(CSV_ALL) or not os.path.exists(CSV_SETUP): return None
    A = load([CSV_ALL], "RA")
    A = A[(A.moves == "reshape") & A.n.isin(TOTAL_NS)].copy()
    S = pd.read_csv(CSV_SETUP)
    S = S.assign(start=S.algo.str.extract(r"^ls-(ck|sweep)-")[0], setup=S["info"].str.extract(r"lsSetup=([0-9.]+)")[0].astype(float))
    A = attach_traces(A.join(S.drop_duplicates(KEY + ["start"], keep="last").set_index(KEY + ["start"]).setup, on=KEY + ["start"]))
    A = A[A.trace.notna() & A.setup.notna()]
    A = A[A.groupby(KEY).start.transform("nunique") == len(START)]   # instances with both starts
    return A.assign(tStart=A.tBuild + A.setup, inp=A.gen.where(~A.gen.str.endswith(".instance"), "CG:SHOP")) if len(A) else None


def total_time_fig(L, fname, size=(6.2, 3.067), title=r"$n={n}$"):
    """Reshape from CK exact and from the sweep against the total time: the time to build the start, plus the setup of
    the local search, plus the time of the local search. One panel per n, one curve per input (the mean over its seeds;
    CG:SHOP: over the image-based instances) and start: color = input, line style = start as in convergence_fig; the dot at the
    beginning of a curve is hollow for CK exact and filled for the sweep. A curve begins at the mean time to build the start (tBuild) and the mean size of the start, is horizontal during
    the setup of the local search (until the mean of tStart), and is then the mean size/reference at equal local-search
    time."""
    from matplotlib.lines import Line2D
    ns = [n for n in TOTAL_NS if (L.n == n).any()]
    e, styles = EPSS[0], {"ck": "--", "sweep": "-"}
    fig, axes = plt.subplots(1, len(ns), figsize=size, squeeze=False)
    ck0 = {}   # for the legend order: input -> size of CK exact in the first panel
    for j, n in enumerate(ns):
        ax = axes[0][j]
        for s in START:
            for inp, (_, c) in INPUTS.items():
                d = L[(L.start == s) & (L.n == n) & (L.eps == e) & (L.inp == inp)]
                if not len(d): continue
                t0 = d.tStart.mean()   # before t0 (during the setup), size_at gives the size of the start
                grid = np.logspace(np.log10(d.tBuild.mean()), np.log10(t0 + BUDGET[n]), 200)
                y = np.mean([np.array([size_at(r.trace, x) for x in grid - t0]) / r.ref for r in d.itertuples()], axis=0)
                ax.plot(grid, y, color=c, ls=styles[s], lw=0.7)
                ax.plot(grid[0], y[0], ls="none", marker="o", ms=2.2, mec=c, mew=0.6, mfc="white" if s == "ck" else c, zorder=3)   # hollow: CK exact
                if s == "ck": ck0.setdefault(inp, y[0])
        ax.set_xscale("log"); ax.grid(alpha=.3)
        ax.set_title(title.format(n=kfmt(n)), loc="left", pad=2)
        ax.set_xlabel("total time (s)")
        if j == 0: ax.set_ylabel("size / reference")
    # legends: the inputs below the panels, by decreasing size of CK exact (the hollow dots) in the first panel; the
    # starts inside the last panel, CK exact first, with long handles (so that the dashes show) that begin with the dot, as the curves do
    h = [Line2D([], [], color="k", ls=styles[s], marker="o", ms=2.2, mew=0.6, mfc="white" if s == "ck" else "k",
                label={"ck": "from CK", "sweep": "from the sweep"}[s]) for s in START]
    leg = axes[0][-1].legend(handles=h, loc="upper right", fontsize=6, numpoints=2, handlelength=3.6, handletextpad=0.5, borderaxespad=0.3)
    for ln in leg.get_lines(): ln.set_markevery([0])
    h = [Line2D([], [], color=INPUTS[i][1], label=INPUTS[i][0]) for i in sorted(ck0, key=lambda k: -ck0[k])]
    fig.legend(handles=h, loc="lower center", ncol=len(h), fontsize=6)
    fig.tight_layout(rect=(0, 0.207 / size[1], 1, 1))   # 0.207 in for the input legend; height: the axes are as tall as fig:partconv (c)
    fig.savefig(fname if os.path.isabs(fname) else os.path.join(RESULTS, "figs", fname)); plt.close(fig)


def summary(L):
    out = ["# Local-search comparison (reshape, LNS), eps = 1/20 (ls_cmp_rs.csv, ls_cmp_big.csv)", ""]
    g = L.groupby(["n", "eps", "start", "moves"]).agg(r=("r", "mean"), rmax=("r", "max"), N=("r", "count"),
                                                        exact=("exact", "all"))
    out.append(g.reset_index().to_markdown(index=False, floatfmt=".3f"))
    return out


def main():
    L = attach_traces(load_samples())
    open(os.path.join(RESULTS, "tables", "ls_time.tex"), "w").write("\n".join(time_table(L)) + "\n")
    convergence_fig(fig_data(L), "ls_convergence.pdf", budget=FIG_BUDGET)
    T = total_data()
    if T is not None: total_time_fig(T, "ls_total_time_inputs.pdf")   # needs runs with n = 1000 or 10000 and ls_setup.csv
    s = summary(L)
    open(os.path.join(RESULTS, "ls_cmp_summary.md"), "w").write("\n".join(s) + "\n")
    print("\n".join(s))


if __name__ == "__main__":
    main()

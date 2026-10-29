"""Tables and figures of the experiments section (main text and appendix) and the numbers quoted in the text.
Run from the repository root:  python3 analysis/make_tables_figs.py   (DATA=runs ... to analyze your own runs)

Inputs (DATA, default data/)
  exp2.csv      n = 200..1000 (step 100): the constructions (sweep = lazyprune, CK box = ck, CK exact = ckexact);
                exp2_opt.csv: exact OPT (CPLEX); optx.csv: exact OPT (optx)
  exp3.csv      n = 2000..10000 (step 1000);  exp4.csv: n = 20000..50000 (step 10000)
  ls_all_rs.csv local search from CK exact and from the sweep, n <= 10000; ls_all_big_5000.csv: n >= 20000 (5000 s)
  traces/RA_sweep_reshape_*  traces of sweep + reshape (time to reach 1.001 OPT)
  memory.csv    peak resident memory of every job
  lb_check.csv  the lower bounds (experiments/certify.py); alpha.csv: the largest sets of pairwise incompatible
                pairs for n = 200..500 (experiments/alpha.py)
Reference of an instance: OPT where proven, else the lower bound LB (lb_check.csv), see reference.py.
Outputs (RESULTS, default results/)
  figs/summary.pdf (Fig. summary), figs/ratio.pdf (Fig. ratio), figs/time.pdf (Fig. time), figs/lb.pdf (Fig. lb),
  tables/results.tex (Table exp, appendix), tables/summary.txt (numbers for the text).
The paper shows eps = 1/20 only (PAPER_EPS)."""
import os, re
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import DATA, RESULTS, KEY, MIN_N, MAX_N, read, trace_files, with_optx, proven_opt, sweep_lb, lower_bounds
from figutil import ordered_handles, top_to_bottom, METHOD, INPUT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAIN = "lazyprune"                                   # the practical polynomial-time sweep
ALGOS = {MAIN: "Sweep", "ck": "CK box", "ckexact": "CK"}      # the constructions compared in the tables
EXTRA = {"lsck": "CK + reshape", "lssweep": "Sweep + reshape", "lscklns": "CK + LNS",
         "lssweeplns": "Sweep + LNS"}
TABLE_ALGOS = {"ck": "CK box", "ckexact": "CK", MAIN: "Sweep", "lsck": "CK + reshape",
               "lssweep": "Sweep + reshape"}   # columns of the table
STYLE = {MAIN: dict(marker="o", ls="-", color=METHOD["blue"]),
         "ck": dict(marker="s", ls="--", color=METHOD["gold"]),
         "ckexact": dict(marker="D", ls=":", color=METHOD["green"]),
         "opt": dict(marker="^", ls="-.", color=METHOD["gray"]),
         "lssweep": dict(marker="x", ls="-", color=METHOD["magenta"]),
         "lsck": dict(marker="+", ls="--", color=METHOD["pink"])}
NEAR = 1.001                                         # time_fig: time for sweep + reshape to reach NEAR * OPT
SYN_GENS = ["uniform", "gauss", "cluster", "grid", "spiral", "annulus"]
PAPER_EPS = 0.05                                     # the eps of the tables and figures (inside the range of thm:sweep)
LABEL = {"uniform": "uniform", "gauss": "Gaussian", "cluster": "clustered", "grid": "perturbed grid", "spiral": "spiral",
         "annulus": "annulus"}                       # names used in the text
TABLE_NS = [1000, 10000, MAX_N]
FIGW = 6.2
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "lines.markersize": 2.2, "lines.linewidth": 0.9, "axes.linewidth": 0.6,
                     "grid.linewidth": 0.4, "legend.frameon": False})


def nominal_n(gen, n):
    m = re.search(r"-(\d{7})", gen)
    return int(m.group(1)) if gen.endswith(".instance") and m else int(n)


def OPTX(d):
    """Rows of optima proven by src/optx: their times are not CPLEX times."""
    return d["info"].astype(str).str.contains("solver=optx")


def load(*names):
    df = read(*names)
    if not len(df): return df
    df = with_optx(df)   # optima proven by src/optx where CPLEX has none (spiral n = 800..1000)
    bad = df[~df.valid.astype(str).str.startswith("ok")]
    if len(bad): print("WARNING: invalid outputs\n", bad)
    df = df[df.valid.astype(str).str.startswith("ok")].drop_duplicates(KEY + ["algo"], keep="last").copy()
    df["cgshop"] = df["gen"].str.endswith(".instance")
    df["family"] = [re.match(r"(.*?)-\d{7}", g).group(1) if c else g for g, c in zip(df.gen, df.cgshop)]
    return df


def with_ratios(df):
    """Per-instance reference: exact OPT when proven; otherwise the lower bound LB (reference.sweep_lb: the anchors of the sweep
    enlarged by src/lb_swap, DATA/lb_check.csv). Column reflb: True where the reference is the lower bound; exact: True
    where it is a proven OPT; steps: the number of steps of the sweep."""
    o = df[df.algo == "opt"].copy()
    if len(o):
        o["optimal"] = o["info"].str.contains("optimal=1")
        o = o.set_index(KEY)
    st = df[df.algo == MAIN].set_index(KEY)
    steps = st["info"].str.extract(r"steps=(\d+)")[0].astype(float).rename("steps")
    d = df[df.algo.isin(list(ALGOS) + list(EXTRA) + ["opt"])].join(sweep_lb(st), on=KEY).join(steps, on=KEY)
    if len(o):
        d = d.join(o["size"].rename("optsize"), on=KEY).join(o["optimal"], on=KEY)
        d["opt"] = d["optsize"].where(d["optimal"] == True)
        d = d.join(proven_opt(df[df.algo == "opt"]).rename("best"), on=KEY)
    else:
        d["opt"] = float("nan"); d["best"] = float("nan"); d["optimal"] = float("nan")
    d["reflb"] = d["best"].isna()
    d["ref"] = d["best"].fillna(d["lb"])
    d["exact"] = d["opt"].notna()
    d["ratio"] = d["size"] / d["ref"]
    return d


def kn(n):
    """n in the tables: 1000 -> 1k, 50000 -> 50k (text mode)."""
    return f"{n // 1000}k" if n >= 1000 and n % 1000 == 0 else str(n)


def refmark(flags):
    """Marker of a group from the reflb flags of its instances: '' (all OPT), dagger (all LB), double dagger (mixed)."""
    return r"$^\dagger$" if flags.all() else r"$^\ddagger$" if flags.any() else r"\phantom{$^\dagger$}"


def table_rows(d, rowkey):
    """Rows (input, n) for eps = PAPER_EPS, columns TABLE_ALGOS; mean of size/ref. Bold: smallest construction (Sweep,
    CK) of the row."""
    d = d[d.eps == PAPER_EPS].dropna(subset=["ratio"])
    g = d.groupby([rowkey, "nn", "algo"])["ratio"].mean()
    lbref = d.groupby([rowkey, "nn"])["reflb"].agg(list)
    out = []
    present = set(d[d.nn.isin(TABLE_NS)][rowkey])     # the inputs with a row
    # CG:SHOP: the families present at every n (euro-night, us-night) first, then the others by the smallest n shown, then by name
    first = ["euro-night", "us-night"]
    rest = sorted((r for r in present if r not in first), key=lambda r: (d[(d[rowkey] == r) & d.nn.isin(TABLE_NS)].nn.min(), r))
    rows = [r for r in (SYN_GENS if rowkey == "gen" else first + rest) if r in present]
    for i, r in enumerate(rows):
        ns = [n for n in TABLE_NS if (d[(d[rowkey] == r) & (d.nn == n)]).shape[0]]
        for k, n in enumerate(ns):
            vals = {a: f"{g.loc[(r, n, a)]:.2f}" for a in TABLE_ALGOS if (r, n, a) in g.index}
            cons = [vals[a] for a in ALGOS if a in vals]
            best = min(cons, key=float) if cons else None   # smallest displayed construction of the row (ties: all bold)
            cells = [(rf"\textbf{{{vals[a]}}}" if a in ALGOS and vals[a] == best else vals[a]) if a in vals else "--"
                     for a in TABLE_ALGOS]
            mark = refmark(np.array(lbref.loc[(r, n)]))
            out.append((LABEL.get(r, r) if k == 0 else "") + f" & {kn(n)}{mark} & " + " & ".join(cells) + r" \\")
        if i < len(rows) - 1 and ns: out.append(r"\addlinespace[2pt]")
    return out


def results_table(syn, cgs, fname, caption, label):
    """Appendix table: synthetic inputs, then the CG:SHOP families (one tabular, full width)."""
    head = [r"\begin{tabular}{@{}lr" + "c" * len(TABLE_ALGOS) + "@{}}", r"\toprule",
            r"Input & $n$ & " + " & ".join(r"\shortstack{\strut " + v.replace(" (", r"\\\strut (").replace(" + ", r"\\\strut + ")
                                             + "}" for v in TABLE_ALGOS.values()) + r" \\", r"\midrule"]
    out = [r"\begin{table}[htbp]", r"\centering\footnotesize\setlength{\tabcolsep}{4pt}", rf"\caption{{{caption}}}\label{{{label}}}"]
    out += head + table_rows(syn, "gen") + [r"\midrule"] + table_rows(cgs, "family") + [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(RESULTS, "tables", fname), "w").write("\n".join(out) + "\n")


def summary_fig(d, fname, size=(FIGW, 2.0)):
    """Main-text figure (fig:summary): mean size/reference over all instances of each n (eps = PAPER_EPS, log
    axes), one line per method of TABLE_ALGOS. Filled markers: reference OPT for all instances; open: the lower bound for
    all; half-filled: mixed. Returns the number of instances per n of TABLE_NS (for summary.txt)."""
    d = d[d.eps == PAPER_EPS].dropna(subset=["ratio"])
    fig, ax = plt.subplots(figsize=size)
    for a in ["ck", "ckexact", MAIN, "lsck", "lssweep"]:
        g = d[d.algo == a].groupby("nn").agg(r=("ratio", "mean"), lb=("reflb", "all"), mix=("reflb", "any"))
        st = STYLE[a]; colr, mk = st["color"], st["marker"]   # as in time_fig
        if a == "lsck":   # almost coincides with sweep + reshape: a dotted line on top of it, no markers
            ax.plot(g.index, g.r, label=TABLE_ALGOS[a], color=colr, ls=":", lw=1.6, zorder=3)
            continue
        ax.plot(g.index, g.r, label=TABLE_ALGOS[a], **st)
        ax.plot(g.index[g.lb], g.r[g.lb], ls="none", marker=mk, mfc="white", mec=colr, mew=0.7)
        m = g.mix & ~g.lb
        if m.any(): ax.plot(g.index[m], g.r[m], ls="none", marker=mk, fillstyle="left", mfc=colr, mfcalt="white", mec=colr, mew=0.7)
    ax.set_xscale("log"); ax.set_yscale("log")
    ticks = [1, 1.25, 1.5, 1.75, 2, 2.5]
    ax.set_yticks(ticks); ax.set_yticklabels([f"{t:g}" for t in ticks]); ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("$n$"); ax.set_ylabel("size / reference"); ax.grid(alpha=.3)
    ax.legend(*ordered_handles(ax), loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)
    return {n: d[(d.nn == n) & (d.algo == MAIN)].shape[0] for n in TABLE_NS if (d.nn == n).any()}


INPUT_STYLE = INPUT
FIG_ALGOS = {"ckexact": "CK", MAIN: "Sweep", "lssweep": "Sweep + reshape"}


def ratio_fig(d, fname, size=(FIGW, 2.2)):
    """eps = PAPER_EPS, n = MIN_N..MAX_N (log). Columns = algorithm; one line per input (mean over its seeds / instances).
    The CG:SHOP line is the mean over the image-based families. Filled markers: reference OPT for all instances; open: the
    lower bound for all; half-filled: mixed."""
    fig, axes = plt.subplots(1, len(FIG_ALGOS), figsize=size, sharey="row", squeeze=False)
    for col, a in enumerate(FIG_ALGOS):
        ax = axes[0][col]
        s = d[(d.eps == PAPER_EPS) & (d.algo == a)].dropna(subset=["ratio"])
        s = s.assign(inp=s.gen.where(~s.cgshop, "CG:SHOP"))
        for inp, colr in INPUT_STYLE.items():
            g = s[s.inp == inp].groupby("nn").agg(r=("ratio", "mean"), lb=("reflb", "all"), mix=("reflb", "any"))
            if not len(g): continue
            ax.plot(g.index, g.r, color=colr, label=LABEL.get(inp, inp), marker="o")
            ax.plot(g.index[g.lb], g.r[g.lb], ls="none", marker="o", mfc="white", mec=colr, mew=0.7)
            m = g.mix & ~g.lb
            if m.any(): ax.plot(g.index[m], g.r[m], ls="none", marker="o", fillstyle="left", mfc=colr, mfcalt="white", mec=colr, mew=0.7)
        ax.set_xscale("log")
        ax.text(0.03, 0.95, FIG_ALGOS[a], transform=ax.transAxes, ha="left", va="top",   # title inside the panel
                fontsize=plt.rcParams["axes.titlesize"], bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
        ax.grid(alpha=.3)
        ax.set_xlabel("$n$")
        if col == 0: ax.set_ylabel("size / reference")
    h, l = ordered_handles(list(axes[0]))
    fig.legend(h, l, loc="lower center", ncol=len(h))
    fig.tight_layout(rect=(0, 0.25 / size[1], 1, 1)); fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)


def per_point_fig(d, fname, algo="lssweep", size=(FIGW, 2.0)):
    """eps = PAPER_EPS, n = MIN_N..MAX_N (log): mean number of dumbbells per point (size / n, with n the number of distinct
    points) of sweep + reshape, one line per input (mean over its seeds; CG:SHOP: over the image-based families)."""
    fig, ax = plt.subplots(figsize=size)
    s = d[(d.eps == PAPER_EPS) & (d.algo == algo)]
    s = s.assign(inp=s.gen.where(~s.cgshop, "CG:SHOP"), pp=s["size"] / s.n)
    out = []
    for inp, colr in INPUT_STYLE.items():
        g = s[s.inp == inp].groupby("nn").pp.mean()
        if not len(g): continue
        ax.plot(g.index, g.values, color=colr, label=LABEL.get(inp, inp), marker="o")
        out.append(f"dumbbells per point, {EXTRA.get(algo, algo)}, {inp}: " + ", ".join(f"{n}: {v:.2f}" for n, v in g.items()))
    ax.set_xscale("log"); ax.grid(alpha=.3); ax.set_xlabel("$n$"); ax.set_ylabel("dumbbells / point")
    h, l = ordered_handles([ax])
    fig.legend(h, l, loc="lower center", ncol=len(h))
    fig.tight_layout(rect=(0, 0.25 / size[1], 1, 1)); fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)
    return out


LB_STYLE = {"steps": dict(marker="o", ls=":", color=METHOD["blue"]), "lb1": dict(marker="s", ls="--", color=METHOD["gold"]),
            "lb": dict(marker="D", ls="-", color=METHOD["green"]), "alpha": dict(marker="^", ls="-.", color=METHOD["gray"])}
LB_LABEL = {"steps": r"$\mathrm{LB}_0$ (steps)", "lb1": r"$\mathrm{LB}_1$ (swaps)", "lb": r"$\mathrm{LB}$ (plateau moves)",
            "alpha": r"$\alpha$ (exact)"}


def alpha_exact():
    """Series indexed by KEY: alpha, the largest set of pairwise incompatible pairs (experiments/alpha.py, DATA/alpha.csv;
    n = 200..500), where CPLEX proved it largest and check_lb accepted it."""
    a = read("alpha.csv")
    a = a[(a.optimal == 1) & (a.check == "ok")].copy()
    a["eps"] = PAPER_EPS
    return a.drop_duplicates(KEY, keep="last").set_index(KEY)["alpha"]


def lb_fig(d, fname, size=(FIGW, 2.1)):
    """fig:lb, eps = PAPER_EPS. The three lower bounds (DATA/lb_check.csv): LB_0 = number of steps of the sweep, LB_1 =
    after the free pairs and the 1-for-2 swaps, LB = after the plateau passes; alpha = largest set of pairwise
    incompatible pairs (n <= 500). (a) mean bound / OPT over the instances with proven OPT, n <= 1000; (b) mean of
    (sweep + reshape) / bound over all instances; (c) mean time: the sweep, the sweep and the swaps, the sweep and
    lb_swap with the plateau passes (the times of the runs of experiments/certify.py)."""
    L = lower_bounds()
    L = L[(L.index.get_level_values("eps") == PAPER_EPS) & (L.index.get_level_values("n") >= MIN_N)]
    sw = d[(d.algo == MAIN) & (d.eps == PAPER_EPS)].copy()
    sw["n"] = sw.nn
    sw = sw.set_index(KEY)
    al = alpha_exact()
    fig, axes = plt.subplots(1, 3, figsize=size)
    ax, bx, cx = axes
    # (a)
    o = L.join(sw["opt"]).join(al.rename("alpha")).dropna(subset=["opt"])
    o = o[o.index.get_level_values("n") <= 1000]
    nn = o.index.get_level_values("n")
    for k in ("steps", "lb1", "lb", "alpha"):
        g = (o[k] / o.opt).groupby(nn).mean().dropna()
        ax.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    ax.set_xticks(range(MIN_N, 1001, 200)); ax.set_xticklabels([kn(n) for n in range(MIN_N, 1001, 200)])
    ax.set_ylabel(r"bound / $\mathrm{OPT}$")
    # (b)
    u = L.join(d[(d.algo == "lssweep") & (d.eps == PAPER_EPS)].assign(n=lambda x: x.nn).set_index(KEY)["size"].rename("ub"))
    assert u.ub.notna().all()
    nn = u.index.get_level_values("n")
    for k in ("steps", "lb1", "lb"):
        g = (u.ub / u[k]).groupby(nn).mean()
        bx.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    bx.set_xscale("log"); bx.set_ylabel("sweep + reshape / bound")
    # (c)
    nn = L.index.get_level_values("n")
    for k, t in (("steps", L.tSweep), ("lb1", L.tSweep + L.tSwap0), ("lb", L.tSweep + L.tSwap)):
        g = t.groupby(nn).mean()
        cx.plot(g.index, g.values, label=LB_LABEL[k], **LB_STYLE[k])
    cx.set_xscale("log"); cx.set_yscale("log"); cx.set_ylabel("time (s)")
    for i, a in enumerate(axes):
        a.grid(alpha=.3); a.set_xlabel("$n$")
        a.text(0.03, 0.05 if i == 0 else 0.95, "(%s)" % "abc"[i], transform=a.transAxes, ha="left", va="bottom" if i == 0 else "top",   # (a): its lines are at the top
               fontsize=plt.rcParams["axes.titlesize"],
               bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
    h, l = ordered_handles(list(axes), [LB_LABEL[k] for k in ("alpha", "lb", "lb1", "steps")])
    fig.legend(h, l, loc="lower center", ncol=len(h))
    fig.tight_layout(rect=(0, 0.22 / size[1], 1, 1)); fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)
    lines = ["lower bounds (fig:lb):"]
    for n, g in o.groupby(level="n"):
        lines.append(f"  n={n} proven={len(g)} steps/OPT mean {(g.steps / g.opt).mean():.4f} min {(g.steps / g.opt).min():.4f}; "
                     f"LB1/OPT mean {(g.lb1 / g.opt).mean():.4f} min {(g.lb1 / g.opt).min():.4f}; LB/OPT mean {(g.lb / g.opt).mean():.4f} "
                     f"min {(g.lb / g.opt).min():.4f}; alpha/OPT mean {(g.alpha / g.opt).mean():.4f} (N={g.alpha.notna().sum()})")
    for n, g in u.groupby(level="n"):
        lines.append(f"  n={n} N={len(g)} LB/steps mean {(g.lb / g.steps).mean():.4f} LB1/steps {(g.lb1 / g.steps).mean():.4f}; "
                     f"sweep+reshape / steps {(g.ub / g.steps).mean():.4f} / LB1 {(g.ub / g.lb1).mean():.4f} / LB {(g.ub / g.lb).mean():.4f} "
                     f"(max {(g.ub / g.lb).max():.4f}); time sweep {L.tSweep[nn == n].mean():.1f} swaps {L.tSwap0[nn == n].mean():.1f} "
                     f"plateau {L.tSwap[nn == n].mean():.1f} passes0 max {L.passes0[nn == n].max()}")
    return lines


def time_to_near(d):
    """Series indexed by n: mean time for sweep + reshape to reach size <= NEAR * OPT, over the instances with proven OPT
    (eps = PAPER_EPS). Time = row time - lsTime (the sweep, the local-search setup and the final pruning) plus
    the first time of the trace (DATA/traces/RA_sweep_reshape_*, time since the start of the local search, logged at every
    size change) with size <= NEAR * OPT. If the sweep is already within NEAR * OPT, the time is that of the sweep alone
    (tBuild)."""
    s = d[(d.eps == PAPER_EPS) & (d.algo == "lssweep") & d.opt.notna()]
    out = []
    for r in s.itertuples():
        inst = f"*_{r.gen.replace('.instance', '')}_instance_0_0_1_20" if r.cgshop else f"{r.gen}_{r.n}_{r.seed}_1_20"
        f = trace_files(f"RA_sweep_reshape_{inst}")
        if len(f) != 1: raise SystemExit(f"trace of {r.gen} {r.n} {r.seed}: {f}")
        t = pd.read_csv(f[0], header=None, names=["t", "size"])
        hit = t[t["size"] <= NEAR * r.opt]
        if not len(hit): raise SystemExit(f"{r.gen} {r.n} {r.seed}: sweep + reshape never within {NEAR} OPT")
        if hit.t.iloc[0] == 0:           # the sweep is already within NEAR * OPT: only the sweep counts
            out.append((r.nn, float(re.search(r"tBuild=([0-9.]+)", r.info).group(1)))); continue
        lst = float(re.search(r"lsTime=([0-9.]+)", r.info).group(1))
        out.append((r.nn, r.time - lst + hit.t.iloc[0]))
    return pd.DataFrame(out, columns=["n", "t"]).groupby("n")["t"].mean()


def peak_memory():
    """Peak resident memory (GB) of the jobs in DATA/memory.csv: DataFrame (algo, n, rss). algo: lazyprune, ck, ckexact,
    lsck (reshape from CK exact), lssweep (reshape from the sweep); n is the nominal n."""
    m = read("memory.csv")
    if not len(m): return pd.DataFrame(columns=["algo", "n", "rss"])
    m = m[m.out.isin(["exp2", "exp3", "exp4", "ls_all_rs", "ls_all_big_5000"])]
    m = m.assign(algo=m.algo.replace({"ls-ck-reshape": "lsck", "ls-sweep-reshape": "lssweep"}))
    m = m[m.algo.isin([MAIN, "ck", "ckexact", "lsck", "lssweep"])]
    return pd.DataFrame({"algo": m.algo, "n": [nominal_n(g, n) for g, n in zip(m.gen, m.n)], "rss": m.rss_gb})


def time_fig(d, fname, size=(FIGW, 2.0)):
    """(a) Mean running time over all inputs for eps = PAPER_EPS, and the time for sweep + reshape to reach NEAR * OPT
    (instances with proven OPT only). (b) Largest peak memory over the inputs of each n >= 5000 (DATA/memory.csv)."""
    fig, (ax, bx) = plt.subplots(1, 2, figsize=size, gridspec_kw=dict(width_ratios=[1.25, 1]))
    for a in [MAIN, "lssweep"] + [a for a in ALGOS if a != MAIN] + ["opt"]:
        if a == "lssweep":
            s = time_to_near(d)
            ax.plot(s.index, s.values, label=f"Sweep + reshape\nto {NEAR:g}$\\,\\mathrm{{OPT}}$", **STYLE[a])
            continue
        s = d[(d.eps == PAPER_EPS) & (d.algo == a) & ~OPTX(d)].groupby("nn")["time"].mean()
        s = s[s > 0]
        if len(s): ax.plot(s.index, s.values, label={"opt": "Exact (CPLEX)"}.get(a, TABLE_ALGOS.get(a)), **STYLE[a])
        if a == "opt" and len(s):
            # n where CPLEX reached the time limit on at least one input: the mean is then only a lower bound
            to = d[(d.eps == PAPER_EPS) & (d.algo == a) & (d.optimal == False)]["nn"].unique()
            t = s[s.index.isin(to)]
            if len(t): ax.plot(t.index, t.values, ls="none", marker="v", color=METHOD["rose"], label="CPLEX time limit")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("$n$"); ax.set_ylabel("time (s)"); ax.grid(alpha=.3)
    ax.text(0.03, 0.95, "(a) time", transform=ax.transAxes, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"])
    mem = peak_memory()
    mem = mem[mem.n >= 5000]   # the log rounds to 0.1 GB
    names = {MAIN: "Sweep", "ck": "CK box", "ckexact": "CK", "lsck": "Reshape from CK",
             "lssweep": "Reshape from the sweep"}
    for a in [MAIN, "lssweep", "lsck", "ckexact"]:
        s = mem[mem.algo == a].groupby("n")["rss"].max()
        if len(s): bx.plot(s.index, s.values, label=names[a], **STYLE[a])
    bx.set_xscale("log"); bx.set_yscale("log"); bx.set_xlabel("$n$"); bx.set_ylabel("peak memory (GB)"); bx.grid(alpha=.3)
    bx.text(0.03, 0.95, "(b) memory", transform=bx.transAxes, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"])
    # (a) top to bottom, then the lines of (b) only, top to bottom in (b)
    h, l = ordered_handles([ax, bx], top_to_bottom(ax) + top_to_bottom(bx))
    fig.legend(h, l, loc="lower center", ncol=4)
    fig.tight_layout(rect=(0, 0.36 / size[1], 1, 1)); fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)
    return mem


# local search: algo names ls-<start>-<moves> -> ls<start> (reshape) or ls<start>lns (LNS); added to the batch of
# matching n, so that it gets the same reference as the other algorithms
ls = load("ls_all_rs.csv", "ls_all_big_5000.csv")
if len(ls):
    m = ls.algo.str.extract(r"^ls-(ck|sweep)-(reshape|lns)$")
    ls = ls[m[0].notna()].copy()
    ls["algo"] = "ls" + m[0] + m[1].map({"reshape": "", "lns": "lns"})
    ls["nn"] = [nominal_n(g, n) for g, n in zip(ls.gen, ls.n)]
parts = []
for fs, lo, hi in ((("exp2.csv", "exp2_opt.csv"), MIN_N, 1000), (("exp3.csv",), 1001, 10000), (("exp4.csv",), 10001, 10**9)):
    x = load(*fs)
    if len(ls): x = pd.concat([x, ls[(ls.nn >= lo) & (ls.nn <= hi)].drop(columns="nn")], ignore_index=True)
    if len(x): parts.append(with_ratios(x))
if not parts: raise SystemExit("no data")
d = pd.concat(parts, ignore_index=True)
d["nn"] = [nominal_n(g, n) for g, n in zip(d.gen, d.n)]
d = d[d.nn <= MAX_N]
d = d[~(d.cgshop & (d.family == "uniform"))]         # the CG:SHOP uniform family is not used
syn, cgs = d[~d.cgshop], d[d.cgshop]

CAP = (r"Mean size divided by the reference, per input: synthetic inputs (three seeds for $n\le10000$, one for"
       r" $n=50000$) and CG:SHOP 2019 families. $\dagger$: the reference is the lower bound. Bold: smallest"
       r" construction of the row.")
results_table(syn, cgs, "results.tex", CAP, "tab:exp")
cnt = summary_fig(d, "summary.pdf")
ratio_fig(d, "ratio.pdf")
mem = time_fig(d, "time.pdf")
lb_lines = lb_fig(d, "lb.pdf")
pp_lines = per_point_fig(d, "per_point.pdf")

# numbers for the text
lines = [f"instances per n: {cnt}"]
e = PAPER_EPS
for name, x in [("all", d), ("synthetic", syn), ("cgshop", cgs)]:
    for n in TABLE_NS:
        for a in list(ALGOS) + list(EXTRA):
            s = x[(x.algo == a) & (x.eps == e) & (x.nn == n)].dropna(subset=["ratio"])
            if len(s): lines.append(f"{name} n={n} {a}: ratio mean {s.ratio.mean():.3f} min {s.ratio.min():.3f} max {s.ratio.max():.3f} "
                                    f"(N={len(s)}, ref {'OPT' if not s.reflb.any() else 'LB' if s.reflb.all() else 'mixed'})")
    s = x[(x.algo == MAIN) & x.exact]
    if len(s): lines.append(f"{name} steps/OPT: min {(s.steps / s.opt).min():.4f} mean {(s.steps / s.opt).mean():.4f}; LB/OPT: min {(s.lb / s.opt).min():.4f} mean {(s.lb / s.opt).mean():.4f}")
for n in TABLE_NS:
    for gname, gs in [("per input", "gen")]:
        s = syn[(syn.algo.isin([MAIN, "ckexact", "lssweep"])) & (syn.nn == n)].groupby(["gen", "algo"]).ratio.mean().unstack()
        lines.append(f"synthetic n={n} per input:\n{s.round(3).to_string()}")
    s = cgs[(cgs.algo.isin([MAIN, "ckexact", "lssweep"])) & (cgs.nn == n)].groupby(["family", "algo"]).ratio.mean().unstack()
    if len(s): lines.append(f"cgshop n={n} per family:\n{s.round(3).to_string()}")
t = d[(d.eps == e) & ~OPTX(d)].groupby(["algo", "nn"]).time.agg(["mean", "median", "max"])
lines.append("times (s):\n" + t.loc[[a for a in [MAIN, "ck", "ckexact", "opt"] if a in t.index.get_level_values(0)]].round(3).to_string())
lines.append("peak memory (GB), max over inputs:\n" + mem.groupby(["algo", "n"]).rss.max().unstack(0).round(2).to_string())
lines += lb_lines + pp_lines
open(os.path.join(RESULTS, "tables", "summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines[:80]))

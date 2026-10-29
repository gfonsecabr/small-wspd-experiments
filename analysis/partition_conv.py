"""Partition WSPDs: conversion of cover WSPDs (src/partition.hpp) and the quadtree construction of thm:partition.
Run from the repository root:  python3 analysis/partition_conv.py   (DATA=runs ... to analyze your own runs)

Inputs (DATA, default data/; eps = 1/20). Every exp row of an algorithm is followed by the row <algo>+part of its
conversion (fewest/min/merge; the size without the merge pass is beforeMerge in the info column).
  partition_exp.csv      CK exact and the exact quadtrees (qtexact, qtexactrand), n = 200..10000; they are partitions,
                         so their conversion is the merge pass only
  partition_ls.csv       sweep + reshape, rerun with the number of rounds of ls_all_rs.csv (same WSPD), n <= 10000
  partition_ckls.csv, partition_ckls_big.csv    CK exact + reshape, as partition_ls.csv and partition_ls_big.csv
  partition_sweep.csv    the sweep (lazyprune), n = 200..10000 and n = 20000
  partition_exp_big.csv, partition_ls_big.csv   n = 20000: CK exact, the exact quadtrees, and the 5000 s sweep +
                         reshape WSPDs of ls_all_big_5000.csv (rows load-ls-sweep-reshape, renamed ls-sweep-reshape)
  partition_opt500.csv   the optimal WSPDs with n = 500 (rows load-opt in the paper's data, opt in new runs)
  exp2.csv ... exp4.csv, exp2_opt.csv, optx.csv, ls_all_rs.csv, ls_all_big_5000.csv   references and reproduction checks
No partition data for n > 20000: there the 5000 s local search is far from converged, so its WSPD keeps most of the
overlap of the sweep, and the conversion would measure that.
Reference of an instance: as in make_tables_figs.py. Since LB <= cover OPT <= partition OPT, size / reference of a
partition WSPD overestimates its ratio to the minimum partition WSPD.
Outputs (RESULTS, default results/): tables/partition.tex (Table partconv), figs/partition_conv.pdf (Fig. partconv),
partition_summary.txt (numbers for the text), partition_wide.csv (one row per instance)."""
import os, re, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import DATA, RESULTS, KEY, with_optx, proven_opt, sweep_lb, read as read_csv
from figutil import ordered_handles, top_to_bottom, METHOD, INPUT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = 0.05
SYN_GENS = ["uniform", "gauss", "cluster", "grid", "spiral", "annulus"]
LABEL = {"uniform": "uniform", "gauss": "Gaussian", "cluster": "clustered", "grid": "perturbed grid", "spiral": "spiral",
         "annulus": "annulus", "CG:SHOP": "CG:SHOP"}
INPUT_STYLE = INPUT
FIGW = 6.2
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.5, "lines.markersize": 2.2, "lines.linewidth": 0.9, "axes.linewidth": 0.6,
                     "grid.linewidth": 0.4, "legend.frameon": False})
OUT = []


def say(s=""):
    OUT.append(s); print(s)


def nominal_n(gen, n):
    m = re.search(r"-(\d{7})", gen)
    return int(m.group(1)) if gen.endswith(".instance") and m else int(n)


def read(*names):
    df = read_csv(*names)
    if not len(df): return df
    bad = df[~df.valid.astype(str).str.startswith("ok")]
    if len(bad): say(f"WARNING: {len(bad)} invalid rows\n{bad.head()}")
    df = df[df.eps == EPS].copy()
    df["cgshop"] = df.gen.str.endswith(".instance")
    df["family"] = [re.match(r"(.*?)-\d{7}", g).group(1) if c else g for g, c in zip(df.gen, df.cgshop)]
    df["inp"] = df.gen.where(~df.cgshop, "CG:SHOP")
    df["nn"] = [nominal_n(g, n) for g, n in zip(df.gen, df.n)]
    df["conv"] = df["info"].astype(str).str.extract(r"conv=([^;]+)")[0]
    df["cover"] = df["info"].astype(str).str.extract(r"cover=(\d+)")[0].astype(float)
    df["overlap"] = df["info"].astype(str).str.extract(r"overlap=([0-9.]+)")[0].astype(float)
    return df


def references():
    """DataFrame indexed by KEY: ref (OPT where proven, else the lower bound LB of reference.sweep_lb), reflb, and lb."""
    refs = []
    for fs in (("exp2.csv", "exp2_opt.csv"), ("exp3.csv",), ("exp4.csv",)):
        x = with_optx(read(*fs))
        if not len(x): continue
        st = x[x.algo == "lazyprune"].drop_duplicates(KEY, keep="last").set_index(KEY)
        lb = sweep_lb(st)
        o = x[x.algo == "opt"]
        best = proven_opt(o).reindex(lb.index) if len(o) else pd.Series(np.nan, index=lb.index)
        r = pd.DataFrame({"ref": best.fillna(lb), "reflb": best.isna(), "lb": lb})
        refs.append(r)
    return pd.concat(refs)


# ---------------------------------------------------------------------------------------------- main runs
DEFAULT = "fewest/min/merge"
NOMERGE = "fewest/min/nomerge"
# per-instance columns: (label, algo, conv) with conv None for the algorithm's own output
COLS = {"rs": ("ls-sweep-reshape", None), "rs_p": ("ls-sweep-reshape", DEFAULT), "rs_pn": ("ls-sweep-reshape", NOMERGE),
        "rc": ("ls-ck-reshape", None), "rc_p": ("ls-ck-reshape", DEFAULT), "rc_pn": ("ls-ck-reshape", NOMERGE),
        "sw": ("lazyprune", None), "sw_p": ("lazyprune", DEFAULT), "sw_pn": ("lazyprune", NOMERGE),
        "opt": ("opt", None), "opt_p": ("opt", DEFAULT), "opt_pn": ("opt", NOMERGE),
        "cke": ("ckexact", None), "cke_m": ("ckexact", DEFAULT), "qte": ("qtexact", None), "qte_m": ("qtexact", DEFAULT),
        "qter": ("qtexactrand", None), "qter_m": ("qtexactrand", DEFAULT)}
NAMES = {"rs": "Sweep + reshape (cover)", "rs_p": "Sweep + reshape, converted", "rc": "CK + reshape (cover)",
         "rc_p": "CK + reshape, converted", "sw": "Sweep (cover)", "sw_p": "Sweep, converted",
         "opt_p": "OPT, converted", "cke": "CK", "cke_m": "CK + merge", "qte": "Quadtree (exact)",
         "qter": "Quadtree (exact, random)", "qte_m": "Quadtree (exact) + merge", "qter_m": "Quadtree (exact, random) + merge"}
TABLE_NS = [1000, 10000]
MAX_N = 20000   # no partition data for n > 20000 (see the docstring)


def load_runs():
    """All rows of the partition runs with n <= MAX_N (rows loaded from dumps renamed: load-<algo> -> <algo>)."""
    e = read("partition_opt500.csv", "partition_exp.csv", "partition_ls.csv", "partition_sweep.csv",
             "partition_exp_big.csv", "partition_ls_big.csv", "partition_ckls.csv", "partition_ckls_big.csv")
    if not len(e): return e
    e["algo"] = e.algo.str.replace(r"^load-", "", regex=True)
    e["conv"] = e["conv"].fillna("")
    e = e[e.nn <= MAX_N]
    return e.drop_duplicates(KEY + ["algo", "conv"], keep="last")


def wide():
    """One row per instance (KEY) with the sizes of COLS, the overlaps, the conversion times, and the reference."""
    e = load_runs()
    if not len(e): return None, e
    base = e[~e.algo.str.endswith("+part")]
    part = e[e.algo.str.endswith("+part")].assign(src=lambda x: x.algo.str[:-5])
    W = base.drop_duplicates(KEY).set_index(KEY)[["cgshop", "family", "inp", "nn"]].copy()
    for c, (a, conv) in COLS.items():
        s = (base[base.algo == a] if conv is None else part[(part.src == a) & (part.conv == conv)]).set_index(KEY)
        W[c] = s["size"].reindex(W.index)
        if conv == NOMERGE:   # runs with the merged conversion only: the size before the merge pass (beforeMerge)
            m = part[(part.src == a) & (part.conv == DEFAULT)].set_index(KEY)
            W[c] = W[c].fillna(m["info"].str.extract(r"beforeMerge=(\d+)")[0].astype(float).reindex(W.index))
        W[c + "_t"] = s["time"].reindex(W.index)
        if conv is not None:
            W[c + "_ov"] = s["overlap"].reindex(W.index)
    R = references()
    W = W.join(R, how="left")
    return W, e


def checks(e):
    say("== reproduction checks (sizes equal to the paper's runs)")
    ref_runs = read("exp2.csv", "exp2_opt.csv", "exp3.csv", "exp4.csv", "ls_all_rs.csv", "ls_all_big_5000.csv")
    base = e[~e.algo.str.endswith("+part")]
    for a in ("lazyprune", "ckexact", "ls-sweep-reshape", "ls-ck-reshape", "opt"):
        x = base[base.algo == a].set_index(KEY)["size"]
        y = ref_runs[ref_runs.algo == a].drop_duplicates(KEY, keep="last").set_index(KEY)["size"].reindex(x.index)
        same = (x == y)[y.notna()]
        say(f"{a}: {len(x)} rows, {int(same.sum())} identical, {int(y.isna().sum())} without a paper row"
            + ("" if same.all() else f"  DIFFERENT: {list(x[(x != y) & y.notna()].index[:5])}"))
    for a in ("ls-sweep-reshape", "ls-ck-reshape"):   # the reproduced runs also match in moves and kills
        if a not in set(base.algo): continue
        x = base[base.algo == a].set_index(KEY)["info"].str.extract(r"(rounds=\d+;reshapeMoves=\d+;reshapeKills=\d+)")[0]
        y = ref_runs[ref_runs.algo == a].drop_duplicates(KEY, keep="last").set_index(KEY)["info"].str.extract(
            r"(rounds=\d+;reshapeMoves=\d+;reshapeKills=\d+)")[0].reindex(x.index)
        x = x.dropna()   # the rows loaded from dumps (n >= 20000) have no counters; their sizes are checked above
        say(f"{a}: rounds, moves, kills identical on {int((x == y.reindex(x.index)).sum())} of {len(x)} reruns")
    q = base[base.algo.str.startswith("qt")]
    if len(q): say(f"quadtree rows: {len(q)}, partition verified on {int(q['info'].str.contains('partition=1').sum())}")


def kn(n):
    """n in the tables: 1000 -> 1k, 50000 -> 50k (text mode)."""
    return f"{n // 1000}k" if n >= 1000 and n % 1000 == 0 else str(n)


def part_table(W, fname):
    """Appendix table: rows (input, n), columns: size/reference."""
    cols = [("rs", r"Sweep\\+ reshape"), ("rs_pn", "Converted"), ("rs_p", r"Converted\\+ merge"), ("cke", "CK"), ("cke_m", r"CK\\+ merge"),
            ("qte_m", r"Quadtree\\+ merge")]
    out = [r"\begin{table}[htbp]", r"\centering\footnotesize\setlength{\tabcolsep}{4pt}",
           r"\caption{Partition WSPDs: mean size divided by the reference ($\dagger$: the lower bound) per input and $n$."
           r" First column: the cover WSPD of the sweep followed by reshape. The others are partition WSPDs: its conversion"
           r" and CK, each without and with the merge pass, and the deterministic quadtree with exact separation"
           r" tests and the merge pass. Bold: smallest partition WSPD of the row.}",
           r"\label{tab:partconv}",
           r"\begin{tabular}{@{}lr" + "c" * len(cols) + "@{}}", r"\toprule",
           r"Input & $n$ & " + " & ".join(r"\shortstack{\strut " + h.replace(r"\\", r"\\\strut ") + "}" for _, h in cols) + r" \\",
           r"\midrule"]
    # CG:SHOP families with several sizes first (by name), then those with one size, by that size (then name)
    cg = W[W.cgshop & W.nn.isin(TABLE_NS)]
    sizes = cg.groupby("family").nn.nunique()
    multi = sorted(f for f in sizes.index if sizes[f] > 1)
    single = sorted((f for f in sizes.index if sizes[f] == 1), key=lambda f: (cg[cg.family == f].nn.min(), f))
    rows = [(g, False) for g in SYN_GENS] + [(f, True) for f in multi + single]
    for i, (r, cg) in enumerate(rows):
        sub = W[(W.family == r) & (W.cgshop == cg)]
        ns = [n for n in TABLE_NS if (sub.nn == n).any()]
        for k, n in enumerate(ns):
            x = sub[sub.nn == n]
            mark = r"$^\dagger$" if x.reflb.all() else r"$^\ddagger$" if x.reflb.any() else r"\phantom{$^\dagger$}"
            cells = []
            vals = {c: (x[c] / x.ref).mean() for c, _ in cols}
            best = min(vals[c] for c, _ in cols[1:])
            cells += [(rf"\textbf{{{vals[c]:.2f}}}" if c != "rs" and f"{vals[c]:.2f}" == f"{best:.2f}" else f"{vals[c]:.2f}") for c, _ in cols]
            out.append((LABEL.get(r, r) if k == 0 else "") + f" & {kn(n)}{mark} & " + " & ".join(cells) + r" \\")
        if i < len(rows) - 1 and ns: out.append(r"\addlinespace[2pt]")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(RESULTS, "tables", fname), "w").write("\n".join(out) + "\n")


def refmarks(ax, g, colr, fill_marker="o", **kw):
    """Markers as in fig:ratio (make_tables_figs.ratio_fig): filled where the reference is OPT for all instances, open
    where it is the lower bound for all, half-filled where mixed. g: DataFrame indexed by n with columns r, lb, mix."""
    ax.plot(g.index[g.lb], g.r[g.lb], ls="none", marker=fill_marker, mfc="white", mec=colr, mew=0.7, **kw)
    m = g.mix & ~g.lb
    if m.any(): ax.plot(g.index[m], g.r[m], ls="none", marker=fill_marker, fillstyle="left", mfc=colr, mfcalt="white", mec=colr, mew=0.7, **kw)


def panel_title(ax, t):   # title inside the panel, as in fig:ratio
    ax.text(0.03, 0.95, t, transform=ax.transAxes, ha="left", va="top", fontsize=plt.rcParams["axes.titlesize"],
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))


def part_fig(W, fname, size=(FIGW, 5.0)):
    """Top: size/reference vs n per input, (a) sweep + reshape (cover), (b) its conversion with the merge pass, same y
    scale (as fig:ratio). Bottom, full width: (c) mean size/reference vs n of the cover and the partition WSPDs."""
    fig = plt.figure(figsize=size)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.45], hspace=0.5, wspace=0.06)
    ax = fig.add_subplot(gs[0, 0]); bx = fig.add_subplot(gs[0, 1], sharey=ax); cx = fig.add_subplot(gs[1, :])
    for p_, col, t in ((ax, "rs", "(a) Sweep + reshape (cover)"), (bx, "rs_p", "(b) Converted + merge (partition)")):
        for inp, colr in INPUT_STYLE.items():
            x = W[W.inp == inp]
            g = pd.DataFrame({"r": (x[col] / x.ref), "lb": x.reflb.astype(bool), "nn": x.nn}).dropna(subset=["r"]).groupby("nn").agg(
                r=("r", "mean"), lb=("lb", "all"), mix=("lb", "any"))
            if not len(g): continue
            p_.plot(g.index, g.r, color=colr, marker="o", label=LABEL[inp])
            refmarks(p_, g, colr)
        p_.set_xscale("log"); p_.grid(alpha=.3); p_.set_xlabel("$n$")
        panel_title(p_, t)
    ax.set_ylabel("size / reference"); plt.setp(bx.get_yticklabels(), visible=False)
    h, l = ordered_handles([ax, bx])
    fig.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.08)
    y = (ax.get_position().y0 - 0.06 + cx.get_position().y1) / 2   # between the x label of (a), (b) and panel (c)
    fig.legend(h, l, loc="center", ncol=len(h), bbox_to_anchor=(0.5, y))
    lines = [("rs", "Sweep + reshape (cover)", dict(color=METHOD["magenta"], ls=":", marker="x")),
             ("sw_pn", "Sweep, converted", dict(color=METHOD["blue"], ls="--", marker="^")),
             ("sw_p", "Sweep, converted + merge", dict(color=METHOD["blue"], ls="-", marker="^")),
             # CK exact + reshape is within 0.4% of sweep + reshape on average for every n: drawn below it, with larger markers
             ("rc_pn", "CK + reshape, converted", dict(color=METHOD["pink"], ls="--", marker="v", ms=3.8)),
             ("rc_p", "CK + reshape, converted + merge", dict(color=METHOD["pink"], ls="-", marker="v", ms=3.8)),
             ("rs_pn", "Sweep + reshape, converted", dict(color=METHOD["magenta"], ls="--", marker="o")),
             ("rs_p", "Sweep + reshape, converted + merge", dict(color=METHOD["magenta"], ls="-", marker="o")),
             ("cke", "CK", dict(color=METHOD["green"], ls="--", marker="D")),
             ("cke_m", "CK + merge", dict(color=METHOD["green"], ls="-", marker="D")),
             ("qte_m", "Quadtree + merge", dict(color=METHOD["gold"], ls="-", marker="s"))]
    for c, lab, st in lines:
        g = pd.DataFrame({"r": W[c] / W.ref, "lb": W.reflb.astype(bool), "nn": W.nn}).dropna(subset=["r"]).groupby("nn").agg(
            r=("r", "mean"), lb=("lb", "all"), mix=("lb", "any"))
        if not len(g): continue
        cx.plot(g.index, g.r, label=lab, **st)
        if st["marker"] != "x": refmarks(cx, g, st["color"], st["marker"], **({"ms": st["ms"]} if "ms" in st else {}))
    cx.set_xscale("log"); cx.grid(alpha=.3); cx.set_xlabel("$n$"); cx.set_ylabel("size / reference")
    cx.legend(*ordered_handles(cx), loc="upper left", bbox_to_anchor=(0.0, 0.88))
    panel_title(cx, "(c) Mean over all inputs")
    fig.savefig(os.path.join(RESULTS, "figs", fname)); plt.close(fig)


def numbers(W):
    say("== numbers for the text")
    say(f"instances: {len(W)}; per n: {W.groupby('nn').size().to_dict()}")
    for c, lab in (("rs", "sweep + reshape"), ("rc", "CK exact + reshape"), ("sw", "sweep"), ("opt", "OPT")):
        for conv in ("_p", "_pn"):
            o = (W[c + conv] / W.ref).dropna()
            if not len(o): continue
            g = o.groupby(W.nn).agg(["mean", "max", "min", "count"])
            say(f"converted / reference, {lab}{' (no merge)' if conv == '_pn' else ''}: all mean {o.mean():.4f} max {o.max():.4f} min {o.min():.4f}")
            say("   by n (mean/max): " + "; ".join(f"{n}: {r['mean']:.4f}/{r['max']:.4f}" for n, r in g.iterrows()))
            say(f"   overlap mean by n: " + "; ".join(f"{n}: {v:.3f}" for n, v in W.groupby('nn')[c + conv + '_ov'].mean().items()))
    o = (W.opt_p / W.opt).dropna()
    if len(o):
        say(f"converted OPT / OPT (an upper bound on partition OPT / cover OPT): mean {o.mean():.4f}, max {o.max():.4f} "
            f"({W.loc[o.idxmax()][['family', 'nn']].to_list()}), = 1 on {int((o == 1).sum())} of {len(o)}")
    for c in ("rs", "rs_p", "rc", "rc_p", "sw", "sw_p", "cke", "cke_m", "qte", "qte_m", "qter", "qter_m"):
        r = (W[c] / W.ref).groupby(W.nn).mean()
        say(f"size/ref {NAMES[c]:34s} " + "  ".join(f"{n}:{v:.3f}" for n, v in r.items() if n in (200, 500, 1000, 2000, 5000, 10000, 20000)))
    for c in ("cke", "cke_m", "qte_m", "qter_m"):
        x = W.dropna(subset=["rs_p", c])
        say(f"converted sweep + reshape <= {NAMES[c]} on {int((x.rs_p <= x[c]).sum())} of {len(x)} instances; by n (count larger): "
            + "; ".join(f"{n}: {int(v)}" for n, v in (x.rs_p > x[c]).groupby(x.nn).sum().items() if v))
    x = W.dropna(subset=["rs_p", "rc_p"])
    for a, b, lab in (("rc", "rs", "cover"), ("rc_pn", "rs_pn", "converted"), ("rc_p", "rs_p", "converted + merge")):
        say(f"CK exact + reshape vs sweep + reshape, {lab}: smaller on {int((x[a] < x[b]).sum())}, equal on {int((x[a] == x[b]).sum())} of {len(x)} "
            f"instances; ratio mean by n " + "; ".join(f"{n}: {v:.3f}" for n, v in (x[a] / x[b]).groupby(x.nn).mean().items()))
    say(f"converted sweep + reshape vs CK exact: ratio mean by n " + "; ".join(f"{n}: {v:.3f}" for n, v in (W.rs_p / W.cke).groupby(W.nn).mean().items()))
    say(f"converted sweep + reshape vs CK exact + merge: ratio mean by n " + "; ".join(f"{n}: {v:.3f}" for n, v in (W.rs_p / W.cke_m).groupby(W.nn).mean().items()))
    for c in ("rs_p", "rc_p", "sw_p", "cke_m", "qter_m"):
        t = W[c + "_t"].groupby(W.nn).agg(["mean", "max"])
        say(f"time of {NAMES[c]} (conversion or merge only): " + "; ".join(f"{n}: {r['mean']:.2f}/{r['max']:.1f}s" for n, r in t.iterrows() if n in (1000, 10000, 20000)))
    for c in ("qte", "qter", "cke"):
        t = W[c + "_t"].groupby(W.nn).mean()
        say(f"time of {NAMES[c]}: " + "; ".join(f"{n}: {v:.3f}s" for n, v in t.items() if n in (1000, 10000, 20000)))
    # per input at n = 10000 and 20000
    for nn in (10000, 20000):
        x = W[W.nn == nn]
        if len(x):
            say(f"n = {nn} per input: overlap; size/ref of the cover, of the converted WSPD, of CK exact, of CK exact + merge")
            say(x.groupby("inp").apply(lambda g: pd.Series({"overlap": g.rs_p_ov.mean(), "cover": (g.rs / g.ref).mean(),
                                                            "conv": (g.rs_p / g.ref).mean(), "cke": (g.cke / g.ref).mean(),
                                                            "cke_m": (g.cke_m / g.ref).mean()}), include_groups=False).round(3).to_string())
    c = W.rs_p_ov.astype(float).corr((W.rs_p / W.rs).astype(float))
    say(f"correlation overlap vs partition / cover (sweep + reshape): {c:.3f}")


def main_runs():
    W, e = wide()
    if W is None: say("no main data"); return
    checks(e)
    W.reset_index().to_csv(os.path.join(RESULTS, "partition_wide.csv"), index=False)
    part_table(W.dropna(subset=["rs", "rs_p"]), "partition.tex")
    part_fig(W, "partition_conv.pdf")
    numbers(W)


if __name__ == "__main__":
    main_runs()
    open(os.path.join(RESULTS, "partition_summary.txt"), "w").write("\n".join(OUT) + "\n")

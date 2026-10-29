"""The numbers quoted in the text of the paper (abstract, introduction, and Section Experiments), computed from the data.
Run from the repository root:  python3 analysis/text_numbers.py   (DATA=runs ... to analyze your own runs)
Output: RESULTS/text_numbers.txt (and printed). Each line names the statement it supports; the paper rounds.
Reference of an instance: OPT where proven, else the lower bound LB (lb_check.csv, reference.py)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import KEY, MIN_N, MAX_N, RESULTS, read, refs, trace_files, with_optx, lower_bounds

BIG = MAX_N   # the largest n

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def nominal(df):
    m = df.gen.str.extract(r"-(\d{7})\.instance$")[0]
    return m.astype(float).fillna(df.n).astype(int)


def pct(x):
    return f"{100 * (x - 1):.3f}%"


def runs():
    """All construction and local-search rows (eps = 1/20) with nn, opt, ref, ratio, and the LB."""
    b = with_optx(read("exp2.csv", "exp2_opt.csv", "exp3.csv", "exp4.csv"))
    ls = read("ls_all_rs.csv", "ls_all_big_5000.csv")
    d = pd.concat([b, ls], ignore_index=True)
    d = d[d.valid == "ok"].drop_duplicates(KEY + ["algo"], keep="last").copy()
    opt, ref = refs()
    steps = d[d.algo == "lazyprune"].set_index(KEY)["info"].str.extract(r"steps=(\d+)")[0].astype(float).rename("steps")
    d = d.join(opt, on=KEY).join(ref.rename("ref"), on=KEY).join(steps, on=KEY)
    d["nn"] = nominal(d)
    L = lower_bounds()[["lb1", "lb"]]
    d = d.join(L, on=["gen", "nn", "seed", "eps"])
    d["nn"] = nominal(d)
    d["cg"] = d.gen.str.endswith(".instance")
    d["r"] = d["size"] / d["ref"]
    d["ropt"] = d["size"] / d["opt"]
    return d


def main():
    d = runs()
    A = lambda a: d[d.algo == a]
    sw, cke, ckb, opt = A("lazyprune"), A("ckexact"), A("ck"), A("opt")
    rs, rsck = A("ls-sweep-reshape"), A("ls-ck-reshape")
    lns_sw, lns_ck = A("ls-sweep-lns"), A("ls-ck-lns")
    lbopt = sw.dropna(subset=["opt"]).assign(q=lambda x: x.lb / x.opt, q0=lambda x: x.steps / x.opt, q1=lambda x: x.lb1 / x.opt)

    say("== Reference value")
    o = opt.assign(proven=opt["info"].str.contains("optimal=1"))
    o = o[~o["info"].str.contains("solver=optx")]
    g = o.groupby("nn").proven.agg(["sum", "count"])
    say("CPLEX proves optimality on all instances with n <= 500; for 600 <= n <= 1000 on some (proven/runs by n; the"
        " runs that ran out of memory have no row): " + "; ".join(f"{n}: {int(r['sum'])}/{int(r['count'])}" for n, r in g.iterrows()))
    say(f"optima proven by optx (spiral, n = 800..1000): {int((opt['info'].str.contains('solver=optx')).sum())}")
    say("instances per n: " + ", ".join(f"{n}: {c}" for n, c in sw.groupby("nn").size().items()))
    for n in (200, 500, 1000):
        x = lbopt[lbopt.nn == n]
        say(f"n = {n} ({len(x)} instances with OPT): LB/OPT mean {x.q.mean():.4f}, min {x.q.min():.4f}; LB1/OPT mean {x.q1.mean():.4f}, "
            f"min {x.q1.min():.4f}; LB0 (steps)/OPT mean {x.q0.mean():.4f}, min {x.q0.min():.4f}")
    say(f"over all instances with OPT ({len(lbopt)}): LB/OPT min {lbopt.q.min():.4f}, mean {lbopt.q.mean():.4f}; LB0/OPT min {lbopt.q0.min():.4f}; "
        f"LB = OPT on {int((lbopt.lb == lbopt.opt).sum())}")
    say("\n== Lower bounds (Section Lower bounds, Fig. lb)")
    L = lower_bounds()
    L = L[L.index.get_level_values("n") >= MIN_N]
    nn = L.index.get_level_values("n")
    a = read("alpha.csv")
    if len(a):
        a = a[(a.optimal == 1) & (a.check == "ok")].assign(eps=0.05).set_index(KEY)
        j = lbopt.set_index(KEY).join(a.alpha)
        j = j[j.alpha.notna()]
        say(f"alpha (n = 200..500, {len(j)} instances): alpha/OPT mean {(j.alpha / j.opt).mean():.4f}, min {(j.alpha / j.opt).min():.4f}, "
            f"alpha = OPT on {int((j.alpha == j.opt).sum())}; LB = alpha on {int((j.lb == j.alpha).sum())}; of the gap alpha - LB0, "
            f"LB1 closes {(j.lb1 - j.steps).sum() / (j.alpha - j.steps).sum():.3f} and LB {(j.lb - j.steps).sum() / (j.alpha - j.steps).sum():.3f}")
    for n in (1000, 10000, 20000, BIG):
        x = L[nn == n]
        say(f"n = {n}: LB/LB0 mean {(x.lb / x.steps).mean():.4f} (min {(x.lb / x.steps).min():.4f}, max {(x.lb / x.steps).max():.4f}), "
            f"LB1/LB0 mean {(x.lb1 / x.steps).mean():.4f}; time: sweep {x.tSweep.mean():.1f} s, lb_swap to LB1 {x.tSwap0.mean():.1f} s, "
            f"to LB {x.tSwap.mean():.1f} s (max {x.tSwap.max():.1f} s); passes to LB1 at most {x.passes0.max()}")
    say(f"passes without plateau moves: {L.passes0.min()} to {L.passes0.max()}")
    for n in (10000, BIG):
        x = rs[rs.nn == n]
        say(f"n = {n}, sweep + reshape: mean {x.r.mean():.3f} LB (max {x.r.max():.3f}), {(x['size'] / x.lb1).mean():.3f} LB1 (max {(x['size'] / x.lb1).max():.3f}), "
            f"{(x['size'] / x.steps).mean():.3f} LB0 (max {(x['size'] / x.steps).max():.3f})")

    say("\n== Quality (Q1)")
    for name, x in (("sweep", sw), ("CK box", ckb), ("CK exact", cke)):
        y = x[x.nn == 500]
        say(f"n = 500, {name}: size/OPT mean {y.ropt.mean():.4f}, max {y.ropt.max():.4f}")
    for n in (10000, BIG):
        say(f"n = {n}: " + ", ".join(f"{name} mean {x[x.nn == n].r.mean():.3f} max {x[x.nn == n].r.max():.3f}"
                                     for name, x in (("sweep", sw), ("CK box", ckb), ("CK exact", cke))) + " (x LB)")
    u = sw[(sw.gen == "uniform") & (sw.nn == 1000)]
    say(f"uniform n = 1000, sweep: {u.ropt.mean():.3f} OPT, {(u['size'] / u.lb).mean():.3f} LB, {(u['size'] / u.steps).mean():.3f} LB0")
    for n in (1000, 10000, BIG):
        say(f"spiral n = {n}: CK exact {cke[(cke.gen == 'spiral') & (cke.nn == n)].r.mean():.3f}, sweep {sw[(sw.gen == 'spiral') & (sw.nn == n)].r.mean():.3f}")
    g = d[(d.gen == "grid") & (d.nn == 1000) & d.algo.isin(["lazyprune", "ck", "ckexact"])].groupby("algo").r.mean()
    say(f"perturbed grid n = 1000: " + ", ".join(f"{a} {v:.3f}" for a, v in g.items()) + f" (max/min {g.max() / g.min():.3f})")

    say("\n== Running times and memory (Q2)")
    t = d[d.algo.isin(["lazyprune", "ck", "ckexact"])].groupby(["algo", "nn"]).time.agg(["mean", "max"])
    for a in ("ck", "ckexact"):
        say(f"{a}: n = 10000 mean {t.loc[(a, 10000), 'mean']:.2f} s, max {t.loc[(a, 10000), 'max']:.2f} s; n = {BIG} mean {t.loc[(a, BIG), 'mean']:.1f} s, max {t.loc[(a, BIG), 'max']:.1f} s")
    say(f"sweep mean: n = 1000 {t.loc[('lazyprune', 1000), 'mean']:.2f} s, n = 10000 {t.loc[('lazyprune', 10000), 'mean']:.1f} s, "
        f"n = {BIG} {t.loc[('lazyprune', BIG), 'mean'] / 60:.1f} min")
    m = read("memory.csv")
    if len(m):
        m = m.assign(nn=nominal(m))
        m = m[(m.nn == BIG) & m.out.isin(["exp4", "ls_all_big_5000"])].groupby("algo").rss_gb.max()
        say(f"peak memory n = {BIG} (GB): " + ", ".join(f"{a} {v:.1f}" for a, v in m.items()))
    p = o[o.proven]
    for n in (200, 1000):
        say(f"OPT time (enumeration + CPLEX), median where proven, n = {n}: {p[p.nn == n].time.median():.2f} s ({int((p.nn == n).sum())} instances)")
    cols = p["info"].str.extract(r"cols=(\d+)")[0].astype(float)
    say(f"maximal dumbbells (set-cover columns) where CPLEX proves OPT: median for n = 1000 {cols[p.nn == 1000].median():.0f}, largest {cols.max():.0f}")

    say("\n== Local search (Q3)")
    for st, R, L in (("ck", rsck, lns_ck), ("sweep", rs, lns_sw)):
        j = R.set_index(KEY)[["size", "nn"]].join(L.set_index(KEY)["size"].rename("lns"), how="inner")
        say(f"from {st}, ls_all_rs (n <= 1000): reshape <= LNS on {int((j['size'] <= j.lns).sum())} of {len(j)} (equal on {int((j['size'] == j.lns).sum())}); "
            "mean reshape/LNS by n " + ", ".join(f"{n}: {v:.4f}" for n, v in (j['size'] / j.lns).groupby(j.nn).mean().items()))
    c = read("ls_cmp_rs.csv", "ls_cmp_big.csv")   # the samples; n = BIG: LNS only, with the budget of ls_all_big_5000
    big = c[nominal(c) == BIG][KEY].drop_duplicates()
    b5 = read("ls_all_big_5000.csv")
    c = pd.concat([c, b5[b5.algo.str.endswith("-reshape")].merge(big, on=KEY)], ignore_index=True)   # their reshape runs
    c = c[c.valid == "ok"].assign(nn=lambda x: nominal(x))
    c[["start", "moves"]] = c.algo.str.extract(r"^ls-(ck|sweep)-(reshape|lns)$")
    w = c.pivot_table(index=KEY + ["nn", "start"], columns="moves", values="size").reset_index()
    say(f"ls_cmp (samples, n = 1000, 10000, {BIG}): reshape <= LNS on "
        f"{int((w.reshape <= w.lns).sum())} of {len(w)}; mean reshape/LNS by n "
        + ", ".join(f"{n}: {v:.4f}" for n, v in (w.reshape / w.lns).groupby(w.nn).mean().items()))
    budget = {1000: 100, 10000: 1000, BIG: 5000}
    reach = []
    for r in w.itertuples():
        tag = f"images_{r.gen.removesuffix('.instance')}_instance_0_0_1_20" if r.gen.endswith(".instance") else f"{r.gen}_{r.n}_{r.seed}_1_20"
        f = trace_files(f"{'RL' if r.nn == BIG else 'R'}_{r.start}_reshape_{tag}")
        if len(f) != 1: continue
        tr = pd.read_csv(f[0], header=None, names=["t", "size"])
        h = tr[tr["size"] <= r.lns]
        if len(h): reach.append((r.nn, h.t.iloc[0], h.t.iloc[0] / budget[r.nn]))
    reach = pd.DataFrame(reach, columns=["nn", "t", "frac"])
    say(f"time for reshape to reach the final size of LNS (same start): max {reach.t.max():.2f} s, max fraction of the budget "
        f"{reach.frac.max():.4f} ({len(reach)} runs with traces); by n: "
        + ", ".join(f"{n}: max {g.t.max():.2f} s, {100 * g.frac.max():.2f}%" for n, g in reach.groupby("nn")))
    ref = refs()[1]
    for n in (10000, BIG):
        for st in ("ck", "sweep"):
            x = c[(c.nn == n) & (c.moves == "lns") & (c.start == st)].join(ref.rename("ref"), on=KEY)
            s0 = x["info"].str.extract(r"start=(\d+)")[0].astype(float)
            calls = x["info"].str.extract(r"lnsCalls=(\d+)")[0].astype(float)
            say(f"LNS from {st}, n = {n}: start {(s0 / x.ref).mean():.3f}, end {(x['size'] / x.ref).mean():.3f} (x reference, {len(x)} instances; "
                f"improvement {100 * (1 - x['size'] / s0).min():.2f}% to {100 * (1 - x['size'] / s0).max():.2f}%); "
                f"LNS moves {calls.min():.0f} to {calls.max():.0f}, dumbbells {s0.min() / 1e6:.1f} to {s0.max() / 1e6:.1f} million")
    x = rs.dropna(subset=["opt"])
    x = x[x.nn <= 1000]
    say(f"sweep + reshape, {len(x)} instances with n <= 1000 and OPT: mean {pct(x.ropt.mean())} above OPT, max {pct(x.ropt.max())}, "
        f"= OPT on {int((x['size'] == x.opt).sum())}; by n (= OPT / with OPT): "
        + ", ".join(f"{n}: {int((g['size'] == g.opt).sum())}/{len(g)}" for n, g in x.groupby("nn")))
    y = lns_sw.dropna(subset=["opt"])
    y = y[y.nn <= 1000]
    say(f"LNS from the sweep, same instances: mean {pct(y.ropt.mean())} above OPT, max {pct(y.ropt.max())}, = OPT on {int((y['size'] == y.opt).sum())}")
    for n in (500,):
        z = x[x.nn == n]
        say(f"n = {n}: sweep + reshape max {pct(z.ropt.max())} above OPT, = OPT on {int((z['size'] == z.opt).sum())} of {len(z)}")
    for n in (10000, BIG):
        say(f"n = {n}: sweep + reshape mean {rs[rs.nn == n].r.mean():.3f} LB, max {rs[rs.nn == n].r.max():.3f}; sweep alone mean {sw[sw.nn == n].r.mean():.3f}")
    j = rs.set_index(KEY)[["size", "nn"]].join(rsck.set_index(KEY)["size"].rename("ck"), how="inner")
    q = (j[["size", "ck"]].max(axis=1) / j[["size", "ck"]].min(axis=1) - 1)
    say(f"two starts of reshape end within {100 * q[j.nn <= 10000].max():.2f}% (n <= 10000) and {100 * q[j.nn >= 20000].max():.2f}% (n >= 20000)")
    e = j.ck / j["size"]
    for lo, hi in ((MIN_N, 1000), (2000, 10000), (20000, BIG)):
        x = e[(j.nn >= lo) & (j.nn <= hi)]
        say(f"reshape from CK exact vs from the sweep, n = {lo}..{hi}: ends smaller on {int((x < 1).sum())}, equal on {int((x == 1).sum())}, "
            f"larger on {int((x > 1).sum())} of {len(x)}; at most {100 * (1 - x.min()):.2f}% smaller, {100 * (x.max() - 1):.2f}% larger")
    say("mean final size from CK exact / from the sweep: " + ", ".join(f"n = {n}: {e[j.nn == n].mean():.4f}" for n in (1000, 10000, BIG)))
    cu = read("catchup.csv")
    if len(cu):
        cu = cu.assign(nn=nominal(cu), frac=lambda x: x.total_s / x.sweep_s)
        never = cu[cu.total_s.isna()]
        say(f"reshape from CK exact is smaller than the sweep's output before the sweep ends (CK + setup + reshape time < sweep time) on "
            f"{int((cu.frac < 1).sum())} of {len(cu)} instances; never smaller on {len(never)} (n = {sorted(set(never.nn))}, {sorted(set(never.gen))}); "
            f"smallest n from which it holds on every instance: {int(cu[~(cu.frac < 1)].nn.max()) + 100 if (~(cu.frac < 1)).any() else MIN_N}")
        say("  time until then / sweep time: " + "; ".join(
            f"n = {n}: mean {cu[cu.nn == n].frac.mean():.2f}, max {cu[cu.nn == n].frac.max():.2f}" for n in (1000, 10000, BIG))
            + f"; max over n >= 300: {cu[cu.nn >= 300].frac.max():.2f}")
        if "sweep_ahead" in cu:
            a = cu.sweep_ahead
            say(f"at equal total time (build + setup + reshape), reshape from the sweep is at most {100 * a.max():.2f}% smaller than reshape from "
                f"CK exact ({len(cu)} instances; it is never smaller for n >= {int(cu[a > 0].nn.max()) + 10000 if (a > 0).any() else MIN_N})")

    say("\n== Partition WSPDs (Q4), abstract, and introduction")
    import partition_conv as pc
    W, _ = pc.wide()
    for nn in (500, 1000, 10000, 20000):
        x = W[W.nn == nn]
        say(f"n = {nn}: converted {(x.rs_p / x.ref).mean():.4f} (max {(x.rs_p / x.ref).max():.3f}), cover {(x.rs / x.ref).mean():.4f}, "
            f"CK exact {(x.cke / x.ref).mean():.3f}, + merge {(x.cke_m / x.ref).mean():.3f}, quadtree {(x.qte / x.ref).mean():.3f}, "
            f"+ merge {(x.qte_m / x.ref).mean():.3f}; conversion time max {x.rs_p_t.max():.0f} s")
    i = (W.rs_p / W.ref).idxmax()
    say(f"largest converted / reference: {(W.rs_p / W.ref).max():.3f} on {W.loc[i, 'inp']} n = {W.loc[i, 'nn']}")
    say(f"CK exact never smaller than the converted WSPD: {bool((W.cke >= W.rs_p).all())}, also after merge: {bool((W.cke_m >= W.rs_p).all())}")
    o5 = (W.opt_p / W.opt).dropna()
    say(f"converted OPT, n = 500: mean {o5.mean():.4f} OPT, max {o5.max():.4f}")
    q = W[W.nn >= 1000].groupby("nn")[["qte", "cke", "qte_m", "cke_m", "ref"]].apply(
        lambda x: pd.Series({"qt": (x.qte / x.ref).mean(), "ck": (x.cke / x.ref).mean(),
                             "qtm": (x.qte_m / x.ref).mean(), "ckm": (x.cke_m / x.ref).mean()}))
    say(f"quadtree (exact) larger than CK exact on average for every n >= 1000: without merge {bool((q.qt > q.ck).all())}, "
        f"with merge {bool((q.qtm > q.ckm).all())}")
    for nn in (500, 1000, 2000, 10000):
        x = W[W.nn == nn]
        say(f"n = {nn}: sweep {(x.sw / x.ref).mean():.3f}, CK exact + merge {(x.cke_m / x.ref).mean():.3f} (x reference); "
            f"CK exact + merge <= sweep on {int((x.cke_m <= x.sw).sum())} of {len(x)}")
    for n in (500, 10000):
        say(f"n = {n}: sweep {sw[sw.nn == n].r.mean():.3f}, CK exact {cke[cke.nn == n].r.mean():.3f}, LB/OPT min "
            f"{lbopt[lbopt.nn == n].q.min() if (lbopt.nn == n).any() else float('nan'):.3f}, sweep max {sw[sw.nn == n].r.max():.3f}, "
            f"sweep + reshape max {rs[rs.nn == n].r.max():.4f} (x reference)")
    x = sw.dropna(subset=["opt"])
    say(f"open problems: wherever OPT is known, LB >= {lbopt.q.min():.3f} OPT and the sweep <= {x.ropt.max():.3f} OPT")
    # lower end of the approximation factor of the greedy sweep (without the redundancy pass, sweep_unpruned.csv):
    # its size divided by the smallest WSPD known for the instance (proven OPT or the smallest valid cover WSPD of any
    # run), which is at most its size divided by OPT
    u = read("sweep_unpruned.csv")
    if len(u):
        c = with_optx(read("exp2.csv", "exp2_opt.csv", "exp3.csv", "exp4.csv", "ls_all_rs.csv", "ls_cmp_rs.csv",
                           "ls_cmp_big.csv", "ls_all_big_5000.csv", "partition_exp.csv", "partition_exp_big.csv",
                           "partition_sweep.csv", "partition_ls.csv", "partition_ls_big.csv", "partition_ckls.csv",
                           "partition_ckls_big.csv"))
        c = c[(c.valid == "ok") & (c.eps == 0.05) & ~c.algo.str.endswith("+part")]
        c = c[(c.algo != "opt") | c["info"].astype(str).str.contains("optimal=1")]
        best = c.groupby(KEY)["size"].min()
        u = u[u.valid == "ok"].drop_duplicates(KEY, keep="last").set_index(KEY)
        q = u["size"] / best.reindex(u.index)
        i = q.idxmax()
        say(f"open problem (factor of the greedy sweep): greedy sweep without the redundancy pass / smallest known WSPD, "
            f"max {q.max():.4f} on {i[0]} (n = {i[1]}, seed {i[2]}; {q.notna().sum()} instances, "
            f"n = {u.index.get_level_values('n').min()}..{u.index.get_level_values('n').max()})")
    open(os.path.join(RESULTS, "text_numbers.txt"), "w").write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()

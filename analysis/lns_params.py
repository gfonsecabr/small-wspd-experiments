#!/usr/bin/env python3
"""Summary of the tests of the LNS parameters (experiments/lns_params.py; the paper's results: data/lns_params.csv).
  python3 analysis/lns_params.py            ->  results/lns_params_summary.md

Per n and start: for every configuration, the mean over the instances of the sample of size / reference at the end of
the budget, its ratio to the configuration "default" (the LNS before the tests) and to the paper's reshape run with
the same budget (ls_cmp_rs.csv; for n = 50000 the reshape runs of ls_all_big_5000.csv have 5000 s, so this column is
left out), the number of instances on which it ends smaller than "default", and the counters of LNS (calls, share of
accepted calls, share of calls with more than UMAX uncovered pairs, share of the time in CPLEX). At the end, the
numbers quoted in the paper: the largest change of the mean size caused by the CPLEX limit (TL0.5, TL10) or by the
number of candidates per pair (cap50, cap2000)."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import DATA, KEY, RESULTS, read, refs


def main():
    L = read("lns_params.csv")
    L = L[L.valid == "ok"].drop_duplicates(["round", "config"] + KEY + ["algo"], keep="last").copy()
    L["start"] = L.algo.str.extract(r"^ls-(ck|sweep)-")[0]
    for k in ("lnsCalls", "lnsOk", "lnsBig", "lnsFail"):
        L[k] = L["info"].str.extract(k + r"=(\d+)")[0].astype(float)
    L["start0"] = L["info"].str.extract(r"start=(\d+)")[0].astype(float)
    L["cplex"] = L["info"].str.extract(r"lnsCplex=([0-9.]+)")[0].astype(float) / L["info"].str.extract(r"lsTime=([0-9.]+)")[0].astype(float)
    _, ref = refs()
    L = L.join(ref.rename("ref"), on=KEY)
    L["r"] = L["size"] / L["ref"]
    P = read("ls_cmp_rs.csv")
    P = P[(P.valid == "ok") & P.algo.str.endswith("-reshape")].drop_duplicates(KEY + ["algo"], keep="last")
    P["start"] = P.algo.str.extract(r"^ls-(ck|sweep)-")[0]
    L = L.join(P.set_index(KEY + ["start"])["size"].rename("reshape"), on=KEY + ["start"])
    L = L.join(L[L.config == "default"].set_index(KEY + ["start"])["size"].rename("def"), on=KEY + ["start"])
    out = ["# Tests of the LNS parameters, eps = 1/20", "",
           "Means over the instances of the sample. `size/ref`: at the end of the budget; `start/ref`: the start; `/default`:",
           "size / size of the configuration default; `/reshape`: size / size of the paper's reshape run with the same budget;",
           "a name ending in `_x10` has 10 times the budget; `< default`: instances on which it ends smaller than default;",
           "`ok`: share of the accepted calls; `big`: share of the calls with more than UMAX uncovered pairs (no new",
           "candidates); `cplex`: share of the time in CPLEX.", ""]
    for (n, start), g in L.groupby(["n", "start"]):
        out += [f"## n = {n}, from {'CK exact' if start == 'ck' else 'the sweep'} ({g.drop_duplicates(KEY).shape[0]} instances)", ""]
        rows = []
        for cfg, x in g.groupby("config", sort=False):
            rows.append({"config": cfg, "K": x.K.iloc[0], "TL": x.tl.iloc[0], "cap": x["cap"].iloc[0],
                         "radius": x.radius.iloc[0], "UMAX": "inf" if x.umax.iloc[0] >= 10 ** 9 else x.umax.iloc[0], "N": len(x),
                         "start/ref": (x.start0 / x.ref).mean(), "size/ref": x.r.mean(), "/default": (x["size"] / x["def"]).mean(),
                         "/reshape": (x["size"] / x.reshape).mean(), "< default": int((x["size"] < x["def"]).sum()),
                         "calls": x.lnsCalls.mean(), "ok": (x.lnsOk / x.lnsCalls).mean(), "big": (x.lnsBig / x.lnsCalls).mean(),
                         "cplex": x.cplex.mean(), "fail": int(x.lnsFail.sum())})
        t = pd.DataFrame(rows).sort_values("size/ref")
        out += [t.to_markdown(index=False, floatfmt=("", "", "", "", "", "", "", ".3f", ".4f", ".4f", ".4f", "", ".0f", ".3f", ".3f", ".2f", "")), ""]
    # the CPLEX limit and the candidates per pair (round 1): largest change of the mean size against default
    r1 = L[L["round"] == 1]
    m = r1.groupby(["n", "start", "config"])["size"].mean().unstack("config")
    ch = {c: (m[c] / m["default"] - 1).abs().max() for c in ("TL0.5", "TL10", "cap50", "cap2000") if c in m}
    out += ["## Numbers of the paper", "",
            "Largest change of the mean size against default, over n and start (round 1): " +
            ", ".join(f"{c} {100 * v:.2f}%" for c, v in ch.items()) +
            (f"; maximum {100 * max(ch.values()):.2f}% (paper: at most 0.3%)." if ch else "."), ""]
    open(os.path.join(RESULTS, "lns_params_summary.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()

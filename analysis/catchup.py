"""Reshape from CK exact against the sweep and against reshape from the sweep, per instance -> DATA/catchup.csv.
Run from the repository root:  python3 analysis/catchup.py   (DATA=runs ... for your own runs)

It needs the trace of every reshape run from both starts (RA_{ck,sweep}_reshape_* for n <= 10000, RL_* above), which
a full rerun writes to runs/traces. The release ships only the result, data/catchup.csv, computed from the paper's runs.
It also needs the setup times of the local search (ls_setup.csv, group ls_setup of the experiments), which the traces
do not include.
Columns: sweep_size and sweep_s (size and running time of the sweep), ck_s (time of CK exact inside the local-search
run, tBuild), setup_s (setup of the local search from CK exact), reshape_s (time of the local search until its size
is below sweep_size; empty if it never is), total_s = ck_s + setup_s + reshape_s, to be compared with sweep_s, and
sweep_ahead: the largest value of (size from CK exact - size from the sweep) / size from the sweep over all total
times (time to build the start, plus setup, plus time of reshape) at which both runs have a WSPD and neither has ended.
It is positive if, at some total time, the run from the sweep is smaller. The setup from the sweep comes from
ls_setup.csv for n <= 10000; above, where ls_setup.csv has no run from the sweep, it is the lsSetup that the reshape
run itself reports (0 if the row has none, which can only increase sweep_ahead)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference import DATA, KEY, MIN_N, MAX_N, read, trace_files


def tag(r):
    return f"images_{r.gen.removesuffix('.instance')}_instance_0_0_1_20" if r.gen.endswith(".instance") else f"{r.gen}_{r.n}_{r.seed}_1_20"


def sweep_ahead(ck, tr, s, ts):
    """Largest (size from CK exact - size from the sweep) / size from the sweep at equal total time. ck, s: build times
    of the two starts; tr, ts: their traces. Both sizes are step functions of the total time, so it is enough to compare
    them at the times where one of them changes, from the end of the later build to the end of the first run."""
    lo, hi = max(ck, s), min(ck + tr.t.iloc[-1], s + ts.t.iloc[-1])
    x = np.unique(np.concatenate([ck + tr.t.values, s + ts.t.values, [lo]]))
    x = x[(x >= lo) & (x <= hi)]
    at = lambda t, b: t["size"].values[np.maximum(np.searchsorted(t.t.values, x - b, side="right") - 1, 0)]
    a, b = at(tr, ck), at(ts, s)
    return ((a - b) / b).max()


def main():
    sw = read("exp2.csv", "exp3.csv", "exp4.csv")
    sw = sw[(sw.algo == "lazyprune") & (sw.valid == "ok")].drop_duplicates(KEY, keep="last").set_index(KEY)
    ls = read("ls_all_rs.csv", "ls_all_big_5000.csv")
    ls = ls[(ls.valid == "ok") & (ls.n >= MIN_N) & (ls.n <= MAX_N)]
    lsw = ls[ls.algo == "ls-sweep-reshape"].drop_duplicates(KEY, keep="last").set_index(KEY)
    su = read("ls_setup.csv")
    if not len(su): sys.exit(f"{DATA}/ls_setup.csv is missing; nothing written")
    su = su.assign(setup=su["info"].str.extract(r"lsSetup=([0-9.]+)")[0].astype(float)).drop_duplicates(KEY + ["algo"], keep="last")
    su = su.set_index(KEY + ["algo"]).setup
    ls = ls[ls.algo == "ls-ck-reshape"].drop_duplicates(KEY, keep="last")
    rows, missing = [], 0
    for r in ls.sort_values(["n", "gen", "seed"]).itertuples():
        pre = "RL" if r.n >= 20000 else "RA"
        f, fs = trace_files(f"{pre}_ck_reshape_{tag(r)}"), trace_files(f"{pre}_sweep_reshape_{tag(r)}")
        if len(f) != 1 or len(fs) != 1:
            missing += 1
            continue
        s = sw.loc[(r.gen, r.n, r.seed, r.eps)]
        tr = pd.read_csv(f[0], header=None, names=["t", "size"])
        assert tr["size"].iloc[-1] == r.size, f[0]
        h = tr[tr["size"] < s["size"]]
        ck = float(r.info.split("tBuild=")[1].split(";")[0])
        t = h.t.iloc[0] if len(h) else float("nan")
        key = (r.gen, r.n, r.seed, r.eps)
        rs, cs, ss = lsw.loc[key], su.get(key + ("ls-ck-reshape",)), su.get(key + ("ls-sweep-reshape",))
        if ss is None:   # no separate measurement (n >= 20000): the setup time recorded by the run itself
            ss = float(rs["info"].split("lsSetup=")[1].split(";")[0]) if "lsSetup=" in rs["info"] else 0.0
        if cs is None:
            missing += 1
            continue
        ts = pd.read_csv(fs[0], header=None, names=["t", "size"])
        assert ts["size"].iloc[-1] == rs["size"], fs[0]
        rows.append(dict(gen=r.gen, n=r.n, seed=r.seed, eps=r.eps, sweep_size=int(s["size"]), sweep_s=s.time, ck_s=ck,
                         setup_s=cs, reshape_s=t, total_s=ck + cs + t,
                         sweep_ahead=sweep_ahead(ck + cs, tr, float(rs["info"].split("tBuild=")[1].split(";")[0]) + ss, ts)))
    if missing:
        sys.exit(f"the reshape traces or the setup time of {missing} of {len(ls)} instances are missing in {DATA}; nothing written")
    out = os.path.join(DATA, "catchup.csv")
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.6g")
    print(f"{out}: {len(rows)} instances")


if __name__ == "__main__":
    main()

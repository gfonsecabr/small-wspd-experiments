#!/usr/bin/env python3
"""Compares the outputs of your runs (default runs/) with the paper's results (data/), row by row.
  python3 experiments/check.py [RUNS]

Every row of RUNS/<group>.csv is matched with the row of data/<group>.csv of the same instance and algorithm (and
conversion). Outputs that do not depend on time are compared exactly:
  the constructions (sizes, and the number of steps of the sweep, with or without the redundancy pass, which is the
  lower bound), CK, the quadtrees, the
  conversions to partition WSPDs, OPT and optx where both runs proved optimality, and reshape when replayed
  (run.py --replay; partition_ls and partition_ckls always replay): size, rounds, moves, and deletions.
The other local-search runs stop after a time budget, so their sizes depend on the speed of the machine; for them the
script reports the mean of (your size) / (paper's size). Every output must also be valid (column valid = ok; exp
verifies coverage and exact separation of every WSPD, and that the partition WSPDs cover every pair exactly once).
Exit status 1 if a time-independent output differs or an output is invalid.
"""
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
KEY = ["gen", "n", "seed", "eps", "algo", "conv"]
COUNTERS = r"(rounds=\d+;reshapeMoves=\d+;reshapeKills=\d+)"


def load(path):
    d = pd.read_csv(path)
    d["algo"] = d.algo.str.replace(r"^load-", "", regex=True)   # WSPDs loaded from dumps
    d["conv"] = d["info"].astype(str).str.extract(r"conv=([^;]+)")[0].fillna("")
    return d.drop_duplicates(KEY, keep="last").set_index(KEY)


def main():
    runs = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "runs"))
    groups = sorted(f[:-4] for f in os.listdir(runs) if f.endswith(".csv") and f != "memory.csv")
    bad = 0
    for g in groups:
        r = load(os.path.join(runs, f"{g}.csv"))
        inv = r[r.valid != "ok"]
        if len(inv):
            bad += len(inv)
            print(f"{g}: {len(inv)} INVALID outputs: {list(inv.index[:3])}")
        if not os.path.exists(os.path.join(DATA, f"{g}.csv")):
            print(f"{g}: no paper data")
            continue
        p = load(os.path.join(DATA, f"{g}.csv"))
        j = r.join(p, how="left", rsuffix="_p")
        for algo, x in j.groupby(level="algo"):
            miss = int(x["size_p"].isna().sum())
            x = x[x["size_p"].notna()]
            if not len(x):
                print(f"{g:18s} {algo:24s} {miss} rows without a paper row")
                continue
            note = f", {miss} without a paper row" if miss else ""
            base = algo.removesuffix("+part")
            if base in ("opt", "optx") and not g.startswith("partition"):
                x = x[x["info"].str.contains("optimal=1") & x["info_p"].str.contains("optimal=1")]
            replayed = base.endswith("-reshape") and x["info"].str.contains("rounds=").all() and \
                (x["info"].str.extract(COUNTERS)[0] == x["info_p"].str.extract(COUNTERS)[0]).any()
            if base.startswith("ls-") and not replayed and not algo.endswith("+part"):
                print(f"{g:18s} {algo:24s} {len(x):4d} rows, time budget: mean size / paper's size {(x['size'] / x['size_p']).mean():.4f}{note}")
                continue
            same = x["size"] == x["size_p"]
            what = "sizes"
            if base in ("lazyprune", "lazy") and not algo.endswith("+part"):
                same &= x["info"].str.extract(r"(steps=\d+)")[0] == x["info_p"].str.extract(r"(steps=\d+)")[0]
                what = "sizes and steps"
            elif replayed and not algo.endswith("+part"):
                same &= x["info"].str.extract(COUNTERS)[0] == x["info_p"].str.extract(COUNTERS)[0]
                what = "sizes and counters (replayed)"
            diff = x[~same]
            bad += len(diff)
            print(f"{g:18s} {algo:24s} {len(x):4d} rows, {int(same.sum())} identical {what}{note}"
                  + (f"  DIFFERENT: {list(diff.index[:3])}" if len(diff) else ""))
    print("all time-independent outputs identical and valid" if not bad else f"{bad} differences or invalid outputs")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""The tests of the LNS parameters (appendix of the paper, paragraph "Local search": "The parameters of LNS are the best
ones among those we tested ..."). Results of the paper's tests: data/lns_params.csv; summary: analysis/lns_params.py.
  python3 experiments/lns_params.py [ROUND] [-P 11] [--ns 1000,10000] [--mem-gb 100] [--out runs]

Each configuration below changes some parameters of LNS against the configuration "default" (the LNS before the tests:
K = 16 neighbors within radius 3 eps * scale, candidates only if at most UMAX = 400 pairs are uncovered, CAP = 300
candidates per pair, CPLEX limit TL = 2 s). The runs are LNS from CK exact and from the sweep on the local-search
samples with the paper's budgets: n = 1000 (6 generators with seed 1 and 4 CG:SHOP instances, 100 s), n = 10000 and
n = 50000 (uniform 1, cluster 1, euro-night, 1000 s), and in round 3 n = 3000 (the same three, 1000 s). The paper's
parameters (LSConfig::lnsDefaults in src/ls.hpp) are those of fast_r2_U1000 for n < 5000 and of fast_r2_K32 (with
UMAX = 400) for n >= 5000.

The configurations whose names begin with "fast" use the cached neighbor search of src/ls.hpp (same moves, 3 to 9
times faster). The others (round 1 and r2_U2000) ran in the paper's tests with the earlier neighbor search, which
recomputed the scale and the centroids of every scanned dumbbell; src/ls.hpp has only the cached one, so these
configurations make more moves within their budgets when rerun and end smaller than in data/lns_params.csv.

Appends one row per run to OUT/lns_params.csv (round, config, K, tl, cap, radius, umax, then the row of exp; the info
column has the counters of LNS) and writes the traces OUT/traces/LP_<config>_<start>_<instance>.csv. It resumes: runs
already in the output are skipped. The paper's tests ran round 1 with --ns 1000,10000, round 2 with
--ns 1000,10000,50000, and round 3 with --ns 1000,3000,10000,50000 (845 runs). Round 3 lost 21 runs when the machine
ran out of memory (UMAX = 4000 and no UMAX for n = 10000, and both configurations for n = 50000); they were not rerun. Memory: about 0.3, 0.6, 2.5, and 25 GB per run for n = 1000, 3000, 10000, 50000,
but much more with a large or no UMAX for n = 10000 (31 GB for Uinf and U4000)."""
import argparse
import os
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jobs as J

HEADER = "round,config,K,tl,cap,radius,umax,gen,n,seed,eps,algo,size,valid,time,info"
BUDGET = {1000: 100, 3000: 1000, 10000: 1000, 50000: 1000}
MEM = {1000: 0.3, 3000: 0.6, 10000: 2.5, 50000: 25}
DEFAULT = dict(K=16, tl=2, cap=300, radius=3, umax=400, tmult=1, ns=None)
INF = 10 ** 9
ROUNDS = {
    # round 1: one parameter at a time around the default
    "1": {"default": {}, "K4": dict(K=4), "K8": dict(K=8), "K32": dict(K=32), "K64": dict(K=64),
          "U100": dict(umax=100), "U2000": dict(umax=2000), "Uinf": dict(umax=INF),
          "cap50": dict(cap=50), "cap2000": dict(cap=2000), "TL0.5": dict(tl=0.5), "TL10": dict(tl=10),
          "r2": dict(radius=2), "r5": dict(radius=5)},
    # round 2: combinations of the best values of round 1, the cached neighbor search, 10 times the budget for n = 1000
    # (_x10), and n = 50000
    "2": {"fast": {}, "fast_r2": dict(radius=2), "fast_r2_U2000": dict(radius=2, umax=2000),
          "fast_r2_U1000": dict(radius=2, umax=1000), "fast_r2_K32": dict(radius=2, K=32),
          "fast_r1.5": dict(radius=1.5), "r2_U2000": dict(radius=2, umax=2000, ns=[1000, 10000]),
          "fast_r2_x10": dict(radius=2, tmult=10, ns=[1000]),
          "fast_r2_U2000_x10": dict(radius=2, umax=2000, tmult=10, ns=[1000])},
}
# round 3: UMAX on a grid, with radius 2, for n = 1000, 3000, 10000 (and two values for n = 50000)
_cap = lambda u, ns: dict(radius=2, umax=u, ns=ns)
ROUNDS["3"] = {**{f"fast_r2_U{'inf' if u == INF else u}": _cap(u, [1000, 3000, 10000] + ([50000] if u == 100 else []))
                  for u in (50, 100, 200, 800, 4000, INF)},
               "fast_r2": _cap(400, [3000]), "fast_r2_U1000": _cap(1000, [3000]), "fast_r2_U2000": _cap(2000, [3000, 50000]),
               "fast_r2_U1000_x10": dict(radius=2, umax=1000, tmult=10, ns=[1000]),
               "fast_r2_Uinf_x10": dict(radius=2, umax=INF, tmult=10, ns=[1000])}
BIG_CONFIGS = {"2": ["fast", "fast_r2", "fast_r2_K32"], "3": ["fast_r2_U100", "fast_r2_U2000"]}   # those run with n = 50000


def instances(n):
    s1000, s10000, s50000 = J.samples()
    if n == 1000: return s1000
    return [J.Inst("uniform", n, 1), J.Inst("cluster", n, 1), J.Inst("euro-night", n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("round", nargs="?", default="1", choices=sorted(ROUNDS))
    ap.add_argument("-P", type=int, default=11)
    ap.add_argument("--ns", default="1000,10000", help="comma-separated n (1000, 3000, 10000, 50000)")
    ap.add_argument("--mem-gb", type=float, default=100)
    ap.add_argument("--out", default=os.path.join(J.ROOT, "runs"))
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(os.path.join(out, "traces"), exist_ok=True)
    OUT, LOG = os.path.join(out, "lns_params.csv"), os.path.join(out, "lns_params.log")
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT).read().splitlines()[1:]:
            c = l.split(",")
            done.add((c[1], c[7], c[11].split("-")[1], c[8]))   # config, gen, start, n (as printed)
    else:
        open(OUT, "w").write(HEADER + "\n")
    todo = []
    for n in (int(x) for x in a.ns.split(",")):
        for name, ch in ROUNDS[a.round].items():
            c = {**DEFAULT, **ch}
            if (c["ns"] is not None and n not in c["ns"]) or (n == 50000 and name not in BIG_CONFIGS.get(a.round, [])): continue
            for i in instances(n):
                for start in ("ck", "sweep"):
                    if (name, i.gen, start, str(n)) not in done:   # n printed by exp = nominal n on these instances
                        todo.append(dict(n=n, name=name, c=c, inst=i, start=start))
    todo.sort(key=lambda j: (j["n"], j["name"] != "default"))
    lock, cond = threading.Lock(), threading.Condition()
    state = dict(mem=0.0, run=0, n=0, bad=0)

    def log(s):
        with lock:
            open(LOG, "a").write(time.strftime("%F %T ") + s + "\n")

    log(f"START round {a.round}: {len(todo)} runs, P = {a.P}")
    print(f"{len(todo)} runs to do", flush=True)

    def work(j):
        c, i = j["c"], j["inst"]
        env = dict(os.environ, WSPD_LS_START=j["start"], WSPD_LS_MOVES="lns", WSPD_LS_TIME=str(BUDGET[j["n"]] * c["tmult"]),
                   WSPD_LS_K=str(c["K"]), WSPD_LS_LNS_TL=str(c["tl"]), WSPD_LS_LNS_CAP=str(c["cap"]),
                   WSPD_LS_LNS_RADIUS=str(c["radius"]), WSPD_LS_LNS_UMAX=str(c["umax"]),
                   WSPD_LS_TRACE=os.path.join(out, "traces", f"LP_{j['name']}_{j['start']}_{i.tag}.csv"))
        t0 = time.time()
        r = subprocess.run([J.EXP] + i.args + ["1", "20", "ls"], env=env, capture_output=True, text=True)
        row = r.stdout.strip()
        ok = r.returncode == 0 and row.count(",") >= 8 and row.split(",")[6] == "ok"
        with cond:
            if ok:
                open(OUT, "a").write(f"{a.round},{j['name']},{c['K']},{c['tl']},{c['cap']},{c['radius']},{c['umax']},{row}\n")
            else:
                state["bad"] += 1
            state["n"] += 1; state["run"] -= 1; state["mem"] -= MEM[j["n"]]
            cond.notify_all()
        log(f"{'OK' if ok else 'FAIL'} {j['name']}:{j['start']}:{i.tag} wall={time.time() - t0:.0f}s"
            + ("" if ok else f" rc={r.returncode} {row[:200]} {r.stderr.strip()[:200]}"))

    with cond:
        while todo or state["run"]:
            while todo and state["run"] < a.P and (state["mem"] + MEM[todo[0]["n"]] <= a.mem_gb or not state["run"]):
                j = todo.pop(0)
                state["run"] += 1; state["mem"] += MEM[j["n"]]
                threading.Thread(target=work, args=(j,)).start()
            cond.wait()
    log(f"END round {a.round}: {state['n']} runs, {state['bad']} failed")
    sys.exit(1 if state["bad"] else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""The largest sets of pairwise incompatible pairs (alpha) of the paper's instances with n = 200..500 (Fig. lb).
  python3 experiments/alpha.py -P 6

For every instance, src/exp ... alpha computes alpha by integer programming over all maximal dumbbells (CPLEX, MIP gap
0, limit WSPD_OPT_TL seconds, default 1800) and writes the chosen pairs with WSPD_ALPHA; src/check_lb checks that they
are pairwise incompatible. One row per instance is appended to OUT/alpha.csv (default OUT = runs/):
  gen,n,seed,eps   the instance (n nominal)
  alpha            the number of chosen pairs, a lower bound on OPT
  optimal          1 if CPLEX proved that no larger set exists
  ub               the bound of CPLEX
  cols             the number of maximal dumbbells
  check            ok if check_lb accepts the pairs
  time             seconds of exp
Instances that OUT/alpha.csv already lists are skipped. The paper's results are in data/alpha.csv.
Options: --min-n, --max-n (default 200, 500), --only REGEX, -P, --out."""
import argparse, csv, os, re, subprocess, sys, tempfile, threading
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jobs as J
import run as R

CHECK = os.environ.get("WSPD_CHECK_LB") or os.path.join(J.ROOT, "src", "check_lb")
HEADER = ["gen", "n", "seed", "eps", "alpha", "optimal", "ub", "cols", "check", "time"]


def one(inst):
    with tempfile.TemporaryDirectory() as tmp:
        pre = os.path.join(tmp, "a")
        r = subprocess.run([J.EXP] + inst.args + ["1", "20", "alpha"], env={**os.environ, "WSPD_ALPHA": pre}, capture_output=True, text=True)
        v = r.stdout.strip().split(",")
        if r.returncode or len(v) < 9 or v[4] != "alpha":
            return None, f"alpha failed on {inst.name}: {r.stdout.strip()} {r.stderr.strip()}"
        c = subprocess.run([CHECK, pre + ".pts", pre + ".alpha", "1", "20"], capture_output=True, text=True)
    m = re.match(r"OK lower bound (\d+):", c.stdout)
    check = "ok" if c.returncode == 0 and m and int(m.group(1)) == int(v[5]) else (c.stdout.strip() or "failed").replace(",", ";")
    info = dict(kv.split("=") for kv in v[8].split(";"))
    return [inst.gen, inst.n, inst.seed, "0.05", v[5], info["optimal"], info["ub"], info["cols"], check, v[7]], None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-n", type=int, default=200)
    ap.add_argument("--max-n", type=int, default=500)
    ap.add_argument("--only", default="")
    ap.add_argument("-P", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(J.ROOT, "runs"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, "alpha.csv")
    done = set()
    if os.path.exists(p):
        done = {(r["gen"], int(r["n"]), int(r["seed"])) for r in csv.DictReader(open(p))}
    else:
        with open(p, "w", newline="") as f: csv.writer(f, lineterminator="\n").writerow(HEADER)
    insts = [i for n in J.SIZES_MID if a.min_n <= n <= a.max_n for i in J.insts(n) if re.search(a.only, i.name) and (i.gen, i.n, i.seed) not in done]
    print(f"{len(insts)} instances, output {p}", flush=True)
    lock, bad = threading.Lock(), [0]

    def work(inst):
        row, err = one(inst)
        with lock:
            if err: bad[0] += 1; print(err, flush=True); return
            bad[0] += row[8] != "ok" or row[5] != "1"
            with open(p, "a", newline="") as f: csv.writer(f, lineterminator="\n").writerow(row)
            print(",".join(map(str, row)), flush=True)

    with ThreadPoolExecutor(a.P) as ex: list(ex.map(work, insts))
    print("all checked and optimal" if not bad[0] else f"{bad[0]} failures or unproven")
    sys.exit(1 if bad[0] else 0)


if __name__ == "__main__":
    main()

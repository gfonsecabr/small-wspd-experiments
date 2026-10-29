#!/usr/bin/env python3
"""Computes and checks the lower bounds of the paper independently of the code that computes the WSPDs.
  python3 experiments/certify.py --max-n 1000 -P 6

For every selected instance of the paper, the script runs the sweep (src/exp ... lazyprune) with WSPD_ANCHORS, which
writes the points and the anchors of the steps, then src/lb_swap, which enlarges the set of anchors by local search
(LB1 after the 1-for-2 swaps, LB after the plateau passes; the lower bound of the paper is LB), and src/check_lb, a
separate program that shares no code with the others, on each of the three sets. check_lb verifies with exact integer
tests that the pairs of a set are distinct pairs of points and that no eps-separated dumbbell covers two of them, so
every eps-WSPD of the points has at least as many dumbbells as the set has pairs. The script also compares the three
numbers with those of the paper's runs (steps= in REF/exp2.csv, exp3.csv, exp4.csv, and REF/lb_check.csv).

It appends one row per instance to OUT/lb_check.csv (default OUT = runs/):
  gen,n,seed,eps   the instance, as in the other CSV files
  steps            the number of anchors of the sweep (LB0 in the paper)
  lb1              the size of the set after the free pairs and the 1-for-2 swaps (LB1)
  lb               the size after the plateau passes (LB, the lower bound of the paper)
  check            ok if check_lb accepts the three sets with these sizes, else its message
  paper            same if the numbers equal those of the paper's runs (LB only with --plateau 30), DIFFERENT if not,
                   none if there is no such row
  tSweep           seconds of the sweep
  tSwap0, tSwap    seconds of lb_swap until the end of the passes without plateau moves, and in all
  passes0          number of passes without plateau moves
  check_s          seconds taken by check_lb on the three sets
Instances that OUT/lb_check.csv already lists are skipped. The paper's instances are in data/lb_check.csv. Exit
status 1 if a check fails or a number differs from the paper's.

Options: --min-n/--max-n (nominal n), --only REGEX (on the instance name, e.g. uniform_1000_1 or euro-night-0001000),
-P (parallel jobs), --plateau (passes with plateau moves, default 30 as in the paper), --ref DIR (compare with the
files of DIR instead of data/, e.g. the output folder of run.py), --mem-gb (memory budget, default 16; n = 50000 needs
up to 15 GB; jobs start by decreasing memory estimate while the estimates fit; the estimate of an instance is the
larger of the peak of its sweep in data/memory.csv and the 5.2 n^2 bytes of lb_swap, with a margin), --out, and --keep
(keep the point and pair files in OUT/anchors/ instead of deleting them after the check).
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jobs as J
import run as R

CHECK = os.environ.get("WSPD_CHECK_LB") or os.path.join(J.ROOT, "src", "check_lb")
LBSWAP = os.environ.get("WSPD_LB_SWAP") or os.path.join(J.ROOT, "src", "lb_swap")
HEADER = "gen,n,seed,eps,steps,lb1,lb,check,paper,tSweep,tSwap0,tSwap,passes0,check_s\n"


def paper_numbers(ref):
    """(gen, seed, nominal n) -> (steps, lb1, lb) of the paper's runs in the folder ref: steps from the sweep rows of
    exp2/3/4.csv, lb1 and lb from lb_check.csv (None where missing)."""
    s = {}
    for g in ("exp2", "exp3", "exp4"):
        p = os.path.join(ref, f"{g}.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                if r["algo"] == "lazyprune":
                    s[r["gen"], int(r["seed"]), R.nominal_n(r["gen"], int(r["n"]))] = [int(re.search(r"steps=(\d+)", r["info"]).group(1)), None, None]
    p = os.path.join(ref, "lb_check.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            k = r["gen"], int(r["seed"]), R.nominal_n(r["gen"], int(r["n"]))
            if k in s and "lb1" in r: s[k][1:] = [int(r["lb1"]), int(r["lb"])]
    return {k: tuple(v) for k, v in s.items()}


def paper_peaks():
    """(gen, seed, nominal n) -> peak memory (GB) of the sweep of the paper's run on this instance."""
    m = {}
    p = os.path.join(J.ROOT, "data", "memory.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            if r["algo"] == "lazyprune":
                k = r["gen"], int(r["seed"]), R.nominal_n(r["gen"], int(r["n"]))
                m[k] = max(m.get(k, 0.0), float(r["rss_gb"]))
    return m


def certify(inst, out, keep, plateau):
    """Returns the CSV row of the instance (without the column paper) and an error message or None."""
    pre = f"{out}/anchors/{inst.name}"
    t0 = time.time()
    r = subprocess.run([J.EXP] + inst.args + ["1", "20", "lazyprune"], env={**os.environ, "WSPD_ANCHORS": pre},
                       capture_output=True, text=True)
    t_sweep = time.time() - t0
    row = r.stdout.strip().split(",")
    if r.returncode or len(row) < 9 or row[6] != "ok":
        return None, f"sweep failed on {inst.name}: {r.stdout.strip()} {r.stderr.strip()}"
    steps = int(re.search(r"steps=(\d+)", row[8]).group(1))
    t0 = time.time()
    w = subprocess.run([LBSWAP, f"{pre}.pts", f"{pre}.anchors", "1", "20", f"{pre}.lb", str(plateau)], capture_output=True, text=True)
    t_swap = time.time() - t0
    v = w.stdout.strip().split(",")
    if w.returncode or len(v) < 10:
        return None, f"lb_swap failed on {inst.name}: {w.stdout.strip()} {w.stderr.strip()}"
    if plateau > 0:
        lb1, lb, passes0, t_swap0 = int(v[5]), int(v[10]), int(v[11]), float(v[14])
        sets = [(f"{pre}.anchors", steps), (f"{pre}.lb.swaps", lb1), (f"{pre}.lb", lb)]
    else:
        lb1 = lb = int(v[5]); passes0, t_swap0 = int(v[6]), float(v[9])
        sets = [(f"{pre}.anchors", steps), (f"{pre}.lb", lb)]
    check, t0 = "ok", time.time()
    if int(v[1]) != steps: check = "lb_swap read %s anchors, the sweep reports %d steps" % (v[1], steps)
    for f, size in sets:
        c = subprocess.run([CHECK, f"{pre}.pts", f, "1", "20"], capture_output=True, text=True)
        m = re.match(r"OK lower bound (\d+):", c.stdout)
        # the certified number must be the one that exp or lb_swap reports
        if check == "ok" and not (c.returncode == 0 and m and int(m.group(1)) == size):
            check = (c.stdout.strip() or c.stderr.strip() or "failed").replace(",", ";")
    dt = time.time() - t0
    if not keep:
        for f in (f"{pre}.pts", f"{pre}.anchors", f"{pre}.lb", f"{pre}.lb.swaps"):
            if os.path.exists(f): os.remove(f)
    return [row[0], row[1], row[2], row[3], str(steps), str(lb1), str(lb), check,
            f"{t_sweep:.2f}", f"{t_swap0:.2f}", f"{t_swap:.2f}", str(passes0), f"{dt:.2f}"], None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-n", type=int, default=0)
    ap.add_argument("--max-n", type=int, default=10 ** 9)
    ap.add_argument("--only", default="")
    ap.add_argument("-P", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--mem-gb", type=float, default=16)
    ap.add_argument("--out", default=os.path.join(J.ROOT, "runs"))
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--ref", default=os.path.join(J.ROOT, "data"))
    ap.add_argument("--plateau", type=int, default=30)
    a = ap.parse_args()
    for f in (J.EXP, LBSWAP, CHECK):
        if not os.access(f, os.X_OK):
            sys.exit(f"missing {f}: run make in src/")
    out = os.path.abspath(a.out)
    os.makedirs(f"{out}/anchors", exist_ok=True)
    csvp = f"{out}/lb_check.csv"
    done = set()
    if os.path.exists(csvp):
        done = {(r["gen"], int(r["seed"]), R.nominal_n(r["gen"], int(r["n"]))) for r in csv.DictReader(open(csvp))}
    else:
        open(csvp, "w").write(HEADER)
    insts = [i for n in J.SIZES_MID + J.SIZES_LARGE + J.SIZES_BIG if a.min_n <= n <= a.max_n for i in J.insts(n)
             if re.search(a.only, i.name) and (*i.key(), i.n) not in done]
    peaks, own, paper = R.peak_table(), paper_peaks(), paper_numbers(a.ref)
    # memory estimate: the larger of the peak of the sweep of the paper's run on the same instance (the anchors add less
    # than 1 GB; without such a run, the estimate of run.py) and the 5.2 n^2 bytes of lb_swap, with a margin
    mem = {i.name: max(1.1 * own[(*i.key(), i.n)] + 1.5 if (*i.key(), i.n) in own else R.estimate(J.job("exp", i, "lazyprune"), peaks),
                       1.1 * 5.2e-9 * i.n ** 2 + 1.5) for i in insts}
    insts.sort(key=lambda i: -mem[i.name])
    print(f"{len(insts)} instances to check, P = {a.P}, memory budget {a.mem_gb:g} GB, output {csvp}")
    cond = threading.Condition()
    state = dict(used=0.0, running=0, bad=0, n=0)

    def work(inst):
        row, err = certify(inst, out, a.keep, a.plateau)
        with cond:
            state["n"] += 1
            if err:
                state["bad"] += 1
                print(err, flush=True)
            else:
                ref = paper.get((*inst.key(), inst.n))
                mine = (int(row[4]), int(row[5]), int(row[6]))
                # compared: steps, lb1 if the reference has it, and lb if it has it and --plateau is the paper's 30
                pairs = [] if ref is None else [(ref[0], mine[0])] + ([(ref[1], mine[1])] if ref[1] is not None else []) + \
                        ([(ref[2], mine[2])] if ref[2] is not None and a.plateau == 30 else [])
                same = "none" if ref is None else "same" if all(x == y for x, y in pairs) else "DIFFERENT"
                row.insert(8, same)
                state["bad"] += row[7] != "ok" or same == "DIFFERENT"
                with open(csvp, "a") as f:
                    f.write(",".join(row) + "\n")
                print(f"{time.strftime('%T')} {state['n']}/{len(insts)} {inst.name}: steps {row[4]}, LB1 {row[5]}, LB {row[6]}, check {row[7]}, paper {same}", flush=True)
            state["used"] -= mem[inst.name]
            state["running"] -= 1
            cond.notify_all()

    pending = list(insts)
    with cond:
        while pending or state["running"]:
            # largest estimate first, as long as the estimates fit in the budget (a job above the budget runs alone)
            for i in list(pending):
                if state["running"] < a.P and (state["used"] + mem[i.name] <= a.mem_gb or not state["running"]):
                    pending.remove(i)
                    state["used"] += mem[i.name]
                    state["running"] += 1
                    threading.Thread(target=work, args=(i,)).start()
            cond.wait()
    print("every lower bound is certified and equal to the paper's" if not state["bad"] else f"{state['bad']} failures or differences")
    sys.exit(1 if state["bad"] else 0)


if __name__ == "__main__":
    main()

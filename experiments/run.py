#!/usr/bin/env python3
"""Runs the jobs of jobs.py in parallel, with a memory budget, and writes their outputs to OUT (default runs/; the
paper's results in data/ are never touched):
  OUT/<group>.csv   the CSV lines printed by exp/optx (header gen,n,seed,eps,algo,size,valid,time,info)
  OUT/memory.csv    peak resident memory and wall time of every job
  OUT/jobs.log      one line per finished job (OK, FAIL, or TIMEOUT); jobs logged OK are skipped when run again
  OUT/traces/       local-search traces (time since the start of the local search, size), OUT/dumps/ dumped WSPDs

Examples (from the repository root):
  python3 experiments/run.py --list --groups exp2                      # show the jobs
  python3 experiments/run.py --groups exp2,exp2_opt --max-n 300 -P 8   # constructions and OPT, n <= 300
  python3 experiments/run.py --groups ls_all_rs --max-n 1000 --replay  # local search replayed exactly (see below)
  python3 experiments/run.py                                           # everything (about 400 CPU hours)

Selection: --groups (comma-separated, default all), --min-n/--max-n (nominal n, the n of the generator or of the
CG:SHOP file name), --only REGEX (on the job id).
--replay: reshape jobs stop after the number of rounds of the same run in data/ (WSPD_LS_ROUNDS) instead of after
their time budget. Reshape is deterministic, so this reproduces the paper's local-search outputs exactly (size and
counters), independently of the speed of the machine; traces and times then differ from the paper's. LNS (which calls
CPLEX under time limits) always runs with its time budget.
partition_ls and partition_ckls always replay (from OUT/ls_all_rs.csv if it has the instance, else from
data/ls_all_rs.csv).
Scheduling: at most P jobs, largest memory estimate first, as long as the estimates fit in --mem-gb. The estimate of a
job is 1.25 x the largest peak memory measured in the paper's runs (data/memory.csv) for the same algorithm and n,
plus 0.5 GB; OPT is estimated as 12 (n/1000)^2 + 0.5 GB. A job whose estimate exceeds the budget runs alone.
Running times depend on the load: the paper ran at most 11 single-threaded jobs on 12 cores, without other load.
"""
import argparse
import csv
import os
import re
import resource
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jobs as J

HEADER = "gen,n,seed,eps,algo,size,valid,time,info\n"
GB = 1e9


def peak_table():
    """(algo, nominal n) -> largest peak memory (GB) in data/memory.csv."""
    m = {}
    p = os.path.join(J.ROOT, "data", "memory.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            k = (r["algo"], nominal_n(r["gen"], int(r["n"])))
            m[k] = max(m.get(k, 0.0), float(r["rss_gb"]))
    return m


def nominal_n(gen, n):
    k = re.search(r"-(\d{7})\.instance$", gen)
    return int(k.group(1)) if k else n


def estimate(j, peaks):
    n, a = j["inst"].n, j["algo"]
    if a == "opt":
        return 12 * (n / 1000) ** 2 + 0.5
    if a == "optx":
        return 4.0
    base = {"load": "ls-sweep-reshape", "qtexact": "ckexact", "qtexactrand": "ckexact", "lazy": "lazyprune"}.get(a, a)
    m = peaks.get((base, n))
    if m is None and base.endswith("-lns") and (base[:-3] + "reshape", n) in peaks:
        m = 1.2 * peaks[base[:-3] + "reshape", n]   # LNS without a measured run: the tables of reshape plus CPLEX
    if m is None:
        m = max([v for (b, nn), v in peaks.items() if b == base and nn <= n] or [1.0]) * 4
    if "WSPD_PARTITION" in j["env"]:
        m *= 2   # the conversion keeps the cover and the partition WSPD
    return 1.25 * m + 0.5


def rounds_of(csvpath, j, algo):
    """Number of rounds logged by the local-search row of the instance of j in csvpath (None if absent)."""
    if not os.path.exists(csvpath):
        return None
    gen, seed = j["inst"].key()
    for r in csv.DictReader(open(csvpath)):
        if r["gen"] == gen and int(r["seed"]) == seed and r["algo"] == algo and nominal_n(gen, int(r["n"])) == j["inst"].n:
            return int(re.search(r"rounds=(\d+)", r["info"]).group(1))
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--groups", default=",".join(J.GROUPS))
    ap.add_argument("--min-n", type=int, default=0)
    ap.add_argument("--max-n", type=int, default=10 ** 9)
    ap.add_argument("--only", default="")
    ap.add_argument("-P", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--mem-gb", type=float, default=16)
    ap.add_argument("--out", default=os.path.join(J.ROOT, "runs"))
    ap.add_argument("--replay", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    groups = a.groups.split(",")
    for g in groups:
        if g not in J.GROUPS:
            sys.exit(f"unknown group {g}; groups: {', '.join(J.GROUPS)}")
    sel = [j for j in J.make_jobs() if j["group"] in groups and a.min_n <= j["inst"].n <= a.max_n
           and re.search(a.only, j["id"])]
    peaks = peak_table()
    if a.list:
        for j in sel:
            env = " ".join(f"{k}={v}" for k, v in j["env"].items())
            print(f"{j['id']:70s} {estimate(j, peaks):6.1f} GB  {env} {' '.join(j['cmd'][1:])}")
        print(f"{len(sel)} jobs")
        return
    for f in (J.EXP, J.OPTX):
        if not os.access(f, os.X_OK):
            sys.exit(f"missing {f}: run make in src/")
    out = os.path.abspath(a.out)
    for d in (out, f"{out}/traces", f"{out}/dumps", f"{out}/tmp"):
        os.makedirs(d, exist_ok=True)
    logp = f"{out}/jobs.log"
    done = set()
    if os.path.exists(logp):
        done = {l.split()[2] for l in open(logp) if len(l.split()) > 2 and l.split()[1] == "OK"}
    pending = [j for j in sel if j["id"] not in done]
    todo = len(pending)
    ids = {j["id"] for j in pending}
    print(f"{len(sel)} jobs selected, {len(sel) - len(pending)} already done, P = {a.P}, memory budget {a.mem_gb:g} GB, output {out}")

    def log(line):
        with open(logp, "a") as f:
            f.write(time.strftime("%F_%T ") + line + "\n")

    def append(group, lines):
        p = f"{out}/{group}.csv"
        new = not os.path.exists(p)
        with open(p, "a") as f:
            if new:
                f.write(HEADER)
            f.writelines(lines)

    def launch(j):
        env = {k: v.replace("{out}", out) for k, v in j["env"].items()}
        if j["replay"] or (a.replay and j["algo"].endswith("-reshape")):
            src_group, src_algo = j["replay"] or (j["group"], j["algo"])
            r = rounds_of(f"{out}/{src_group}.csv", j, src_algo) if j["replay"] else None
            if r is None:
                r = rounds_of(os.path.join(J.ROOT, "data", f"{src_group}.csv"), j, src_algo)
            if r is None:
                return None
            env.update(WSPD_LS_ROUNDS=str(r), WSPD_LS_TIME="1000000")
        tag = j["id"].replace(":", "_")
        o, e = open(f"{out}/tmp/{tag}.out", "w"), open(f"{out}/tmp/{tag}.err", "w")
        p = subprocess.Popen(j["cmd"], env={**os.environ, **env}, stdout=o, stderr=e, start_new_session=True)
        o.close(), e.close()
        return p

    running = {}   # pid -> (job, proc, t0, estimate)
    finished = set()
    while pending or running:
        used = sum(r[3] for r in running.values())
        ready = [j for j in pending if j["after"] is None or j["after"] not in ids or j["after"] in finished]
        ready.sort(key=lambda j: -estimate(j, peaks))
        for j in ready:
            if len(running) >= a.P:
                break
            m = estimate(j, peaks)
            if used + m <= a.mem_gb or not running:
                pending.remove(j)
                p = launch(j)
                if p is None:
                    log(f"FAIL {j['id']} (no local-search row to replay)")
                    finished.add(j["id"])
                    continue
                running[p.pid] = (j, p, time.time(), m)
                used += m
        if not running:
            if pending:   # waiting for jobs that are not running: should not happen
                sys.exit(f"{len(pending)} jobs blocked")
            break
        pid, status, ru = os.wait4(-1, 0) if not any(r[0]["timeout"] for r in running.values()) else os.wait4(-1, os.WNOHANG)
        if pid == 0:
            for j, p, t0, _ in running.values():
                if j["timeout"] and time.time() - t0 > j["timeout"]:
                    os.killpg(p.pid, 9)
            time.sleep(1)
            continue
        if pid not in running:
            continue
        j, p, t0, _ = running.pop(pid)
        finished.add(j["id"])
        wall, rc, rss = time.time() - t0, os.waitstatus_to_exitcode(status), ru.ru_maxrss * 1024 / GB
        tag = j["id"].replace(":", "_")
        lines = [l if l.endswith("\n") else l + "\n" for l in open(f"{out}/tmp/{tag}.out") if l.count(",") >= 8]
        if rc == 0 and lines:
            if j["part_out"]:   # the conversion rows go to part_out only
                append(j["part_out"], lines)
                lines = [l for l in lines if not l.split(",")[4].endswith("+part")]
            append(j["group"], lines if j["algo"] != "optx" else lines[-1:])
            c = lines[0].split(",")
            mem = f"{out}/memory.csv"
            new = not os.path.exists(mem)
            with open(mem, "a") as f:
                if new:
                    f.write("out,gen,n,seed,eps,algo,rss_gb,wall_s\n")
                f.write(f"{j['group']},{c[0]},{c[1]},{c[2]},{c[3]},{c[4]},{rss:.2f},{wall:.0f}\n")
            os.remove(f"{out}/tmp/{tag}.out"), os.remove(f"{out}/tmp/{tag}.err")
            bad = "" if all(l.split(",")[6] == "ok" for l in lines) else " INVALID"
            log(f"OK {j['id']} wall={wall:.0f}s rss={rss:.2f}GB{bad}")
        else:
            log(f"{'TIMEOUT' if rc == -9 else 'FAIL'} {j['id']} rc={rc} wall={wall:.0f}s rss={rss:.2f}GB (see {out}/tmp/{tag}.err)")
        print(f"{time.strftime('%T')} {len(finished)}/{todo} {j['id']} rc={rc} {wall:.0f}s", flush=True)


if __name__ == "__main__":
    main()

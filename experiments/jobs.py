"""The experiments of the paper as a list of jobs. Every job is one call of src/exp (or src/optx) with eps = 1/20, and
appends its CSV output to <out>/<group>.csv. The groups are the files of data/:

  exp2, exp2_opt    n = 200..1000 (step 100): the constructions (sweep, CK box, CK exact); exact OPT (CPLEX,
                    1800 s)
  exp3              n = 2000..10000 (step 1000): the constructions
  exp4              n = 20000..50000 (step 10000, seed 1): the constructions
  optx              exact OPT by optx for the spirals with n = 800..1000 (CPLEX alone runs out of memory)
  ls_all_rs         reshape from CK exact and from the sweep on every instance with n <= 10000 (budget 100 s for
                    n <= 1000, 1000 s above), LNS from both starts for n <= 1000; traces RA_*
  ls_cmp_rs         reshape and LNS from both starts on the samples of Fig. ls and Table lstime: 10 instances
                    with n = 1000 (100 s) and 3 with n = 10000 (1000 s); traces R_*
  ls_cmp_big        LNS from both starts on 3 instances with n = 50000 (5000 s, as reshape in ls_all_big_5000); traces RL_*
  ls_all_big_5000   reshape from both starts for n = 20000..50000 (5000 s); traces RL_*; the WSPDs with n = 20000 are
                    dumped for partition_ls_big and partition_ckls_big
  ls_setup          the setup time of the local search (lsSetup in the info column: the time between the end of the
                    construction of the start and the first round), which the traces do not include: every instance
                    from CK exact, and from the sweep for n <= 10000, stopped after one round
  partition_exp     CK exact and the two exact quadtrees, each converted (merge pass), n <= 10000
  partition_exp_big the same for n = 20000
  partition_ls      sweep + reshape converted to a partition WSPD, n <= 10000: the local search is rerun with the
                    number of rounds of the ls_all_rs job of the same instance, which reproduces its WSPD exactly
  partition_ls_big  the dumped 5000 s sweep + reshape WSPDs with n = 20000 (from ls_all_big_5000), converted
  partition_ckls, partition_ckls_big   the same two groups for CK exact + reshape
  partition_sweep   the sweep and its conversion, n <= 10000 and n = 20000
  partition_opt500  the optimal WSPDs with n = 500 and their conversion (run inside the exp2_opt jobs with n = 500)
  sweep_unpruned    the greedy sweep without the redundancy pass (exp ... lazy), n <= 10000 and n = 20000: the lower
                    end 1.503 of the approximation factor of the greedy sweep in the open problem (analysis/text_numbers.py)

Instances: for n <= 10000 the six generators with seeds 1, 2, 3, for n >= 20000 seed 1 only, plus every image-based
CG:SHOP 2019 instance of size n.
"""
import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXP = os.environ.get("WSPD_EXP") or os.path.join(ROOT, "src", "exp")   # WSPD_EXP: a frozen copy of the binary, for long runs
OPTX = os.path.join(ROOT, "src", "optx")
CG = os.path.join(ROOT, "instances", "cgshop2019", "challenge_instances", "data", "images")
GENS = ["uniform", "gauss", "cluster", "grid", "spiral", "annulus"]
SIZES_MID = list(range(200, 1001, 100))
SIZES_LARGE = list(range(2000, 10001, 1000))
SIZES_BIG = list(range(20000, 50001, 10000))
CONSTRUCTIONS = ["lazyprune", "ck", "ckexact"]   # lazyprune = the (greedy) sweep of the paper
QUADTREES = ["qtexact", "qtexactrand"]
CONVERSION = "fewest/min/merge"
OPT_TL = 1800
GROUPS = ["exp2", "exp2_opt", "exp3", "exp4", "optx", "ls_all_rs", "ls_cmp_rs", "ls_cmp_big", "ls_all_big_5000", "ls_setup",
          "partition_exp", "partition_exp_big", "partition_ls", "partition_ls_big", "partition_ckls", "partition_ckls_big",
          "partition_sweep", "sweep_unpruned"]


class Inst:
    """An instance: a generator with n and seed, or a CG:SHOP file (exp gets n = 0, seed = 0)."""

    def __init__(self, gen, n, seed=0):
        self.n, self.seed = n, seed
        if gen in GENS:
            self.gen, self.args, self.name = gen, [gen, str(n), str(seed)], f"{gen}_{n}_{seed}"
        else:   # CG:SHOP family name
            f = f"{gen}-{n:07d}.instance"
            self.gen, self.args, self.name = f, [f"file:{os.path.join(CG, f)}", "0", "0"], f"{gen}-{n:07d}"
        # file name of traces and dumps (as in the paper's runs)
        self.tag = f"{gen}_{n}_{seed}_1_20" if gen in GENS else f"images_{gen}-{n:07d}_instance_0_0_1_20"

    def key(self):   # (gen, seed) as printed by exp; n is printed after the removal of duplicate points
        return self.gen, self.seed


def cgshop(n):
    return sorted(re.sub(r"-\d{7}\.instance$", "", os.path.basename(f)) for f in glob.glob(os.path.join(CG, f"*-{n:07d}.instance")))


def insts(n):
    seeds = (1, 2, 3) if n <= 10000 else (1,)
    return [Inst(g, n, s) for g in GENS for s in seeds] + [Inst(f, n) for f in cgshop(n)]


def budget(n):
    return 100 if n <= 1000 else 1000 if n <= 10000 else 5000


def job(group, inst, algo, env=None, prog=EXP, part_out=None, after=None, replay=None, timeout=None):
    """after: id of a job that must finish first. replay: (group, algo) of the ls job whose number of rounds is
    passed as WSPD_LS_ROUNDS (read when the job starts)."""
    jid = f"{group}:{algo}:{inst.name}"
    env = dict(env or {})
    if algo.startswith("ls-"):
        _, start, moves = algo.split("-")
        env = {"WSPD_LS_START": start, "WSPD_LS_MOVES": moves, **env}
        cmd_algo = "ls"
    else:
        cmd_algo = algo
    cmd = [prog] + inst.args + ["1", "20"] + ([cmd_algo] if prog == EXP else [])
    return dict(id=jid, group=group, inst=inst, algo=algo, cmd=cmd, env=env, part_out=part_out, after=after,
                replay=replay, timeout=timeout)


def ls_job(group, inst, start, moves, T, trace_prefix, env=None, **kw):
    e = {"WSPD_LS_TIME": str(T), "WSPD_LS_TRACE": f"{{out}}/traces/{trace_prefix}_{start}_{moves}_{inst.tag}.csv"}
    return job(group, inst, f"ls-{start}-{moves}", {**e, **(env or {})}, **kw)


def samples():
    s1000 = [Inst(g, 1000, 1) for g in GENS] + [Inst(f, 1000) for f in ("euro-night", "us-night", "paris", "skylake")]
    s10000 = [Inst("uniform", 10000, 1), Inst("cluster", 10000, 1), Inst("euro-night", 10000)]
    s50000 = [Inst("uniform", 50000, 1), Inst("cluster", 50000, 1), Inst("euro-night", 50000)]
    return s1000, s10000, s50000


def make_jobs():
    J = []
    part = {"WSPD_PARTITION": CONVERSION}
    small = SIZES_MID + SIZES_LARGE
    # constructions and exact OPT
    for n in SIZES_MID + SIZES_LARGE + SIZES_BIG:
        group = "exp2" if n <= 1000 else "exp3" if n <= 10000 else "exp4"
        for i in insts(n):
            for a in CONSTRUCTIONS:
                J.append(job(group, i, a))
            if n <= 1000:
                og = "exp2_opt"
                if n == 500:   # the optimal WSPDs with n = 500 are also converted to partition WSPDs
                    J.append(job(og, i, "opt", {"WSPD_OPT_TL": str(OPT_TL), **part}, part_out="partition_opt500",
                                 timeout=OPT_TL + 2200))
                else:
                    J.append(job(og, i, "opt", {"WSPD_OPT_TL": str(OPT_TL)}, timeout=OPT_TL + 2200))
    for n in (800, 900, 1000):
        for s in (1, 2, 3):
            J.append(job("optx", Inst("spiral", n, s), "optx", {"OPTX_THREADS": "3", "OPTX_TL": "2400"}, prog=OPTX))
    # local search
    for n in small:
        for i in insts(n):
            for st in ("ck", "sweep"):
                J.append(ls_job("ls_all_rs", i, st, "reshape", budget(n), "RA"))
                if n <= 1000:
                    J.append(ls_job("ls_all_rs", i, st, "lns", budget(n), "RA"))
    s1000, s10000, s50000 = samples()
    for S, T in ((s1000, 100), (s10000, 1000)):
        for i in S:
            for st in ("ck", "sweep"):
                for mv in ("reshape", "lns"):
                    J.append(ls_job("ls_cmp_rs", i, st, mv, T, "R"))
    for i in s50000:   # LNS with the budget of reshape for n >= 20000; the reshape runs of these instances are in ls_all_big_5000
        for st in ("ck", "sweep"):
            J.append(ls_job("ls_cmp_big", i, st, "lns", 5000, "RL"))
    for n in SIZES_BIG:
        for i in insts(n):
            for st in ("ck", "sweep"):
                env = {"WSPD_DUMP": f"{{out}}/dumps/{i.tag}_{st}"} if n == 20000 else None
                J.append(ls_job("ls_all_big_5000", i, st, "reshape", 5000, "RL", env=env))
    for n in small + SIZES_BIG:
        for i in insts(n):
            for st in ("ck", "sweep") if n <= 10000 else ("ck",):
                J.append(job("ls_setup", i, f"ls-{st}-reshape", {"WSPD_LS_TIME": "1", "WSPD_LS_ROUNDS": "1"}))
    # partition WSPDs (n <= 20000)
    for n in small + [20000]:
        for i in insts(n):
            g = "partition_exp" if n <= 10000 else "partition_exp_big"
            for a in ["ckexact"] + QUADTREES:
                J.append(job(g, i, a, part))
            J.append(job("partition_sweep", i, "lazyprune", part))
            for st, g in (("sweep", "partition_ls"), ("ck", "partition_ckls")):
                a = f"ls-{st}-reshape"
                if n <= 10000:
                    J.append(job(g, i, a, {"WSPD_LS_TIME": "1000000", **part}, after=f"ls_all_rs:{a}:{i.name}",
                                 replay=("ls_all_rs", a)))
                else:
                    J.append(job(g + "_big", i, "load", {"WSPD_LOAD": f"{{out}}/dumps/{i.tag}_{st}.{a}", **part},
                                 after=f"ls_all_big_5000:{a}:{i.name}"))
    # the greedy sweep without the redundancy pass (open problem)
    for n in small + [20000]:
        for i in insts(n):
            J.append(job("sweep_unpruned", i, "lazy"))
    ids = [j["id"] for j in J]
    assert len(ids) == len(set(ids)), "duplicate job ids"
    return J

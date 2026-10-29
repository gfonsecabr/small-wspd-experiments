# Results of the paper's runs

All runs use ε = 1/20. Each CSV file is the output of the experiment group of the same name (see
[`experiments/README.md`](../experiments/README.md) for the instances, budgets, and time limits). Columns:

| Column | Meaning |
| --- | --- |
| `gen` | generator (`uniform`, `gauss`, `cluster`, `grid`, `spiral`, `annulus`) or CG:SHOP file name (`euro-night-0001000.instance`) |
| `n` | number of distinct points (for CG:SHOP, after removing duplicates; the nominal n is in the file name) |
| `seed` | seed of the generator (0 for CG:SHOP) |
| `eps` | ε (0.05) |
| `algo` | algorithm (see [`src/README.md`](../src/README.md)): `lazyprune` (the sweep), `lazy` (the greedy sweep without the redundancy pass), `ck`, `ckexact`, `opt`, `optx`, `ls-<start>-<moves>`, `qtexact`, `qtexactrand`; `<algo>+part` is the conversion of the previous row to a partition WSPD; `load-<algo>` is a WSPD read from a dump |
| `size` | number of dumbbells |
| `valid` | `ok` if the verifier accepts the WSPD (all rows are `ok`) |
| `time` | seconds (for `+part` rows: the conversion only) |
| `info` | details: `steps=` (the lower bound) for the sweep; `optimal=1` and `lb=` for OPT; start size, build time `tBuild`, `rounds`, and move counters for the local search (in `ls_setup.csv` also `lsSetup`, the time between the end of the construction of the start and the first round, `lsFinish`, the relabeling of the result, and `tPrune`, the final redundancy pass); `cover`, `overlap`, `beforeMerge` for the conversions |

| File | Rows |
| --- | --- |
| `exp2.csv`, `exp2_opt.csv` | n = 200..1000: the three constructions; OPT (CPLEX). The OPT runs that ran out of memory (spiral and annulus, n ≥ 800) have no row; a row without `optimal=1` reached the time limit and is not used |
| `exp3.csv`, `exp4.csv` | n = 2000..10000 and n = 20000..50000: the three constructions |
| `optx.csv` | OPT proven by `optx` for the spirals with n = 800, 900, 1000 |
| `ls_all_rs.csv` | local search, every instance with n ≤ 10000 |
| `ls_all_big_5000.csv` | reshape with a 5000 s budget, n ≥ 20000 |
| `ls_cmp_rs.csv`, `ls_cmp_big.csv` | reshape and LNS on the samples of Fig. ls and Table lstime with n = 1000 and 10000; LNS on the sample with n = 50000 (5000 s; its reshape runs are in `ls_all_big_5000.csv`) |
| `ls_setup.csv` | the setup time of the local search (`lsSetup`), which the traces do not include and the runs from CK exact did not record: one run per instance and start, stopped after one round (from the sweep only for n ≤ 10000). The runs from CK exact were measured on 2026-10-11 on the same machine, at most 6 jobs in parallel for n ≤ 10000 and one at a time above; the runs from the sweep are part of the rerun below. Used by Fig. lstotal and `catchup.csv` |
| `partition_*.csv` | partition WSPDs (`partition_ls*`: the sweep followed by reshape and its conversion; `partition_ckls*`: the same for CK exact followed by reshape, run on 2026-10-11); `partition_opt500.csv` holds the optimal WSPDs with n = 500 (dumped by the OPT runs and read back, hence `load-opt`) and their conversion |
| `lb_check.csv` | the lower bounds and their independent check by `src/check_lb` (written by `experiments/certify.py`), one row per instance: `steps` is the number of anchors of the sweep (LB0), `lb1` and `lb` the sizes of the sets of `src/lb_swap` after the swaps and after 30 plateau passes (LB1 and LB; LB is the reference of the paper where OPT is unknown), `check` is `ok` if `check_lb` accepts the three sets, `paper` is `same` if `steps` equals `steps=` of the instance in `exp2.csv`, `exp3.csv`, or `exp4.csv`, `tSweep`, `tSwap0`, `tSwap` are the times in seconds of the sweep and of `lb_swap` until LB1 and in all (Fig. lb), `passes0` the passes before the plateau passes, and `check_s` the time of `check_lb` on the anchors and on LB. All 436 rows are `ok` and `same`. The anchors were checked on 2026-10-12, the sets of `lb_swap` computed and checked on 2026-10-13/14 |
| `alpha.csv` | the largest set of pairwise incompatible pairs (α in the paper) for the 88 instances with n = 200..500 (written by `experiments/alpha.py`, 2026-10-14): `alpha`, `optimal` (1: proven by CPLEX), `ub`, `cols` (maximal dumbbells), `check` (`ok`: accepted by `check_lb`), `time` |
| `catchup.csv` | per instance, the time until reshape from CK exact is smaller than the output of the sweep (`total_s` = `ck_s` + `setup_s` + `reshape_s`, to compare with `sweep_s`; empty if it never is; `setup_s` is the setup of the local search from `ls_setup.csv`), and `sweep_ahead`: the largest (size from CK exact − size from the sweep) / size from the sweep of the two reshape runs at equal total time (time to build the start, plus setup, plus time of reshape; positive if the run from the sweep is ever smaller; for n ≥ 20000, where `ls_setup.csv` has no run from the sweep, the setup from the sweep is the `lsSetup` that the reshape run itself reports). Computed by `analysis/catchup.py` from `ls_setup.csv` and the traces of all reshape runs, most of which are not shipped |
| `sweep_unpruned.csv` | the greedy sweep without the redundancy pass (`lazy`), every instance with n ≤ 10000 and n = 20000 (run on 2026-10-15 with the sources of `src/`; the number of steps of every row equals that of the sweep in `exp2.csv`, `exp3.csv`, `exp4.csv`). Its largest ratio to the smallest WSPD known for the instance, 1.503, is the lower end of the open problem |
| `lns_params.csv` | the tests of the LNS parameters (`experiments/lns_params.py`, 2026-10-12/13; 845 runs): the configuration (`round`, `config`, `K`, `tl`, `cap`, `radius`, `umax`) followed by the row of `exp`, whose `info` also has the counters `lnsOk` (accepted calls), `lnsBig` (calls with more than `umax` uncovered pairs), `lnsFail`, and `lnsCplex` (time in CPLEX), which the sources of `src/` do not print. Round 1 used an earlier neighbor search that recomputed the scale and the centroids of every scanned dumbbell (same moves, slower) |
| `memory.csv` | peak resident memory (`rss_gb`) and wall time of every job of the groups `exp2`–`exp4`, `ls_*` (not OPT) |
| `traces/` | gzipped local-search traces (lines `time,size`: seconds since the start of the local search, size after each change), only those used by the figures and tables: `R_*` for `ls_cmp_rs`, `RL_*` for reshape and LNS on the three instances with n = 50000 of Fig. ls and Table lstime, `RA_sweep_reshape_*` with n ≤ 1000 for the time to reach 1.001 OPT in Fig. time, the other `RA_*` (reshape from CK exact, LNS from both starts) with n = 200 for the first panel of Fig. ls, and `RA_ck_reshape_*` and `RA_sweep_reshape_*` with n = 1000 and n = 10000 for Fig. lstotal. Name: `<prefix>_<start>_<moves>_<instance>_1_20.csv.gz` |

Machine: AMD Ryzen 9 9900X (12 cores, 128 GB of RAM), Fedora Linux 44, GCC 16.2 with `-O2`, CPLEX 22.2; at most
11 single-threaded jobs in parallel (`optx`: 3 threads) and no other load. The runs took place between 2026-09-27 and
2026-10-09, with builds of an earlier version of the sources in `src/`, which had more options; the partition code was
added after the construction, OPT, and local-search runs, and the later builds reproduce those runs exactly (sizes,
steps, and local-search counters at the same number of rounds). The sources in `src/` are that code without the
options the paper does not use, and give identical outputs (see [`src/README.md`](../src/README.md)).
The rows that depend on the sweep (the sweep itself, the local search started from it, the conversions of both, and
`lb_check.csv`) and the local-search sample with n = 50000 come from a rerun of 2026-10-11 and 2026-10-12 on the same
machine, with the sources in `src/` as they are: since then, a step of the greedy sweep takes the closest uncovered
pair of its must-cover set as the next seed, instead of the first one. The LNS rows (`ls-ck-lns`, `ls-sweep-lns`) come from a rerun of 2026-10-13, after the parameters of LNS were tuned (see `src/README.md`).

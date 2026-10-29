# Experiments

- `jobs.py`: every experiment of the paper as a list of jobs, one call of `src/exp` or `src/optx` each, ε = 1/20.
- `run.py`: runs a selection of the jobs in parallel and writes the outputs to `runs/` (never to `data/`).
- `check.py`: compares the outputs in `runs/` with the paper's results in `data/`, row by row.
- `certify.py`: computes the lower bounds (`src/lb_swap`) and checks them with the independent program `src/check_lb` (below).
- `alpha.py`: the largest sets of pairwise incompatible pairs for n ≤ 500 (Fig. lb), by CPLEX, checked by `src/check_lb`.
- `lns_params.py`: the tests of the LNS parameters (below).

```
python3 experiments/run.py --list --groups exp2 --max-n 300               # show the jobs and their memory estimates
python3 experiments/run.py --groups exp2,exp2_opt --max-n 400 -P 8        # run them
python3 experiments/check.py                                              # compare with data/
DATA=runs python3 analysis/make_tables_figs.py                            # tables and figures of your runs
```

Options of `run.py`: `--groups` (comma-separated, default all), `--min-n`, `--max-n` (the nominal n: the n of the
generator or of the CG:SHOP file name), `--only REGEX` (on the job id), `-P` (parallel jobs), `--mem-gb` (memory
budget, default 16), `--out` (default `runs/`), and `--replay` (below). A rerun skips the jobs that `runs/jobs.log`
lists as OK, so an interrupted run can be resumed. The runner starts the jobs with the largest memory estimates first,
as long as the estimates fit in the budget. The estimates come from the peaks measured in the paper's runs
(`data/memory.csv`).

## Groups

Each group appends to the file of the same name, in `runs/` and in `data/`. Instances: for n ≤ 10000 the six
generators with seeds 1, 2, 3, for n ≥ 20000 seed 1 only, plus every image-based CG:SHOP 2019 instance of size n
(n = 200 to 1000 in steps of 100, 2000 to 10000 in steps of 1000, and 20000 to 50000 in steps of 10000).

| Group | Jobs | Contents |
| --- | ---: | --- |
| `exp2`, `exp2_opt` | 594, 198 | n = 200..1000: the three constructions (the sweep (`lazyprune`), CK box, CK exact); OPT (CPLEX, 1800 s) |
| `exp3` | 594 | n = 2000..10000: the three constructions |
| `exp4` | 120 | n = 20000..50000: the three constructions |
| `optx` | 9 | OPT by `optx` for the spirals with n = 800, 900, 1000 |
| `ls_all_rs` | 1188 | reshape from CK exact and from the sweep, every instance with n ≤ 10000 (100 s for n ≤ 1000, 1000 s above); LNS from both starts for n ≤ 1000; traces `RA_*` |
| `ls_cmp_rs` | 52 | reshape and LNS from both starts on the samples of Fig. ls and Table lstime (10 instances with n = 1000, 3 with n = 10000; the n = 200 panel of Fig. ls uses the `ls_all_rs` runs); traces `R_*` |
| `ls_cmp_big` | 6 | LNS from both starts on 3 instances with n = 50000, 5000 s (the reshape runs of these instances are in `ls_all_big_5000`); traces `RL_*` |
| `ls_all_big_5000` | 80 | reshape from both starts, n = 20000..50000, 5000 s; traces `RL_*` |
| `ls_setup` | 832 | the setup time of the local search (`lsSetup` in the `info` column: the time between the end of the construction of the start and the first round, which the traces do not include): the local search stopped after one round, from CK exact on every instance and from the sweep for n ≤ 10000 |
| `partition_exp` | 1188 | CK exact and the two exact quadtrees, each followed by the merge pass, n ≤ 10000 |
| `partition_exp_big` | 30 | the same, n = 20000 |
| `partition_ls` | 396 | sweep + reshape converted to a partition WSPD, n ≤ 10000 (the local search is replayed, see below) |
| `partition_ls_big` | 10 | the 5000 s sweep + reshape WSPDs with n = 20000, converted (needs the `ls_all_big_5000` jobs of these instances in the same output folder, which dump their WSPDs) |
| `partition_ckls` | 396 | CK exact + reshape converted to a partition WSPD, n ≤ 10000 (replayed in the same way) |
| `partition_ckls_big` | 10 | the 5000 s CK exact + reshape WSPDs with n = 20000, converted (dumps as for `partition_ls_big`) |
| `partition_sweep` | 406 | the sweep and its conversion, n ≤ 10000 and n = 20000 |
| `partition_opt500` | (in `exp2_opt`) | the optimal WSPDs with n = 500 and their conversion |
| `sweep_unpruned` | 406 | the greedy sweep without the redundancy pass (`lazy`), n ≤ 10000 and n = 20000: the lower end of the open problem |

Some OPT jobs fail in the paper's setting: CPLEX did not prove optimality within 1800 s for some instances with
n ≥ 600, and the spiral and annulus instances with n ≥ 800 ran out of memory (10 GB per job). These instances use the
lower bound as reference, except the spirals, which `optx` solves.

## Checking the lower bounds

```
python3 experiments/certify.py --max-n 1000 -P 6          # every instance with n <= 1000
python3 experiments/certify.py --only uniform_10000_1     # one instance
```

For every selected instance, `certify.py` runs the sweep with `WSPD_ANCHORS`, which writes the points and the anchors
of the steps, then `src/lb_swap`, which enlarges the set of anchors (LB1 after the swaps, LB after 30 passes with
plateau moves), and `src/check_lb` on the three sets, which verifies that each is pairwise incompatible (see
[`src/README.md`](../src/README.md)). It appends one row per instance to `runs/lb_check.csv`: the three numbers
(`steps`, `lb1`, `lb`), the verdict of `check_lb` (`check`), whether the numbers equal those of the paper's runs in
`data/` (`paper`), the times of the sweep and of `lb_swap`, the number of passes without plateau moves, and the time
of `check_lb`. The exit status is 1 if a check fails or a number differs. Options: `--min-n`, `--max-n`,
`--only REGEX` (on the instance name), `-P`, `--plateau` (default 30), `--mem-gb` (default 16), `--out`, and `--keep`
(keep the point and pair files in `runs/anchors/`). `lb_swap` takes on average 2.7 min for n = 10000 and 82 min for
n = 50000, and 5.2 n² bytes of memory (13 GB for n = 50000); `check_lb` takes seconds. `data/lb_check.csv` holds the
result for the 436 instances of the paper (all `ok` and `same`); its times come from runs with up to 9 instances in
parallel for n ≤ 10000 and 4 above. A rerun of `certify.py --max-n 3000` (242 instances) gave the same three numbers
on every instance.

```
python3 experiments/alpha.py -P 6        # the 88 instances with n = 200..500, at most 2 min each
```

`alpha.py` runs `src/exp ... alpha` (integer program over all maximal dumbbells, at most one chosen pair per dumbbell)
with `WSPD_ALPHA`, which writes the chosen pairs, and checks them with `src/check_lb`. It appends one row per instance
to `runs/alpha.csv`; `data/alpha.csv` holds the paper's results (all proven optimal and checked).

## Exact reproduction and `--replay`

Everything except the local search under a time budget is deterministic: the constructions and their numbers of steps
(the lower bound), the quadtrees, the conversions, and the optima must match `data/` exactly, and `check.py` reports
any difference. The local search stops after a time budget, so the number of moves it makes depends on the machine.
Given the number of rounds, however, reshape is deterministic. With `--replay`, every reshape job stops after the
number of rounds of the same run in `data/` instead of after its time budget, which reproduces the paper's WSPD (size
and counters) exactly. `partition_ls` and `partition_ckls` always do this, with the rounds of `ls_all_rs`, so they convert
exactly the WSPDs of `ls_all_rs`. LNS calls CPLEX under time limits, so it is not replayed; `check.py` reports the ratio of its sizes
to the paper's.

## Cost

The paper's runs took about 400 CPU hours, most of them in the local search for n ≥ 20000 (5000 s per run), in the
local search for n ≤ 10000 and its replays for the conversions, and in the OPT runs that reach the CPLEX time limit. Peak memory grows quadratically
with n. With n ≤ 10000, every job needs at most 2.3 GB, except OPT (more than 10 GB for some instances with n ≥ 800).
With n = 20000, the largest job needs 6 GB. With n = 50000, local search needs 22 GB, the sweep 15 GB, and CK
8 GB.

## Tests of the LNS parameters

```
python3 experiments/lns_params.py 1 --ns 1000,10000 -P 11            # round 1: one parameter at a time
python3 experiments/lns_params.py 2 --ns 1000,10000,50000             # round 2: combinations
python3 experiments/lns_params.py 3 --ns 1000,3000,10000,50000        # round 3: the bound on the uncovered pairs
python3 analysis/lns_params.py                                        # DATA=runs for your runs
```

LNS from CK exact and from the sweep on the local-search samples (n = 1000, 10000, 50000, and three instances with
n = 3000), with the budgets of the paper, for each configuration of the round (the configurations are listed in
`lns_params.py`). The output `runs/lns_params.csv` has the columns of `data/lns_params.csv`. The paper's parameters are
the best ones found. These runs stop after a time budget, so their sizes depend on the machine; round 1 also used an
earlier, slower neighbor search than the one in `src/ls.hpp` (see `lns_params.py`). With 11 parallel jobs, rounds 1
and 2 took 3 and 4 hours. In round 3, the configurations with a bound of 4000 or no bound on the uncovered pairs need
up to 31 GB per run for n = 10000; the machine ran out of memory after 3.5 hours, and the 21 runs that were lost
(listed in `lns_params.py`) were not rerun.

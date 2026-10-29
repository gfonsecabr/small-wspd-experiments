# Certified Small Well-Separated Pair Decompositions in Practice

Code and data for the paper *Certified Small Well-Separated Pair Decompositions in Practice* by Guilherme D. da Fonseca.

Given n points in the plane and ε > 0, a well-separated pair decomposition (WSPD) is a collection of *dumbbells*
{X, Y}, where X and Y have diameter at most ε times their distance, such that every pair of points has one point in X
and the other in Y for some dumbbell. The code computes small WSPDs, certifies their quality, and reproduces every
experiment of the paper (ε = 1/20), in the plane and in three dimensions ([`three_d/`](three_d/README.md)):

- **The sweep**: the greedy sweep of the paper, a certified 40-approximation for ε ≤ 1/16. Its number of steps is a
  lower bound on the minimum size. The lower bound comes with a certificate, the *anchors* of the steps, which the
  separate program `check_lb` verifies.
- **The lower bound LB** of the paper: `lb_swap` enlarges the set of anchors by a local search (1-for-2 swaps and
  plateau moves). The result is again a certificate that `check_lb` verifies, and it certifies the quality of any WSPD.
  Where the optimum is known, LB is at least 0.989 times the optimum.
- **Callahan–Kosaraju** (CK) on a fair split tree, in two variants: *CK box* accepts two tree nodes if their bounding
  boxes are ε-separated, and *CK exact*, or just *CK* as in the paper's tables and figures, if their point sets are.
- **Exact optimum** by set cover over all maximal dumbbells, solved with CPLEX (`exp ... opt`), and `optx`, which adds
  exact set-cover reductions and solves each component separately.
- **Local search**: *reshape*, the heuristic of the paper, and *LNS*, a large neighborhood search with CPLEX used for
  comparison.
- **Partition WSPDs**: conversion of a cover WSPD into a partition WSPD, and the quadtree construction of the paper.

Points have integer coordinates and ε = a/b is rational, so every separation test is exact. Every output is checked
by a separate brute-force verifier, for coverage and exact separation (and, for partition WSPDs, that every pair is covered
exactly once). The lower bounds are checked by `check_lb`, which shares no code with the algorithms (see
[Checking the lower bounds](#checking-the-lower-bounds)).

## Repository layout

| Folder | Contents |
| --- | --- |
| [`src/`](src/README.md) | C++ sources of the four programs: `exp` (all algorithms), `optx` (exact optimum for larger n), `lb_swap` (the lower bound), and `check_lb` (independent check of a lower bound) |
| [`experiments/`](experiments/README.md) | the list of every experiment of the paper, a runner, a check against the paper's results, and the computation and check of the lower bounds |
| [`analysis/`](analysis/README.md) | Python scripts that compute the tables, the figures, and the numbers of the text from the results |
| [`theory/`](theory/README.md) | computer checks used in proofs: admissibility of the piece families, constants of the partition theorem |
| [`data/`](data/README.md) | the raw results of the paper's runs (CSV) and the local-search traces the figures need |
| [`three_d/`](three_d/README.md) | the experiments in three dimensions (appendix *Experiments in Three Dimensions*; Section 5.3 of the arXiv version): separate code, data, and analysis |

## Requirements

- Linux on x86-64, a C++23 compiler (the paper used GCC 16.2 with `-O2`), `git`, `curl`, and `unzip`
- [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) 22.2 (free for academic use through
  the IBM Academic Initiative). Only the exact optimum and LNS call CPLEX, but `exp` and `optx` link against it.
  `lb_swap` and `check_lb` need only a C++17 compiler (`make -C src lb_swap check_lb`).
- Python 3 with `numpy`, `pandas`, `matplotlib`, and `tabulate` (the paper used Python 3.14, pandas 2.3, matplotlib 3.10)
- [pgl](https://github.com/gfonsecabr/pgl), a header-only geometry library, and the
  [CG:SHOP 2019](https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2019/) instances: `./get_deps.sh` downloads both
  (pgl at the commit used for the paper, the instances with a checksum) into `deps/` and `instances/`

## Quick start

```
./get_deps.sh
make -C src                    # set CPLEX=/path/to/CPLEX_Studio222 if CPLEX is not in /opt/ibm/ILOG/CPLEX_Studio222

# tables, figures, and numbers of the paper, from the paper's results in data/  ->  results/
python3 analysis/make_tables_figs.py
python3 analysis/ls_compare.py
python3 analysis/partition_conv.py
python3 analysis/text_numbers.py
python3 analysis/lns_params.py
python3 analysis/inputs_fig.py                # runs src/exp
python3 three_d/analysis/exp3d.py             # the experiments in three dimensions (see three_d/README.md)

# computer checks of the theory part
python3 theory/verify_families.py
python3 theory/partition_bound.py

# rerun the constructions, OPT, and the partition WSPDs with n <= 400 (under a minute with 6 parallel jobs) and the
# local search with n = 200, replayed (about 25 minutes), and compare with the paper's results
python3 experiments/run.py --groups exp2,exp2_opt,partition_exp,partition_sweep --max-n 400 -P 6
python3 experiments/run.py --groups ls_all_rs,partition_ls,partition_ckls --max-n 200 --only reshape --replay -P 6
python3 experiments/check.py

# compute and check the lower bounds of the instances with n <= 1000 (about a minute)
python3 experiments/certify.py --max-n 1000 -P 6
```

One run of `exp` prints one CSV line per algorithm:

```
$ src/exp uniform 1000 1 1 20 lazyprune ckexact
uniform,1000,1,0.05,lazyprune,186537,ok,0.34,steps=173070;...
uniform,1000,1,0.05,ckexact,223185,ok,0.015,
```

(gen, n, seed, ε, algorithm, number of dumbbells, verification, time in seconds, details). Here the sweep finds 186537
dumbbells and certifies that every ε-WSPD of these points has at least 173070, and CK exact finds 223185.

## Where each result of the paper comes from

| Paper | Produced by | From |
| --- | --- | --- |
| Fig. summary, Fig. ratio, Fig. time, Fig. lb, Table exp (appendix) | `analysis/make_tables_figs.py` | `exp2`–`exp4`, `exp2_opt`, `optx`, `ls_all_rs`, `ls_all_big_5000`, `memory`, traces `RA_*` |
| Fig. ls, Fig. lstotal, Table lstime (appendix) | `analysis/ls_compare.py` | `ls_cmp_rs`, `ls_cmp_big`, `ls_all_big_5000`, `ls_all_rs` (n = 200 for Fig. ls; n = 1000 and 10000 for Fig. lstotal), `ls_setup` (Fig. lstotal), traces `R_*`, `RL_*`, `RA_*` |
| Fig. partconv, Table partconv (appendix) | `analysis/partition_conv.py` | `partition_*` |
| Fig. inputs | `analysis/inputs_fig.py` | the generators of `src/exp` and the CG:SHOP instances |
| Numbers quoted in the abstract, the introduction, and Section Experiments | `analysis/text_numbers.py` | all of the above |
| Lower end 1.503 of the open problem (factor of the greedy sweep) | `analysis/text_numbers.py` | `sweep_unpruned` and all cover WSPDs above |
| Choice of the LNS parameters (appendix, paragraph "Local search") | `analysis/lns_params.py` | `lns_params` |
| Table exp3d, Fig. inputs3d, lb3d, summary3d, ratio3d, perpoint3d, time3d, ls3d, lstotal3d, partconv3d, and the numbers of the three-dimensional experiments | `three_d/analysis/exp3d.py` | `three_d/data/` |
| Proposition families, Table families (and the 40 pieces of the sweep theorem) | `theory/verify_families.py` | `theory/families.json` |
| Constants 529 and 252 of the partition theorem | `theory/partition_bound.py` | |

The files of `data/` are the outputs of the experiment groups of the same name (see
[`experiments/README.md`](experiments/README.md)).

## Checking the lower bounds

The sweep proceeds in steps. Each step starts at a pair of points that no earlier dumbbell covers, its *anchor*, and
the number of steps is the lower bound. The bound holds because the anchors are pairwise *incompatible*: no
ε-separated dumbbell covers two of them, so every ε-WSPD needs one dumbbell per anchor. This property can be checked
without trusting the sweep:

```
$ WSPD_ANCHORS=/tmp/u src/exp uniform 1000 1 1 20 lazyprune       # writes /tmp/u.pts and /tmp/u.anchors
uniform,1000,1,0.05,lazyprune,186537,ok,0.34,steps=173070;...
$ src/check_lb /tmp/u.pts /tmp/u.anchors 1 20
OK lower bound 173070: the anchors are pairwise incompatible, n = 1000, eps = 1/20
```

`check_lb` ([`src/check_lb.cpp`](src/check_lb.cpp), one file without dependencies) reads the points and the anchors and
verifies, with exact integer tests, that the anchors are distinct pairs of points and that no two of them are
compatible.

The lower bound LB of the paper is a larger set of pairwise incompatible pairs, computed by `lb_swap`
([`src/lb_swap.cpp`](src/lb_swap.cpp), one file without dependencies) from the anchors, and `check_lb` checks it in
the same way:

```
$ src/lb_swap /tmp/u.pts /tmp/u.anchors 1 20 /tmp/u.lb 30      # writes /tmp/u.lb.swaps (LB1) and /tmp/u.lb (LB)
1000,173070,1409,6554,30,179219,34,1,0,0.59,181063,4,4713,27,0.09
$ src/check_lb /tmp/u.pts /tmp/u.lb 1 20
OK lower bound 181063: the anchors are pairwise incompatible, n = 1000, eps = 1/20
```

`experiments/certify.py` does all this for the instances of the paper and compares the three numbers (steps, LB1,
LB) with those of the paper's runs. The results are in `data/lb_check.csv`: all 436 instances of the paper pass, with
the numbers of the paper (`check_lb` calls every checked pair an anchor). For n ≤ 500, `experiments/alpha.py` computes the largest set of pairwise incompatible pairs
exactly (CPLEX) and checks it with `check_lb` (`data/alpha.csv`).

The sweep itself stops with an error if a step leaves a pair uncovered that it has to cover, which is what would
make its number of steps wrong.

## Reproducing the experiments

`experiments/run.py` runs any subset of the 6515 jobs of the paper and writes the outputs to `runs/` (never to
`data/`). `experiments/check.py` then compares them row by row with `data/`. All outputs except the local search
under a time budget are deterministic and must match exactly: the sizes of all constructions, the lower bounds, the
quadtrees, the conversions, and the optima. The local search is deterministic for a given number of rounds:
`--replay` stops each reshape run after the number of rounds logged in `data/` instead of after its time budget, which
reproduces the paper's local-search outputs exactly on any machine. Running times, the memory figure, and traces
depend on the machine.

The full set takes about 400 CPU hours. The paper's runs used an AMD Ryzen 9 9900X (12 cores, 128 GB of RAM, Fedora
Linux 44, GCC 16.2), with at most 11 single-threaded jobs in parallel and no other load. The largest jobs (n = 50000)
need up to 22 GB of RAM, and CPLEX needs more than 10 GB for some instances with n ≥ 800.

Two sets of runs are not in this job list: the tests of the LNS parameters (`experiments/lns_params.py`, see
[`experiments/README.md`](experiments/README.md#tests-of-the-lns-parameters)) and the experiments in three dimensions
(`three_d/run3d.py`, see [`three_d/README.md`](three_d/README.md)).

## License

MIT (see [`LICENSE`](LICENSE)).

## Citation

```
@misc{Fon26practical,
  author = {Guilherme D. da Fonseca},
  title  = {Certified Small Well-Separated Pair Decompositions in Practice},
  year   = {2026}
}
```

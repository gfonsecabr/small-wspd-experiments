# Experiments in three dimensions

Code, data, and analysis of the experiments in three dimensions of the paper (appendix *Experiments in Three
Dimensions*; Section 5.3 of the arXiv version). Nothing here is shared with the planar code in `../src/`; only the plot style (`../analysis/figutil.py`) is
reused. The paper uses ε = 1/20, as in the plane.

| Path | Contents |
| --- | --- |
| `src/wspd3.hpp` | instances in R^3, exact separation tests, a static kd-tree, the greedy sweep (with the 64 + 32 grid pieces of the higher-dimensional sweep theorem as fixed collection), CK box and CK exact on a fair split tree, the enumeration of the maximal dumbbells, set cover (exact OPT) and the largest set of pairwise incompatible pairs (α) with CPLEX, the redundancy pass, and the verifier |
| `src/ls3.hpp` | local search: reshape (the heuristic of the paper) and LNS (comparison); points relabeled in Morton order |
| `src/partition3.hpp` | conversion of a cover WSPD into a partition WSPD, the merge pass, the octree construction, and the partition verifier |
| `src/exp3d.cpp` | the driver `exp3d` |
| `src/check_lb3.cpp`, `src/lb_swap3.cpp` | the independent check of a lower bound and the lower bound LB of the paper, in R^3 (no dependency) |
| `prep_scans.py` | downloads the four Stanford scans (bunny, armadillo, XYZ RGB dragon, happy Buddha) from the mirror [alecjacobson/common-3d-test-models](https://github.com/alecjacobson/common-3d-test-models) and writes integer point files `data/scans/*.xyz` (largest side of the bounding box 10^6) |
| `run3d.py` | every experiment in three dimensions, as a resumable job queue with a memory budget |
| `analysis/exp3d.py` | the table, the nine figures, and every number of the text (`results/` and `data/summary.md`) |
| `data/` | the results of the paper's runs, the local-search traces the figures need (gzipped), and `summary.md`, every number of the text, written by `analysis/exp3d.py` |
| `results/` | the table and the figures, written by `analysis/exp3d.py` from `data/` (not in the repository; identical to the paper's) |

## Build and quick start

```
make -C three_d/src                 # CPLEX as for ../src (set CPLEX=...); check_lb3 and lb_swap3 need only a compiler
cd three_d && python3 prep_scans.py && cd ..
python3 three_d/analysis/exp3d.py   # results/ from data/ (needs src/exp3d and data/scans for the figure of the inputs)
```

`exp3d` has the interface of `../src/exp`:

```
three_d/src/exp3d GEN n seed a b ALGO [ALGO ...]      (eps = a/b)
```

`GEN` is `uniform`, `gauss`, `cluster`, `grid`, `spiral` (the conical spiral of the paper), `shell` (the spherical
shell), or `file:PATH` of an `x y z` file, of which a random sample of n points is taken with the seed. `ALGO` is
`sweep` (the greedy sweep and the redundancy pass), `ck` (CK box), `ckexact` (CK exact, the paper's "CK"), `opt`,
`alpha`, `ls` (`WSPD_LS_START=ck|sweep`, `WSPD_LS_MOVES=reshape|lns`, `WSPD_LS_TIME`, `WSPD_LS_TRACE`,
`WSPD_LS_ROUNDS`), `qtexact` (the octree construction with exact tests; also `qt`, `qtrand`, `qtexactrand`), or
`load`. `WSPD_PARTITION=fewest/min/merge` converts each output into a partition WSPD (one more line, `<algo>+part`),
`WSPD_ANCHORS=PREFIX` writes the anchors of the sweep for `check_lb3` and `lb_swap3`, and `WSPD_DUMP=PREFIX` writes
the WSPDs. Every output is checked by the verifier (`ok` in the seventh column).

## Rerunning the experiments

The paper uses seed 1 of each generator and one sample of each scan, for n = 200..1000 (step 100) and 2000..10000
(step 1000): 180 instances. From `three_d/`:

```
FIRST= SEEDS=1 PART=1 ONLY=sep,opt,ls,part,lns,alpha DATA3D=runs P=10 MEM_GB=70 python3 run3d.py
DATA3D=runs RESULTS3D=runs/results python3 analysis/exp3d.py    # from the repository root: three_d/analysis/exp3d.py
```

(`DATA3D` needs a link `scans -> ../data/scans`.) This takes about 7 hours on 10 cores and up to 9 GB per job for
n = 10000 (the OPT runs on the conical spirals with n = 1000 take up to 37 GB and 1 hour). Reshape and LNS run for a
fixed time (100 s for n ≤ 1000, 1000 s above), so their results depend on the speed of the machine; all other outputs
are deterministic. A test with n = 200 and 2 s of local search (`MAX_N=200 LS_BUDGET=2`) reproduced the 110 sizes of
the sweep, CK box, CK exact, OPT, α, and the partition WSPDs of CK exact, the sweep, and the octree in `data/`.

The paper's runs were made in two stages (2026-10-28/29), which explains the files in `data/`. The first run
(n = 200, 500, 1000, 2000, 5000, 10000, three seeds) computed the sweep and CK in one process with the lower bound
(`exp3d.csv`, `lb3d.csv`), OPT (`opt3d.csv`), and reshape (`ls3d.csv`). The second one (seed 1) computed, for every
n, the sweep, CK box, and CK exact in separate processes, for their running times and memory (`exp3d_sep.csv`,
`mem3d.csv`), the lower bounds with the times of their stages (`lb3d_sep.csv`), reshape for the other values of n,
converted in the same run (`ls3d.csv`, `part3d.csv`), replays of the first run's reshape runs with the same number of
rounds, converted (`part3d.csv`; they reproduce the first run exactly), the partition WSPDs of CK exact, the sweep, and
the octree (`part3d.csv`), LNS (`lns3d.csv`), α (`alpha3d.csv`), and OPT for the other values of n. Where both stages
computed a value, the analysis checks that they agree. The seeds 2 and 3 of the first run are in the files but not
used. `run3d.log` has every job with its peak memory.

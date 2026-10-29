# C++ sources

`make` builds four programs (see the [top-level README](../README.md) for the requirements):

| File | Contents |
| --- | --- |
| `wspd.hpp` | instances, exact separation tests, the sweep, Callahan–Kosaraju, enumeration of the maximal dumbbells, set cover with CPLEX (exact OPT), the redundancy pass, and the verifier |
| `ls.hpp` | local search on any cover WSPD: the reshape move (the heuristic of the paper) and the LNS move (comparison) |
| `partition.hpp` | conversion of a cover WSPD into a partition WSPD, the merge pass, the quadtree construction, and the partition verifier |
| `inputs.hpp` | the synthetic generators and the CG:SHOP 2019 loader |
| `exp.cpp` | the experiment driver `exp` |
| `optx.cpp` | `optx`: exact OPT with set-cover reductions, for the instances where `exp ... opt` runs out of memory |
| `lb_swap.cpp` | `lb_swap`: the lower bound of the paper, a local search that enlarges the set of anchors of the sweep; it includes none of the files above and needs neither pgl nor CPLEX |
| `check_lb.cpp` | `check_lb`: independent check of a lower bound; it includes none of the files above and needs neither pgl nor CPLEX |

The results in `data/` were produced by an earlier version of these sources that also had options the paper does not
use (other sweep and conversion variants, among them the fixed sweep of the sweep theorem, which is of theoretical
interest only). These options were removed afterwards. The remaining code
gives identical outputs: on a test set of 10 instances with n = 100 to 3000, the dumped WSPDs of every algorithm, of
the conversions, of OPT, and of both local-search moves after a fixed number of rounds are byte-identical, and
`experiments/check.py` finds no difference from `data/`.

## `exp`

```
exp GEN n seed a b ALGO [ALGO ...]        (eps = a/b; the paper uses a = 1, b = 20)
```

`GEN` is a generator (`uniform`, `gauss`, `cluster`, `grid`, `spiral`, `annulus`; the paper's *perturbed grid* is
`grid`, *Gaussian* is `gauss`, and *clustered* is `cluster`) or `file:PATH` for a CG:SHOP 2019 instance (then n and
seed are ignored: pass `0 0`; duplicate points are removed).

| ALGO | Paper |
| --- | --- |
| `lazyprune` | the sweep (greedy sweep followed by the redundancy pass); `info` gives the number of steps, the lower bound (`steps=`) |
| `lazy` | the greedy sweep without the redundancy pass (the open problem); `info` as for `lazyprune` |
| `ck`, `ckexact` | Callahan–Kosaraju, CK box and CK exact |
| `opt` | exact OPT: CPLEX on the set cover over all maximal dumbbells (time limit `WSPD_OPT_TL`, default 1800 s); `info` says whether it is proven optimal |
| `ls` | local search (below) |
| `qtexact`, `qtexactrand` | the quadtree partition WSPD, deterministic or randomly shifted and scaled |
| `load` | reads the WSPD in the file `WSPD_LOAD` (written by `WSPD_DUMP`, for the same instance) |
| `alpha` | not a WSPD: the largest set of pairwise incompatible pairs, a lower bound on OPT, by CPLEX over all maximal dumbbells (at most one chosen pair per dumbbell; time limit `WSPD_OPT_TL`); `size` is its number of pairs and `valid` is `unchecked`: check the pairs written with `WSPD_ALPHA` with `check_lb` (`experiments/alpha.py` does this). Added on 2026-10-14; the other algorithms give the same outputs as before |

Every algorithm prints one line `gen,n,seed,eps,algo,size,valid,time,info`: `n` is the number of distinct points,
`size` the number of dumbbells, `valid` is `ok` if the verifier accepts the WSPD (every pair covered, every dumbbell
exactly ε-separated), and `time` is in seconds, without the verification.

Environment variables:

- `WSPD_LS_START` (`ck`: CK exact, or `sweep`), `WSPD_LS_MOVES` (`reshape` or `lns`), `WSPD_LS_TIME` (budget in
  seconds), `WSPD_LS_TRACE` (file of `time,size` lines, written at every change of the size), `WSPD_LS_ROUNDS` (stop
  after this many rounds: with the same instance and seed, a run with the number of rounds of a timed run ends with
  the same WSPD). The local-search row is `ls-<start>-<moves>`; its `info` gives the start size, the time to build the
  start (`tBuild`), the time to set up the local search (`lsSetup`: relabeling the points and building its tables,
  before the first round; the trace times start after it), and the counters.
- LNS parameters: `WSPD_LS_K` (number of neighbors removed with the target), `WSPD_LS_LNS_RADIUS` (radius of the
  neighborhood, in units of ε times the scale of the target), `WSPD_LS_LNS_UMAX` (new candidate dumbbells are
  generated only if at most this many pairs are left uncovered). The defaults are the values of the paper: radius 2,
  and 16 neighbors and 1000 pairs for n < 5000, 32 neighbors and 400 pairs for n ≥ 5000 (the paper runs LNS for
  n ≤ 1000, n = 10000, and n = 50000). They are the best values of our tests on the local-search samples (`experiments/lns_params.py`), which also
  set `WSPD_LS_LNS_TL` (CPLEX time limit per call, default 2 s) and `WSPD_LS_LNS_CAP` (candidate dumbbells enumerated
  per uncovered pair, default 300).
- `WSPD_PARTITION` (any value) converts the output of every algorithm into a partition WSPD and prints one more row
  `<algo>+part` (`info`: size before the merge pass `beforeMerge`, overlap of the cover WSPD, conversion time).
- `WSPD_DUMP=PREFIX` writes the points to `PREFIX.pts` (one `x y` line per point; the point index is the line
  number, from 0) and every WSPD to `PREFIX.<algo>`, one dumbbell per line: `x1 x2 ...|y1 y2 ...` (point indices).
- `WSPD_ALPHA=PREFIX` writes the points to `PREFIX.pts` and the pairs of `alpha` to `PREFIX.alpha` (one `p q` line per
  pair), the input of `check_lb`.
- `WSPD_ANCHORS=PREFIX` writes the points to `PREFIX.pts` and the anchors of the steps of `lazyprune` to
  `PREFIX.anchors`, one `p q` line per step (point indices), in the order of the steps. These two files are the input
  of `check_lb`.

The sweep stops with an error message (`abort`) if a step leaves uncovered a pair of its must-cover set, the
uncovered pairs after p that are compatible with the anchor {p, q}. Covering this set in every step is what makes the
anchors pairwise incompatible, so a run that ends reports a number of steps that `check_lb` can confirm.

## `check_lb`

```
check_lb POINTS ANCHORS a b [brute]        (eps = a/b)
```

Checks the lower bound of a run of the sweep from the two files written with `WSPD_ANCHORS`. Two pairs of points are
*compatible* if one ε-separated dumbbell covers both. The program verifies that the anchors are distinct pairs of
distinct points and that no two anchors are compatible. Every dumbbell of an ε-WSPD then covers at most one anchor, so
every ε-WSPD of the points, cover or partition, has at least as many dumbbells as there are anchors. It prints
`OK lower bound <number of anchors>: ...` (exit status 0) or `FAIL ...` with the offending anchors (exit status 1).

The pairs {p, q} and {p', q'} are compatible if and only if {{p, p'}, {q, q'}} or {{p, q'}, {q, p'}} is ε-separated,
and the test is exact (128-bit integers). With `brute`, every pair of anchors is tested, which takes time quadratic
in the number of anchors. By default, an anchor {p, q} is tested only against the anchors with one point within
ε|pq| of p and the other within ε|pq| of q, which are the only ones that can be compatible with it. These are found
with a uniform grid on the points and a table of n² bits (0.31 GB for n = 50000). The header of `check_lb.cpp`
gives the argument. The program is independent of the rest: it has its own input code and separation test. It
checks any set of pairs in this format, in particular the output of `lb_swap`.

## `lb_swap`

```
lb_swap POINTS ANCHORS a b OUT [PLATEAU]        (eps = a/b < 1/2; the paper uses PLATEAU = 30)
```

Computes the lower bound LB of the paper. Any set of pairwise incompatible pairs is a lower bound, and the anchors of
the sweep are one such set. Starting from them (the two files written with `WSPD_ANCHORS`), the program adds every
pair that is compatible with no pair of the set, and then repeats *1-for-2 swaps* until a pass over the set finds
none: a pair e of the set is replaced by two pairs that are incompatible with each other and compatible with no pair
of the set other than e, and the pairs that become free are added. This set (LB₁ in the paper) is written to
`OUT.swaps`. Then PLATEAU more passes follow in which a pair without a swap is replaced, with probability 1/2 and a
fixed random seed, by a pair compatible with no other pair of the set, so the size stays the same and later passes
find new swaps. The final set (LB) is written to `OUT`. Both files have the format of the anchors, so `check_lb`
certifies them. The program prints one line
`n,anchors,free,swaps,swapFree,final,passes,converged,capHits,seconds,LB,passes0,swaps0,swapFree0,seconds0`, where
`final` is LB₁ and the fields after `seconds` describe the end of the passes without plateau moves. It keeps the
points sorted by distance from every point (4n² bytes) and one counter per pair (n² bytes): 13 GB for n = 50000. The
search for the two pairs of a swap tests at most 200000 candidate pairs per pair of the set and pass (`capHits`
counts the pairs where the limit was reached). The header of `lb_swap.cpp` describes the algorithm in more detail.

## `optx`

```
optx GEN n seed a b
```

The same set-cover instance as `exp ... opt`, reduced to a fixpoint by exact reductions (forced columns, column
dominance, row dominance) and split into connected components, each solved by CPLEX with MIP gap 0, from a start
computed by the sweep followed by 60 s of reshape. Environment: `OPTX_THREADS` (enumeration and CPLEX threads),
`OPTX_TL` (CPLEX seconds for the whole run). It prints one CSV line with algorithm `optx`; `info` says whether the
result is proven optimal. The paper's references for the spirals with n = 800, 900, 1000 come from `optx`
(`data/optx.csv`).

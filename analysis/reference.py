"""Shared paths and the reference value of an instance.

DATA (environment variable, default data/) is the folder with the experiment outputs; RESULTS (default results/) is
where the tables, figures, and summaries are written. To analyze your own runs: DATA=runs python3 analysis/...

Reference of an instance: OPT where it is proven optimal (CPLEX in exp2_opt.csv, or optx in optx.csv),
otherwise the lower bound LB of the paper (column lb of lb_check.csv: the anchors of the sweep enlarged by src/lb_swap,
certified by src/check_lb; experiments/certify.py). A CPLEX solution that is not proven optimal is never used."""
import glob
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.abspath(os.environ.get("DATA", os.path.join(ROOT, "data")))
RESULTS = os.path.abspath(os.environ.get("RESULTS", os.path.join(ROOT, "results")))
for sub in ("tables", "figs"):
    os.makedirs(os.path.join(RESULTS, sub), exist_ok=True)
KEY = ["gen", "n", "seed", "eps"]
MIN_N = 200   # smallest n of the paper
MAX_N = 50000   # largest n of the paper


def read(*names):
    """The CSV files of DATA that exist, concatenated."""
    dfs = [pd.read_csv(os.path.join(DATA, f)) for f in names
           if os.path.exists(os.path.join(DATA, f)) and os.path.getsize(os.path.join(DATA, f)) > 0]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def trace_files(pattern):
    """Trace files of DATA/traces matching pattern + '.csv' or pattern + '.csv.gz'."""
    return sorted(glob.glob(os.path.join(DATA, "traces", pattern + ".csv")) +
                  glob.glob(os.path.join(DATA, "traces", pattern + ".csv.gz")))


def with_optx(df):
    """df plus the optima proven by optx (DATA/optx.csv, as rows with algo "opt" and info "solver=optx;...") for the
    instances of df that have no proven CPLEX optimum. Their running times are not CPLEX times: rows with "solver=optx"
    in info are excluded from the timings of the exact solver."""
    x = read("optx.csv")
    if not len(df) or not len(x): return df
    x = x[(x.algo == "optx") & (x.valid == "ok") & x["info"].str.contains("optimal=1")].copy()
    x["info"] = "solver=optx;" + x["info"]
    x["algo"] = "opt"
    x = x.drop_duplicates(KEY, keep="last")
    inst = df.drop_duplicates(KEY).set_index(KEY).index
    o = df[df.algo == "opt"]
    proven = o[o["info"].astype(str).str.contains("optimal=1")].set_index(KEY).index
    k = x.set_index(KEY).index
    x = x[k.isin(inst) & ~k.isin(proven)]
    if not len(x): return df
    drop = (df.algo == "opt") & df.set_index(KEY).index.isin(x.set_index(KEY).index)   # unproven CPLEX rows
    return pd.concat([df[~drop], x], ignore_index=True)


def proven_opt(opt_rows):
    """opt_rows: rows with algo == 'opt'. Series indexed by KEY: the size where proven optimal, else NaN."""
    o = opt_rows.drop_duplicates(KEY).set_index(KEY)
    return o["size"].where(o["info"].str.contains("optimal=1")).rename("opt")


def lower_bounds():
    """DataFrame indexed by KEY (CG:SHOP: nominal n, seed 0), from DATA/lb_check.csv (experiments/certify.py): steps =
    number of steps of the sweep (LB0 in the paper), lb1 = after the 1-for-2 swaps (LB1), lb = after the plateau passes
    (LB, the lower bound of the paper), and the times tSweep, tSwap0 (lb1), tSwap (lb), and passes0. Only certified rows."""
    L = read("lb_check.csv")
    L = L[L.check == "ok"].drop_duplicates(KEY, keep="last").set_index(KEY)
    return L[["steps", "lb1", "lb", "tSweep", "tSwap0", "tSwap", "passes0"]]


def sweep_lb(st):
    """st: the lazyprune rows of the batch files, indexed by KEY (unique). Series indexed by KEY: the lower bound LB
    (lower_bounds().lb), which must exist for every instance of the paper (eps = 1/20, MIN_N <= n <= MAX_N), and where
    the number of steps of the row must equal that of lb_check.csv; elsewhere (not in the paper) the number of steps."""
    steps = st["info"].str.extract(r"steps=(\d+)")[0].astype(float)
    L = lower_bounds().reindex(st.index)
    idx = st.index.to_frame(index=False)
    paper = ((idx.eps == 0.05) & (idx.n >= MIN_N) & (idx.n <= MAX_N)).values
    miss = paper & L.lb.isna().values
    if miss.any(): raise SystemExit("lb_check.csv has no certified lower bound for %d instances, e.g. %s" % (miss.sum(), st.index[miss][:3].tolist()))
    have = L.lb.notna().values
    bad = have & (L.steps.values != steps.values)
    if bad.any(): raise SystemExit("the steps of lb_check.csv differ from the sweep rows on %s" % st.index[bad][:3].tolist())
    return L.lb.where(have, steps).rename("lb")


def refs():
    """(opt, ref) indexed by KEY: proven OPT (NaN if unknown) and the reference (OPT, else the lower bound LB)."""
    b = with_optx(read("exp2.csv", "exp2_opt.csv", "exp3.csv", "exp4.csv"))
    b = b[b.valid.astype(str).str.startswith("ok")]
    s = b[b.algo == "lazyprune"].drop_duplicates(KEY).set_index(KEY)
    lb = sweep_lb(s)
    opt = proven_opt(b[b.algo == "opt"])
    ref = opt.reindex(opt.index.union(lb.index)).fillna(lb).rename("lb")
    return opt, ref

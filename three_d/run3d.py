"""Runs every experiment in R^3 (run from three_d/; resumes: jobs logged OK are not rerun). Outputs in data/:
  exp3d.csv  sweep, CK box, CK exact (exp3d ... sweep ck ckexact, with WSPD_ANCHORS)
  lb3d.csv   lower bounds: steps of the sweep (LB0), lb_swap3 without (LB1) and with 30 plateau passes (LB), each set
             certified by check_lb3:  inst,gen,n,seed,steps,lb1,lb,check,tSwap
  opt3d.csv  exact OPT (CPLEX, 1800 s) for n <= OPT_MAX_N; a run is stopped after OPT_WALL seconds (enumeration of the
             maximal dumbbells included) and then logged OK with "timeout" and no row (the reference stays LB)
  ls3d.csv   reshape from CK exact and from the sweep, 100 s for n <= 1000 and 1000 s above (traces in data/traces/)
  run3d.log  one line per job (START / OK / FAIL), with the peak memory
Kinds added for the second run (2026-10-28/29: seed 1 only, and a 3D analogue of every planar plot; binary exp3d
with LNS, alpha, the partition conversion, and the octree, which reproduces the first run's outputs):
  sep        sweep (with anchors), CK box, CK exact as three processes (time and peak memory of each, mem3d.csv:
             inst,algo,peakGB,wallS), then lb_swap3 (its whole output line) and check_lb3 on the three sets
             -> exp3d_sep.csv, lb3d_sep.csv (inst,gen,n,seed,check,<the output line of lb_swap3>)
  ls         as above; with PART, the result is also converted (WSPD_PARTITION, its row goes to part3d.csv)
  replay     reshape rerun with the rounds of its ls3d.csv row (reproduces it exactly; checked) and converted -> part3d.csv
  part       CK exact, the sweep, and the exact octree, each converted -> part3d.csv
  lns        LNS from both starts on the local-search samples (all instances with n = 200; seed 1 of every generator and
             the four scans with n = 1000; uniform, cluster, bunny with n = 10000), budgets as reshape -> lns3d.csv,
             traces lns_<start>_<inst>.csv
  alpha      the largest set of pairwise incompatible pairs, n <= 500 -> alpha3d.csv (stopped after OPT_WALL seconds)
Environment: P (parallel jobs, default 10), MEM_GB (memory budget for the reservations, default 80), DRY=1 (list the
jobs), ONLY=kind[,kind] (base, opt, ls, sep, replay, part, lns, alpha), SEEDS (seeds of the synthetic inputs, default
1,2,3), PART=1 (convert in the ls jobs), FIRST (see below), MAX_N; for tests: DATA3D (output folder, with a link scans -> data/scans),
LS_BUDGET (seconds).
Long runs: systemd-run --user --unit=wspd3d -p MemoryMax=95G --working-directory=$PWD python3 run3d.py   (from three_d/)"""
import os
import resource
import subprocess
import sys
import threading
import time

D = os.environ.get('DATA3D', 'data')   # another folder for tests
WORK = f'{D}/work'
TR = f'{D}/traces'
GENS = ['uniform', 'gauss', 'cluster', 'grid', 'spiral', 'shell']
SCANS = ['bunny', 'armadillo', 'dragon', 'buddha']
# the values of n of the planar experiments up to 10000 (first run, 2026-10-28: 200, 500, 1000, 2000, 5000, 10000)
SIZES = list(range(200, 1001, 100)) + list(range(2000, 10001, 1000))
OPT_MAX_N = int(os.environ.get('OPT_MAX_N', 1000))
OPT_WALL = int(os.environ.get('OPT_WALL', 7200))
MAX_N = int(os.environ.get('MAX_N', 10000))
P = int(os.environ.get('P', 10))
MEM = float(os.environ.get('MEM_GB', 80))
SEEDS = [int(x) for x in os.environ.get('SEEDS', '1,2,3').split(',')]
PART = os.environ.get('PART') == '1'
SPEC = 'fewest/min/merge'
# the sizes of the first run (their reshape runs are replayed for the conversion); FIRST= (empty): none, so that the ls
# jobs of a complete rerun convert at every size
FIRST = tuple(int(x) for x in os.environ.get('FIRST', '200,500,1000,2000,5000,10000').split(',') if x)
LNS_SAMPLE = {200: None, 1000: [f'{g}_1000_1' for g in GENS] + [f'{f}_1000_1' for f in SCANS],
              10000: ['uniform_10000_1', 'cluster_10000_1', 'bunny_10000_1']}   # None: every instance
EXP, CHECK, SWAP = 'src/exp3d', 'src/check_lb3', 'src/lb_swap3'


def instances():
    for n in SIZES:
        if n > MAX_N:
            continue
        for g in GENS:
            for s in SEEDS:
                yield dict(name=f'{g}_{n}_{s}', gen=g, n=n, seed=s)
        for f in SCANS:
            yield dict(name=f'{f}_{n}_1', gen=f'file:{D}/scans/{f}.xyz', n=n, seed=1)


def budget(n):
    if os.environ.get('LS_BUDGET'):   # tests only
        return float(os.environ['LS_BUDGET'])
    return 100 if n <= 1000 else 1000


def mem_gb(kind, n, gen=''):
    # measured at n = 10000 (first run): base 4.5 GB on uniform points, reshape up to 8.4 GB, lb_swap3 0.7 GB; OPT at
    # n = 1000: 9.4 GB on cluster, 30 and 37 GB on the spirals (its enumeration of maximal dumbbells is far larger)
    # second run: the conversion at n = 10000 (grid) and LNS measured before the launch (see run3d.log)
    if kind in ('opt', 'alpha'):
        return 0.3 + (38 * (n / 1000) ** 3 if gen.startswith('spiral') else 12 * (n / 1000) ** 2)
    if kind == 'ls' and PART:
        kind = 'replay'
    return 0.3 + {'base': 5.0, 'ls': 9.0, 'sep': 9.0, 'replay': 14.0, 'part': 14.0, 'lns': 10.0}[kind] * (n / 10000) ** 2


def jobs():
    J = []
    only = os.environ.get('ONLY', 'base,opt,ls').split(',')
    rounds = {}   # (start, inst) -> the row of ls3d.csv (for replay)
    if 'replay' in only:
        for line in open(f'{D}/ls3d.csv'):
            f = line.strip().split(',')
            g = f[0]
            rounds[(f[4].split('-')[1], f'{g}_{f[1]}_{f[2]}')] = line.strip()
    for i in instances():
        n = i['n']
        if 'sep' in only:
            J.append(dict(id=f"sep:{i['name']}", kind='sep', inst=i, est=n * n))
        if 'part' in only:
            J.append(dict(id=f"part:{i['name']}", kind='part', inst=i, est=n * n))
        if 'alpha' in only and n <= 500:
            J.append(dict(id=f"alpha:{i['name']}", kind='alpha', inst=i, est=n ** 3))
        if 'replay' in only and n in FIRST:
            for st in ('ck', 'sweep'):
                J.append(dict(id=f"replay:{st}:{i['name']}", kind='replay', start=st, inst=i, row=rounds[(st, i['name'])],
                              est=budget(n) * 1e6 + n * n))
        if 'lns' in only and n in LNS_SAMPLE and (LNS_SAMPLE[n] is None or i['name'] in LNS_SAMPLE[n]):
            for st in ('ck', 'sweep'):
                J.append(dict(id=f"lns:{st}:{i['name']}", kind='lns', start=st, inst=i, est=budget(n) * 1e6 + n * n))
        if 'base' in only:
            J.append(dict(id=f"base:{i['name']}", kind='base', inst=i, est=n * n))
        if 'opt' in only and n <= OPT_MAX_N:
            J.append(dict(id=f"opt:{i['name']}", kind='opt', inst=i, est=n ** 3))
        if 'ls' in only and not (PART and n in FIRST):   # the second run: only the sizes not in the first
            for st in ('ck', 'sweep'):
                J.append(dict(id=f"ls:{st}:{i['name']}", kind='ls', start=st, inst=i, est=budget(n) * 1e6 + n * n))
    J.sort(key=lambda j: -j['est'])   # longest first
    return J


lock = threading.Lock()


def log(msg):
    with lock, open(f'{D}/run3d.log', 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S ') + msg + '\n')


def append(path, line):
    with lock, open(path, 'a') as f:
        f.write(line.rstrip('\n') + '\n')


def run(cmd, env=None):
    """runs cmd, returns (returncode, stdout, peak RSS of the child in GB)"""
    e = dict(os.environ)
    e.update(env or {})
    p = subprocess.run(['/usr/bin/time', '-f', 'PEAKKB %M'] + cmd, env=e, capture_output=True, text=True)
    peak = 0.0
    for line in p.stderr.splitlines():
        if line.startswith('PEAKKB'):
            peak = int(line.split()[1]) / 2 ** 20
    return p.returncode, p.stdout, peak, p.stderr


def do_job(j):
    i = j['inst']
    args = [i['gen'], str(i['n']), str(i['seed']), '1', '20']
    if j['kind'] == 'base':
        pre = f"{WORK}/{i['name']}"
        rc, out, peak, err = run([EXP] + args + ['sweep', 'ck', 'ckexact'], {'WSPD_ANCHORS': pre})
        if rc != 0 or out.count('\n') != 3 or ',ok,' not in out:
            return False, f'exp3d rc={rc} {err[-300:]}', peak
        steps = int(out.split('steps=')[1].split(';')[0])
        rc, sw, pk2, err = run([SWAP, pre + '.pts', pre + '.anchors', '1', '20', pre + '.lb', '30'])
        peak = max(peak, pk2)
        if rc != 0:
            return False, f'lb_swap3 rc={rc} {sw} {err[-300:]}', peak
        f = sw.strip().split(',')
        if int(f[1]) != steps:
            return False, f'anchors {f[1]} != steps {steps}', peak
        lb1, lb, tswap = int(f[5]), int(f[10]), float(f[9])
        checks = []
        for path, val in ((pre + '.anchors', steps), (pre + '.lb.swaps', lb1), (pre + '.lb', lb)):
            rc, ck, pk3, err = run([CHECK, pre + '.pts', path, '1', '20'])
            peak = max(peak, pk3)
            checks.append(rc == 0 and ck.startswith(f'OK lower bound {val}:'))
        for line in out.splitlines():
            append(f'{D}/exp3d.csv', line)
        append(f'{D}/lb3d.csv', f"{i['name']},{args[0].split('/')[-1].split('.')[0]},{i['n']},{i['seed']},{steps},{lb1},{lb},"
                                 f"{'ok' if all(checks) else 'FAIL'},{tswap}")
        for ext in ('.pts', '.anchors', '.lb', '.lb.swaps'):
            try:
                os.remove(pre + ext)
            except FileNotFoundError:
                pass
        return all(checks), f'steps={steps} lb1={lb1} lb={lb}', peak
    if j['kind'] == 'opt':
        rc, out, peak, err = run(['timeout', str(OPT_WALL), EXP] + args + ['opt'], {'WSPD_OPT_TL': '1800'})
        if rc == 124:
            return True, f'timeout after {OPT_WALL} s', peak
        if rc != 0 or ',ok,' not in out:
            return False, f'rc={rc} {err[-300:]}', peak
        append(f'{D}/opt3d.csv', out)
        return True, out.split(',')[5] + ' ' + out.split(',')[8].split(';')[1], peak
    if j['kind'] == 'sep':
        pre = f"{WORK}/{i['name']}"
        rows, mem = [], []
        for a in ('sweep', 'ck', 'ckexact'):
            rc, out, peak, err = run([EXP] + args + [a], {'WSPD_ANCHORS': pre} if a == 'sweep' else {})
            if rc != 0 or out.count('\n') != 1 or ',ok,' not in out:
                return False, f'exp3d {a} rc={rc} {err[-300:]}', peak
            rows.append(out.strip()); mem.append(f"{i['name']},{a},{peak:.3f},{float(out.split(',')[7]):.3f}")
        steps = int(rows[0].split('steps=')[1].split(';')[0])
        rc, sw, pk2, err = run([SWAP, pre + '.pts', pre + '.anchors', '1', '20', pre + '.lb', '30'])
        if rc != 0:
            return False, f'lb_swap3 rc={rc} {sw} {err[-300:]}', pk2
        f = sw.strip().split(',')
        if int(f[1]) != steps:
            return False, f'anchors {f[1]} != steps {steps}', pk2
        mem.append(f"{i['name']},lb_swap3,{pk2:.3f},{float(f[9]):.3f}")
        checks = []
        for path, val in ((pre + '.anchors', steps), (pre + '.lb.swaps', int(f[5])), (pre + '.lb', int(f[10]))):
            rc, ck, pk3, err = run([CHECK, pre + '.pts', path, '1', '20'])
            checks.append(rc == 0 and ck.startswith(f'OK lower bound {val}:'))
        for r in rows:
            append(f'{D}/exp3d_sep.csv', r)
        for r in mem:
            append(f'{D}/mem3d.csv', r)
        append(f'{D}/lb3d_sep.csv', f"{i['name']},{args[0].split('/')[-1].split('.')[0]},{i['n']},{i['seed']},{'ok' if all(checks) else 'FAIL'},{sw.strip()}")
        for ext in ('.pts', '.anchors', '.lb', '.lb.swaps'):
            try:
                os.remove(pre + ext)
            except FileNotFoundError:
                pass
        return all(checks), f'steps={steps} lb1={f[5]} lb={f[10]}', max(m for m in [float(x.split(',')[2]) for x in mem])
    if j['kind'] == 'part':
        rc, out, peak, err = run([EXP] + args + ['ckexact', 'sweep', 'qtexact'], {'WSPD_PARTITION': SPEC})
        if rc != 0 or out.count('\n') != 6 or out.count(',ok,') != 6 or 'partition=1' not in out:
            return False, f'rc={rc} {out[-300:]} {err[-300:]}', peak
        for line in out.splitlines():
            append(f'{D}/part3d.csv', line)
        return True, ' '.join(l.split(',')[5] for l in out.splitlines()), peak
    if j['kind'] == 'alpha':
        rc, out, peak, err = run(['timeout', str(OPT_WALL), EXP] + args + ['alpha'], {'WSPD_OPT_TL': '1800'})
        if rc == 124:
            return True, f'timeout after {OPT_WALL} s', peak
        if rc != 0 or ',alpha,' not in out:
            return False, f'rc={rc} {err[-300:]}', peak
        append(f'{D}/alpha3d.csv', out)
        return True, out.split(',')[5] + ' ' + out.split(',')[8].split(';')[1], peak
    st = j['start']
    if j['kind'] == 'lns':
        env = {'WSPD_LS_START': st, 'WSPD_LS_MOVES': 'lns', 'WSPD_LS_TIME': str(budget(i['n'])), 'WSPD_LS_TRACE': f"{TR}/lns_{st}_{i['name']}.csv"}
        rc, out, peak, err = run([EXP] + args + ['ls'], env)
        if rc != 0 or ',ok,' not in out:
            return False, f'rc={rc} {err[-300:]}', peak
        append(f'{D}/lns3d.csv', out)
        return True, out.split(',')[5], peak
    if j['kind'] == 'replay':
        import re
        ref = j['row'].split(',')
        rnd = re.search(r'rounds=(\d+)', ref[8]).group(1)
        env = {'WSPD_LS_START': st, 'WSPD_LS_TIME': '1000000', 'WSPD_LS_ROUNDS': rnd, 'WSPD_PARTITION': SPEC}
        rc, out, peak, err = run([EXP] + args + ['ls'], env)
        lines = out.strip().splitlines()
        if rc != 0 or len(lines) != 2 or out.count(',ok,') != 2:
            return False, f'rc={rc} {err[-300:]}', peak
        key = lambda info: re.sub(r';(tBuild|lsSetup|tPrune|lsTime)=[0-9.]+', '', info)
        got = lines[0].split(',')
        if got[5] != ref[5] or key(got[8]) != key(ref[8]):
            return False, f'replay differs: {got[5]} {got[8]} vs {ref[5]} {ref[8]}', peak
        for line in lines:
            append(f'{D}/part3d.csv', line)
        return True, f'reproduced {ref[5]}, partition {lines[1].split(",")[5]}', peak
    env = {'WSPD_LS_START': st, 'WSPD_LS_TIME': str(budget(i['n'])), 'WSPD_LS_TRACE': f"{TR}/{st}_{i['name']}.csv"}
    if PART:
        env['WSPD_PARTITION'] = SPEC
    rc, out, peak, err = run([EXP] + args + ['ls'], env)
    lines = out.strip().splitlines()
    if rc != 0 or len(lines) != (2 if PART else 1) or out.count(',ok,') != len(lines):
        return False, f'rc={rc} {err[-300:]}', peak
    append(f'{D}/ls3d.csv', lines[0])
    if PART:
        append(f'{D}/part3d.csv', lines[1])
    return True, out.split(',')[5], peak


def main():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(TR, exist_ok=True)
    J = jobs()
    done = set()
    if os.path.exists(f'{D}/run3d.log'):
        for line in open(f'{D}/run3d.log'):
            w = line.split()
            if len(w) >= 4 and w[2] == 'OK':
                done.add(w[3])
    J = [j for j in J if j['id'] not in done]
    if os.environ.get('DRY'):
        for j in J:
            print(j['id'], f"{mem_gb(j['kind'], j['inst']['n'], j['inst']['gen']):.1f} GB")
        print(len(J), 'jobs')
        return
    log(f'RUN {len(J)} jobs, P = {P}, MEM_GB = {MEM}')
    used = [0.0]
    cv = threading.Condition()
    pending = list(J)

    def worker():
        while True:
            with cv:
                while True:
                    if not pending:
                        return
                    k = next((k for k, j in enumerate(pending) if used[0] + mem_gb(j['kind'], j['inst']['n'], j['inst']['gen']) <= MEM), None)
                    if k is not None or used[0] == 0:
                        j = pending.pop(k or 0)
                        m = mem_gb(j['kind'], j['inst']['n'], j['inst']['gen'])
                        used[0] += m
                        break
                    cv.wait()
            log(f"START {j['id']} ({m:.1f} GB)")
            t0 = time.time()
            try:
                ok, info, peak = do_job(j)
            except Exception as ex:   # noqa: BLE001
                ok, info, peak = False, repr(ex), 0
            log(f"{'OK' if ok else 'FAIL'} {j['id']} {time.time() - t0:.0f}s peak={peak:.2f}GB {info}")
            with cv:
                used[0] -= m
                cv.notify_all()

    T = [threading.Thread(target=worker) for _ in range(P)]
    for t in T:
        t.start()
    for t in T:
        t.join()
    log('END')


if __name__ == '__main__':
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    main()

// Exact minimum WSPD for larger n: the set-cover ILP of `exp ... opt` (all maximal separated pairs), with
//   1. enumeration of the maximal pairs on several threads,
//   2. exact set-cover reductions to a fixpoint: forced columns (an element with one column), column dominance
//      (a column whose live elements lie in another column), row dominance (if every live column of e covers f,
//      drop f), and empty columns,
//   3. the remaining instance split into connected components, each solved by CPLEX (MIP gap 0),
//   4. a MIP start from the sweep improved by reshape (60 s), mapped to maximal columns.
// usage: optx GEN n seed a b        (same instances as exp; eps = a/b)
// Environment: OPTX_THREADS (enumeration and CPLEX threads, default 8), OPTX_TL (CPLEX seconds for the whole run,
//   default 3600).
// Prints one CSV line: gen,n,seed,eps,optx,size,valid,time,info  (info: optimal, lb, columns, reductions, components).
#include <atomic>
#include <iomanip>
#include <numeric>
#include <thread>
#include "wspd.hpp"
#include "ls.hpp"
#include "inputs.hpp"
using namespace wspd;

static double T0() { static auto t = std::chrono::steady_clock::now(); return seconds_since(t); }
#define LOG(x) (std::cerr << "[" << std::fixed << std::setprecision(1) << T0() << "s] " << x << std::endl)

// all maximal separated pairs, anchors (p,q) split over threads; each thread deduplicates, then a global merge
static std::vector<Pair> allMaximalPairsPar(const Instance& I, const Spatial& S, int threads, bool& capHit) {
    int n = I.n(); std::atomic<int> next{0}; std::atomic<bool> hit{false};
    std::vector<std::vector<Pair>> part(threads);
    auto work = [&](int t) {
        std::unordered_set<Pair, PairHash, PairEq> seen; std::vector<Pair> tmp;
        for (int p; (p = next++) < n; )
            for (int q = 0; q < n; ++q) if (p != q) {
                tmp.clear(); anchoredMaximalPairs(I, S, p, q, tmp);
                if (tmp.size() >= 200000) hit = true;
                for (auto& pr : tmp) seen.insert(canonical(pr));
            }
        part[t].assign(seen.begin(), seen.end());
    };
    std::vector<std::thread> th; for (int t = 0; t < threads; ++t) th.emplace_back(work, t);
    for (auto& x : th) x.join();
    std::unordered_set<Pair, PairHash, PairEq> seen; std::vector<Pair> out;
    for (auto& v : part) { for (auto& pr : v) if (seen.insert(pr).second) out.push_back(std::move(pr)); v.clear(); v.shrink_to_fit(); }
    capHit = hit;
    return out;
}

int main(int argc, char** argv) {
    if (argc < 6) { std::cerr << "usage: optx GEN n seed a b" << std::endl; return 1; }
    std::string gen = argv[1]; int nArg = std::atoi(argv[2]); unsigned seed = std::atoi(argv[3]);
    Eps eps{std::atoll(argv[4]), std::atoll(argv[5])};
    auto env = [](const char* k, const char* d) { const char* v = std::getenv(k); return std::string(v ? v : d); };
    int threads = std::atoi(env("OPTX_THREADS", "8").c_str());
    double TL = std::atof(env("OPTX_TL", "3600").c_str());
    T0();
    Instance I = generate(gen, nArg, seed, eps);
    int n = I.n(); ElemIndex EI(n); size_t M = EI.size();
    Spatial S(I);

    // ---- heuristic start (runs first: it is cheap and gives the reference size)
    std::vector<Pair> H;
    {
        SweepStats st; H = pruneRedundant(I, sweep(I, st));
        LSConfig cf; cf.timeLimit = 60; cf.seed = seed + 777; LSResult R;
        H = pruneRedundant(I, runLocalSearch(I, std::move(H), cf, R));
    }
    if (!verify(I, H)) { LOG("start WSPD is not valid"); return 1; }
    LOG("start size " << H.size());

    // ---- columns and incidence
    bool capHit = false;
    auto cols = allMaximalPairsPar(I, S, threads, capHit);
    size_t C = cols.size();
    LOG("columns " << C << (capHit ? " (anchor cap hit: NOT exact)" : ""));
    std::vector<size_t> cstart(C + 1, 0);
    for (size_t c = 0; c < C; ++c) cstart[c + 1] = cstart[c] + cols[c].X.size() * cols[c].Y.size();
    std::vector<uint32_t> celem(cstart[C]);
    std::vector<uint32_t> rcnt(M + 1, 0);
    for (size_t c = 0; c < C; ++c) { size_t k = cstart[c];
        for (int u : cols[c].X) for (int v : cols[c].Y) { uint32_t e = (uint32_t)EI(u, v); celem[k++] = e; ++rcnt[e]; } }
    std::vector<size_t> rstart(M + 1, 0);
    for (size_t e = 0; e < M; ++e) rstart[e + 1] = rstart[e] + rcnt[e];
    std::vector<uint32_t> rcol(rstart[M]);
    { std::vector<size_t> pos(rstart.begin(), rstart.end() - 1);
      for (size_t c = 0; c < C; ++c) for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) rcol[pos[celem[k]]++] = (uint32_t)c; }
    LOG("incidence " << cstart[C] << " nonzeros");
    for (size_t e = 0; e < M; ++e) if (rcnt[e] == 0) { LOG("uncoverable element"); return 1; }

    // ---- reductions
    std::vector<char> eAlive(M, 1), cAlive(C, 1);
    std::vector<uint32_t> ccnt(C), rc(M);     // live elements per column, live columns per element
    for (size_t c = 0; c < C; ++c) ccnt[c] = (uint32_t)(cstart[c + 1] - cstart[c]);
    for (size_t e = 0; e < M; ++e) rc[e] = rcnt[e];
    std::vector<int64_t> dom(C, -1);          // dominating column of a column removed by column dominance
    std::vector<uint32_t> forced;
    long nForced = 0, nColDom = 0, nRowDom = 0, nEmpty = 0;
    auto killElem = [&](size_t e) { eAlive[e] = 0; for (size_t k = rstart[e]; k < rstart[e + 1]; ++k) --ccnt[rcol[k]]; };
    auto killCol = [&](size_t c) { cAlive[c] = 0; for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) --rc[celem[k]]; };
    std::vector<uint32_t> stampE(M, 0), stampC(C, 0); uint32_t stamp = 0;
    for (int round = 0; ; ++round) {
        long before = nForced + nColDom + nRowDom + nEmpty;
        // forced columns
        for (size_t e = 0; e < M; ++e) if (eAlive[e] && rc[e] == 1) {
            size_t c = 0; for (size_t k = rstart[e]; k < rstart[e + 1]; ++k) if (cAlive[rcol[k]]) { c = rcol[k]; break; }
            forced.push_back((uint32_t)c); ++nForced;
            for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) if (eAlive[celem[k]]) killElem(celem[k]);
            killCol(c);
        }
        for (size_t c = 0; c < C; ++c) if (cAlive[c] && ccnt[c] == 0) { killCol(c); ++nEmpty; }
        // column dominance: c's live elements within a live column c2 (equal sets: keep the smaller index)
        for (size_t c = 0; c < C; ++c) if (cAlive[c]) {
            size_t best = SIZE_MAX; uint32_t bc = UINT32_MAX;
            for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) { uint32_t e = celem[k]; if (eAlive[e] && rc[e] < bc) { bc = rc[e]; best = e; } }
            ++stamp; for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) if (eAlive[celem[k]]) stampE[celem[k]] = stamp;
            for (size_t k = rstart[best]; k < rstart[best + 1]; ++k) {
                size_t c2 = rcol[k];
                if (c2 == c || !cAlive[c2] || ccnt[c2] < ccnt[c] || (ccnt[c2] == ccnt[c] && c2 > c)) continue;
                uint32_t hits = 0;
                for (size_t j = cstart[c2]; j < cstart[c2 + 1]; ++j) if (stampE[celem[j]] == stamp && eAlive[celem[j]]) ++hits;
                if (hits == ccnt[c]) { dom[c] = (int64_t)c2; killCol(c); ++nColDom; break; }
            }
        }
        // row dominance: if every live column of e covers f, then f is redundant (equal sets: keep the smaller index)
        for (size_t e = 0; e < M; ++e) if (eAlive[e]) {
            size_t cmin = SIZE_MAX; uint32_t best = UINT32_MAX;
            ++stamp;
            for (size_t k = rstart[e]; k < rstart[e + 1]; ++k) { uint32_t c = rcol[k]; if (!cAlive[c]) continue; stampC[c] = stamp; if (ccnt[c] < best) { best = ccnt[c]; cmin = c; } }
            for (size_t k = cstart[cmin]; k < cstart[cmin + 1]; ++k) {
                size_t f = celem[k];
                if (f == e || !eAlive[f] || rc[f] < rc[e] || (rc[f] == rc[e] && f < e)) continue;
                uint32_t hits = 0;
                for (size_t j = rstart[f]; j < rstart[f + 1]; ++j) if (stampC[rcol[j]] == stamp && cAlive[rcol[j]]) ++hits;
                if (hits == rc[e]) { killElem(f); ++nRowDom; }
            }
        }
        long after = nForced + nColDom + nRowDom + nEmpty;
        size_t le = 0, lc = 0; for (size_t e = 0; e < M; ++e) le += eAlive[e]; for (size_t c = 0; c < C; ++c) lc += cAlive[c];
        LOG("round " << round << ": forced " << nForced << " colDom " << nColDom << " rowDom " << nRowDom << " empty " << nEmpty
            << " | live elements " << le << " columns " << lc);
        if (after == before) break;
    }

    // ---- components of the live bipartite graph (union-find on elements through columns)
    std::vector<uint32_t> uf(M); std::iota(uf.begin(), uf.end(), 0);
    std::function<uint32_t(uint32_t)> find = [&](uint32_t a) { while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; } return a; };
    for (size_t c = 0; c < C; ++c) if (cAlive[c]) {
        int64_t first = -1;
        for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) { uint32_t e = celem[k]; if (!eAlive[e]) continue;
            if (first < 0) first = find(e); else { uint32_t r = find(e); if (r != first) uf[r] = (uint32_t)first; } }
    }
    std::unordered_map<uint32_t, int> compId; std::vector<std::vector<uint32_t>> compE;
    for (size_t e = 0; e < M; ++e) if (eAlive[e]) { uint32_t r = find(e); auto it = compId.find(r);
        if (it == compId.end()) { compId[r] = (int)compE.size(); compE.push_back({}); compE.back().push_back((uint32_t)e); } else compE[it->second].push_back((uint32_t)e); }
    std::vector<std::vector<uint32_t>> compC(compE.size());
    for (size_t c = 0; c < C; ++c) if (cAlive[c]) {
        for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) if (eAlive[celem[k]]) { compC[compId[find(celem[k])]].push_back((uint32_t)c); break; }
    }
    size_t K = compE.size(); size_t big = 0;
    for (size_t k = 0; k < K; ++k) big = std::max(big, compE[k].size());
    LOG("components " << K << ", largest " << big << " elements");

    // ---- start solution mapped to live columns: each start dumbbell -> an original column containing it -> follow
    // the dominance pointers to a live column (skipped if it reaches a forced column)
    std::vector<char> inStart(C, 0);
    std::vector<char> isForced(C, 0); for (auto c : forced) isForced[c] = 1;
    long unmapped = 0;
    for (auto& h : H) {
        Pair p = canonical(h); size_t e0 = EI(p.X[0], p.Y[0]);
        int64_t found = -1;
        for (size_t k = rstart[e0]; k < rstart[e0 + 1] && found < 0; ++k) {
            const Pair& q = cols[rcol[k]];
            auto sub = [](const std::vector<int>& a, const std::vector<int>& b) { return std::includes(b.begin(), b.end(), a.begin(), a.end()); };
            if ((sub(p.X, q.X) && sub(p.Y, q.Y)) || (sub(p.X, q.Y) && sub(p.Y, q.X))) found = rcol[k];
        }
        if (found < 0) { ++unmapped; continue; }
        while (dom[found] >= 0) found = dom[found];
        if (cAlive[found]) inStart[found] = 1;
    }
    if (unmapped) LOG("start dumbbells not contained in any column: " << unmapped << " (anchor cap?)");

    // ---- solve components, smallest first. Pass 1: every component with TL1 = 30 seconds.
    // Pass 2: the unsolved ones, smallest first, each with an equal share of the remaining time, warm-started.
    std::vector<size_t> order(K); std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return compE[a].size() < compE[b].size(); });
    const double TL1 = 30;
    struct CompRes { std::vector<char> pick; size_t size = 0; double lb = 1; bool opt = false; };
    std::vector<CompRes> res(K);
    std::vector<int64_t> eLocal(M, -1);
    auto localize = [&](size_t k, bool on) { auto& E = compE[k]; for (size_t i = 0; i < E.size(); ++i) eLocal[E[i]] = on ? (int64_t)i : -1; };
    auto local = [&](uint32_t e) { return eAlive[e] ? eLocal[e] : (int64_t)-1; };
    long startSum = 0, compTrivial = 0;
    // start: start columns of the component, completed greedily (most uncovered live elements first)
    auto startPick = [&](size_t k) {
        auto& E = compE[k]; auto& Cc = compC[k];
        std::vector<char> cov(E.size(), 0), pick(Cc.size(), 0); size_t uncovered = E.size();
        auto take = [&](size_t j) { pick[j] = 1;
            for (size_t t = cstart[Cc[j]]; t < cstart[Cc[j] + 1]; ++t) { int64_t l = local(celem[t]); if (l >= 0 && !cov[l]) { cov[l] = 1; --uncovered; } } };
        for (size_t j = 0; j < Cc.size(); ++j) if (inStart[Cc[j]]) take(j);
        while (uncovered) {
            size_t bj = 0, bg = 0;
            for (size_t j = 0; j < Cc.size(); ++j) if (!pick[j]) { size_t g = 0;
                for (size_t t = cstart[Cc[j]]; t < cstart[Cc[j] + 1]; ++t) { int64_t l = local(celem[t]); if (l >= 0 && !cov[l]) ++g; }
                if (g > bg) { bg = g; bj = j; } }
            take(bj);
        }
        return pick;
    };
    auto solve = [&](size_t k, double tl) {
        auto& E = compE[k]; auto& Cc = compC[k]; auto& R = res[k];
        IloEnv envC;
        try {
            IloModel model(envC); IloBoolVarArray x(envC, (IloInt)Cc.size());
            model.add(IloMinimize(envC, IloSum(x)));
            std::vector<std::vector<int>> rows(E.size());
            for (size_t j = 0; j < Cc.size(); ++j) for (size_t t = cstart[Cc[j]]; t < cstart[Cc[j] + 1]; ++t) {
                int64_t l = local(celem[t]); if (l >= 0) rows[l].push_back((int)j); }
            for (auto& r : rows) { IloExpr ex(envC); for (int j : r) ex += x[j]; model.add(ex >= 1); ex.end(); }
            IloCplex cp(model);
            cp.setOut(envC.getNullStream()); cp.setWarning(envC.getNullStream());
            cp.setParam(IloCplex::Param::TimeLimit, std::max(1.0, tl));
            cp.setParam(IloCplex::Param::Threads, threads);
            cp.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
            IloNumVarArray sv(envC); IloNumArray vals(envC);
            for (size_t j = 0; j < Cc.size(); ++j) { sv.add(x[j]); vals.add(R.pick[j] ? 1.0 : 0.0); }
            cp.addMIPStart(sv, vals, IloCplex::MIPStartCheckFeas);
            double t1 = T0();
            cp.solve();
            auto st = cp.getStatus();
            if (st == IloAlgorithm::Optimal || st == IloAlgorithm::Feasible) {
                std::vector<char> pk(Cc.size(), 0); size_t got = 0;
                for (size_t j = 0; j < Cc.size(); ++j) if (cp.getValue(x[(IloInt)j]) > 0.5) { pk[j] = 1; ++got; }
                if (got <= R.size) { R.pick = pk; R.size = got; }
                R.lb = std::max(R.lb, std::ceil(cp.getBestObjValue() - 1e-6));
                if (st == IloAlgorithm::Optimal) { R.opt = true; R.lb = (double)R.size; }
                if (R.lb >= (double)R.size) R.opt = true;
            }
            if (E.size() >= 300 || !R.opt)
                LOG("component " << E.size() << " elements, " << Cc.size() << " columns: " << R.size << ", lb " << R.lb << (R.opt ? " optimal" : "") << " (" << T0() - t1 << "s)");
        } catch (IloException& ex) { std::cerr << "CPLEX: " << ex << std::endl; envC.end(); exit(1); }
        envC.end();
    };
    for (size_t k : order) {
        localize(k, true);
        auto& R = res[k]; R.pick = startPick(k); R.size = std::count(R.pick.begin(), R.pick.end(), 1); startSum += (long)R.size;
        if (R.size == 1) { R.opt = true; ++compTrivial; }
        else solve(k, std::min(TL1, std::max(5.0, TL - T0())));   // small components get 5 s even after TL
        localize(k, false);
    }
    std::vector<size_t> hard; for (size_t k : order) if (!res[k].opt) hard.push_back(k);
    LOG("pass 1 done: " << hard.size() << " unsolved components");
    for (size_t i = 0; i < hard.size(); ++i) {
        double rem = TL - T0(); if (rem <= 1) break;
        localize(hard[i], true); solve(hard[i], rem / (double)(hard.size() - i)); localize(hard[i], false);
    }
    std::vector<uint32_t> chosen(forced.begin(), forced.end());
    double lbSum = (double)forced.size(); bool allOpt = !capHit;
    for (size_t k = 0; k < K; ++k) {
        for (size_t j = 0; j < compC[k].size(); ++j) if (res[k].pick[j]) chosen.push_back(compC[k][j]);
        lbSum += res[k].lb; if (!res[k].opt) allOpt = false;
    }
    std::vector<Pair> W; for (auto c : chosen) W.push_back(cols[c]);
    std::string why; bool ok = verify(I, W, &why);
    std::string gname = gen.rfind("file:", 0) == 0 ? gen.substr(gen.find_last_of('/') + 1) : gen;
    std::cout << gname << "," << n << "," << seed << "," << eps.val() << ",optx," << W.size() << "," << (ok ? "ok" : "FAIL:" + why) << "," << T0()
              << ",optimal=" << (allOpt && (double)W.size() == lbSum) << ";lb=" << (long)lbSum << ";start=" << H.size() << ";startMapped=" << startSum + (long)forced.size()
              << ";cols=" << C << ";forced=" << nForced << ";colDom=" << nColDom << ";rowDom=" << nRowDom << ";components=" << K << ";largest=" << big
              << ";trivial=" << compTrivial << ";capHit=" << capHit << std::endl;
    return 0;
}

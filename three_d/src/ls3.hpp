// Local search in R^3 (the planar ls.hpp), one move per run. Targets are taken in random-permutation passes over the dumbbells.
//   reshape: the target is replaced by the subdumbbell spanned by its private pairs, which then grows by free growth
//            (points join while the dumbbell stays eps-separated, largest number of newly covered private pairs of other
//            dumbbells first, ties at random), and the dumbbells that became redundant are deleted. The paper's heuristic.
//   lns    : the target and its k nearest neighbours are removed and the pairs left uncovered are re-covered optimally
//            (CPLEX set cover over the removed dumbbells and the maximal dumbbells anchored at those pairs); accepted if
//            smaller. The comparison. Neighbourhood, parameters, and the cached neighbour search as in the planar ls.hpp.
// Both moves keep the cover valid and never increase the size.
// cnt[e] = number of dumbbells covering e, own[e] = XOR of their ids (packed upper triangles, 4 bytes each); for LNS
// also point -> dumbbell incidence lists. The points are relabeled in Morton (Z-order) order first, for memory locality.
#pragma once
#include <fstream>
#include <random>
#include "wspd3.hpp"

namespace wspd3 {

struct LSConfig {
    double timeLimit = 10;
    std::string moves = "reshape";  // reshape | lns
    int lnsK = 16;                  // LNS: neighbours removed with the target
    double lnsTL = 2.0;             // LNS: CPLEX time limit per call
    size_t lnsCap = 300;            // LNS: maximal dumbbells enumerated per pair
    double radius = 2.0;            // LNS: neighbourhood radius in units of eps * scale
    size_t lnsUmax = 1000;          // LNS: candidate dumbbells are enumerated only if at most this many pairs are uncovered
    void lnsDefaults(int n) { lnsK = n >= 5000 ? 32 : 16; lnsUmax = n >= 5000 ? 400 : 1000; radius = 2.0; }   // as the planar ls.hpp
    unsigned seed = 1;
    std::string traceFile;          // if nonempty: "time,size" lines
    double traceEvery = 0.25;
    long maxRounds = 0;             // stop after this many rounds (0: no limit)
};

struct LSResult { long rounds = 0, reshapeMoves = 0, reshapeKills = 0, extraKills = 0, lnsCalls = 0, lnsGain = 0, lnsBig = 0, lnsFail = 0, startSize = 0, finalSize = 0; double time = 0, setup = 0; };

class LocalSearch {
public:
    const Instance& I; int n; ElemIndex EI; Spatial SP;
    std::vector<Pair> W; std::vector<char> alive; std::vector<uint32_t> cnt, own; std::vector<std::vector<int>> inc;
    bool useInc; long aliveCount = 0; std::mt19937_64 rng;
    std::vector<int> perm; size_t permPos = 0;
    std::vector<int> touched; bool track = false;
    std::unordered_map<size_t, std::vector<Pair>> anchorCache; size_t cacheSize = 0;

    LocalSearch(const Instance& I_, std::vector<Pair> W0, unsigned seed, bool needInc)
        : I(I_), n(I_.n()), EI(I_.n()), SP(I_), W(std::move(W0)), alive(W.size(), 1), cnt(EI.size(), 0), own(EI.size(), 0),
          inc(needInc ? I_.n() : 0), useInc(needInc), rng(seed) {
        for (size_t k = 0; k < W.size(); ++k) addIdx((int)k);
        aliveCount = (long)W.size();
    }
    void cover(int k, int s) {
        const Pair& zz = W[k];
        for (int u : zz.X) for (int v : zz.Y) {
            size_t e = EI(u, v);
            if (s > 0 && track && cnt[e] == 1) touched.push_back((int)own[e]);
            cnt[e] += s; own[e] ^= (uint32_t)k;
        }
    }
    void addIdx(int k) { cover(k, +1); if (useInc) { for (int u : W[k].X) inc[u].push_back(k); for (int u : W[k].Y) inc[u].push_back(k); } }
    void removeIdx(int k) {
        cover(k, -1);
        if (useInc) for (auto* S : {&W[k].X, &W[k].Y}) for (int u : *S) { auto& L = inc[u]; L.erase(std::remove(L.begin(), L.end(), k), L.end()); }
    }
    int newPair(const Pair& zz) { W.push_back(zz); alive.push_back(1); int k = (int)W.size() - 1; addIdx(k); ++aliveCount; return k; }
    void release(int k) { alive[k] = 0; --aliveCount; }
    void kill(int k) { removeIdx(k); release(k); }
    void replace(int k, const Pair& zz) { removeIdx(k); W[k] = zz; addIdx(k); }
    bool redundant(const Pair& zz) const { for (int u : zz.X) for (int v : zz.Y) if (cnt[EI(u, v)] == 1) return false; return true; }
    ll dmin2(const std::vector<int>& X, const std::vector<int>& Y) const { ll m = LLONG_MAX; for (int u : X) for (int v : Y) m = std::min(m, I.d2(u, v)); return m; }

    int pickTarget() {
        for (int pass = 0; pass < 2; ++pass) {
            while (permPos < perm.size()) { int k = perm[permPos++]; if (alive[k]) return k; }
            perm.clear(); for (size_t k = 0; k < W.size(); ++k) if (alive[k]) perm.push_back((int)k);
            std::shuffle(perm.begin(), perm.end(), rng); permPos = 0;
        }
        return -1;
    }
    Pair privSub(const Pair& Z) const {
        std::set<int> sx, sy;
        for (int u : Z.X) for (int v : Z.Y) if (cnt[EI(u, v)] == 1) { sx.insert(u); sy.insert(v); }
        return Pair{std::vector<int>(sx.begin(), sx.end()), std::vector<int>(sy.begin(), sy.end())};
    }
    // free growth (as in the planar ls.hpp): candidates within eps * d(Z) of a point of each side; a candidate that does not
    // fit never fits later (diameters only grow, the closest distance only shrinks)
    void freeGrow(Pair& Z) {
        ll D2 = dmin2(Z.X, Z.Y);
        struct C { int v; int side; long g; uint64_t r; ll ms, mc; };
        std::vector<C> cs; ll dS[2] = {0, 0}, dm = D2;
        for (int s = 0; s < 2; ++s) { auto& S = s == 0 ? Z.X : Z.Y; for (size_t i = 0; i < S.size(); ++i) for (size_t j = i + 1; j < S.size(); ++j) dS[s] = std::max(dS[s], I.d2(S[i], S[j])); }
        for (int side = 0; side < 2; ++side) {
            auto& S = side == 0 ? Z.X : Z.Y; auto& O = side == 0 ? Z.Y : Z.X;
            for (int v : SP.ball(S[0], D2)) {
                if (std::find(S.begin(), S.end(), v) != S.end() || std::find(O.begin(), O.end(), v) != O.end()) continue;
                C c{v, side, 0, rng(), 0, LLONG_MAX};
                for (int u : S) c.ms = std::max(c.ms, I.d2(u, v));
                for (int u : O) { c.mc = std::min(c.mc, I.d2(u, v)); c.g += cnt[EI(u, v)] == 1; }
                cs.push_back(c);
            }
        }
        auto fits = [&](const C& c) { return I.withinEps(std::max(std::max(dS[c.side], c.ms), dS[1 - c.side]), std::min(dm, c.mc)); };
        auto filter = [&]() { size_t k = 0; for (auto& w : cs) if (fits(w)) cs[k++] = w; cs.resize(k); };
        filter();
        bool changed = false;
        while (!cs.empty()) {
            size_t bi = 0;
            for (size_t i = 1; i < cs.size(); ++i) if (cs[i].g > cs[bi].g || (cs[i].g == cs[bi].g && cs[i].r > cs[bi].r)) bi = i;
            C c = cs[bi]; cs.erase(cs.begin() + bi);
            (c.side == 0 ? Z.X : Z.Y).push_back(c.v); changed = true;
            dS[c.side] = std::max(dS[c.side], c.ms); dm = std::min(dm, c.mc);
            for (auto& w : cs) {
                ll d = I.d2(c.v, w.v);
                if (w.side == c.side) w.ms = std::max(w.ms, d); else { w.mc = std::min(w.mc, d); w.g += cnt[EI(c.v, w.v)] == 1; }
            }
            filter();
        }
        if (changed) { std::sort(Z.X.begin(), Z.X.end()); std::sort(Z.Y.begin(), Z.Y.end()); }
    }
    bool reshapeRound(int t, LSResult& R) {
        Pair Z = privSub(W[t]);
        if (Z.X.empty()) { kill(t); ++R.extraKills; return true; }
        freeGrow(Z);
        if (Z.X == W[t].X && Z.Y == W[t].Y) return false;
        track = true; touched.clear(); replace(t, Z); track = false; ++R.reshapeMoves;
        std::sort(touched.begin(), touched.end()); touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
        long c = 0;
        for (int k : touched) if (alive[k] && redundant(W[k])) { kill(k); ++c; }
        R.reshapeKills += c;
        return c > 0;
    }
    // ------------------------------------------------------------ LNS (comparison; the dumbbells are never modified in
    // place, so the scale and the centroids of a dumbbell are computed once)
    struct Geo { double s = -1, a[3], b[3]; };
    std::vector<Geo> geo; std::vector<int> stamp; int stampNow = 0; std::vector<int> scan;
    const Geo& geoOf(int k) {
        if ((size_t)k >= geo.size()) geo.resize(W.size());
        Geo& g = geo[k];
        if (g.s < 0) {
            g.s = std::sqrt((double)dmin2(W[k].X, W[k].Y));
            for (int t = 0; t < 2; ++t) {
                auto& S = t == 0 ? W[k].X : W[k].Y; double* c = t == 0 ? g.a : g.b; c[0] = c[1] = c[2] = 0;
                for (int u : S) { c[0] += I.x[u]; c[1] += I.y[u]; c[2] += I.z[u]; }
                for (int q = 0; q < 3; ++q) c[q] /= S.size();
            }
        }
        return g;
    }
    static double dd(const double* p, const double* q) { return std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2])); }
    struct Cand { int k; double dist; };
    // neighbourhood of T: dumbbells with a point within radius * eps * d(T) of a point of T, scale within a factor
    // (1+2eps)^2, and side centroids within radius * eps * d(T) (in one orientation); nearest first
    std::vector<Cand> neighbours(int t, double radius) {
        const Pair& T = W[t];
        Geo gT = geoOf(t);
        double sT = gT.s;
        ll qD2 = (ll)std::ceil(sT * sT); long rnum = (long)std::lround(radius * 1000);
        if (stamp.size() < W.size()) stamp.resize(W.size(), 0);
        ++stampNow; scan.clear();
        for (auto* side : {&T.X, &T.Y}) for (int u : *side) for (int zz : SP.ball(u, qD2, rnum, 1000))
            for (int k : inc[zz]) if (k != t && alive[k] && stamp[k] != stampNow) { stamp[k] = stampNow; scan.push_back(k); }
        double e = I.eps.val(), lo = 1.0 / ((1 + 2 * e) * (1 + 2 * e)), hi = 1.0 / lo;
        std::vector<Cand> out;
        for (int k : scan) {
            const Geo& g = geoOf(k);
            if (g.s < lo * sT || g.s > hi * sT) continue;
            double d1 = std::max(dd(g.a, gT.a), dd(g.b, gT.b)), d2 = std::max(dd(g.a, gT.b), dd(g.b, gT.a));
            double dist = std::min(d1, d2) / (e * sT);
            if (dist <= radius) out.push_back({k, dist});
        }
        std::sort(out.begin(), out.end(), [](const Cand& a, const Cand& b) { return a.dist < b.dist || (a.dist == b.dist && a.k < b.k); });
        return out;
    }
    bool lnsRound(int t, const LSConfig& cf, LSResult& R) {
        auto cands = neighbours(t, cf.radius);
        std::vector<int> N = {t};
        for (auto& c : cands) { if ((int)N.size() > cf.lnsK) break; if (alive[c.k]) N.push_back(c.k); }
        for (int k : N) removeIdx(k);
        std::vector<std::pair<int,int>> U; std::unordered_set<size_t> Us;
        for (int k : N) for (int u : W[k].X) for (int v : W[k].Y) { size_t e = EI(u, v); if (cnt[e] == 0 && Us.insert(e).second) U.push_back({u, v}); }
        ++R.lnsCalls;
        if (U.empty()) { for (int k : N) release(k); R.lnsGain += (long)N.size(); return true; }
        std::unordered_set<Pair, PairHash, PairEq> seen; std::vector<Pair> cols;
        for (int k : N) { auto c = canonical(W[k]); if (seen.insert(c).second) cols.push_back(c); }
        if (U.size() > cf.lnsUmax) ++R.lnsBig;
        else {
            for (auto& e : U) {
                size_t key = (size_t)std::min(e.first, e.second) * n + std::max(e.first, e.second);
                auto it = anchorCache.find(key);
                if (it == anchorCache.end()) {
                    std::vector<Pair> tmp; anchoredMaximalPairs(I, SP, e.first, e.second, tmp, cf.lnsCap);
                    cacheSize += tmp.size();
                    if (cacheSize > 3000000) { anchorCache.clear(); cacheSize = tmp.size(); }
                    it = anchorCache.emplace(key, std::move(tmp)).first;
                }
                for (auto& pr : it->second) {
                    bool touches = false;
                    for (int u : pr.X) { for (int v : pr.Y) if (Us.count(EI(u, v))) { touches = true; break; } if (touches) break; }
                    if (!touches) continue;
                    auto c = canonical(pr); if (seen.insert(c).second) cols.push_back(c);
                }
            }
        }
        CoverResult cr;
        try { cr = setCoverCplex(n, U, cols, cf.lnsTL); }
        catch (...) { ++R.lnsFail; for (int k : N) addIdx(k); return false; }
        if (!cr.chosen.empty() && cr.chosen.size() < N.size()) {
            for (int k : N) release(k);
            for (int c : cr.chosen) newPair(cols[c]);
            R.lnsGain += (long)(N.size() - cr.chosen.size());
            return true;
        }
        for (int k : N) addIdx(k);   // restore
        return false;
    }

    std::vector<Pair> run(const LSConfig& cf, LSResult& R) {
        auto t0 = std::chrono::steady_clock::now();
        if (cf.moves != "reshape" && cf.moves != "lns") throw std::runtime_error("unknown local-search move " + cf.moves);
        R.startSize = aliveCount;
        std::ofstream tr; if (!cf.traceFile.empty()) tr.open(cf.traceFile);
        double lastTrace = 0; long lastLogged = aliveCount;
        if (tr) tr << "0," << aliveCount << "\n";
        while (true) {
            double el = seconds_since(t0);
            if (el >= cf.timeLimit || (cf.maxRounds > 0 && R.rounds >= cf.maxRounds)) break;
            if (tr && ((aliveCount != lastLogged && el >= lastTrace + 0.001) || el >= lastTrace + cf.traceEvery)) {
                tr << el << "," << aliveCount << "\n"; lastTrace = el; lastLogged = aliveCount;
            }
            int t = pickTarget(); if (t < 0) break;
            ++R.rounds;
            if (redundant(W[t])) { kill(t); ++R.extraKills; continue; }
            if (cf.moves == "reshape") reshapeRound(t, R); else lnsRound(t, cf, R);
        }
        R.time = seconds_since(t0);
        std::vector<Pair> out; for (size_t k = 0; k < W.size(); ++k) if (alive[k]) out.push_back(std::move(W[k]));
        R.finalSize = (long)out.size();
        if (tr) tr << R.time << "," << out.size() << "\n";
        return out;
    }
};

// Morton (Z-order) order of the points, 21 bits per coordinate
inline std::vector<int> mortonOrder(const Instance& I) {
    int n = I.n(); if (n == 0) return {};
    ll mn[3] = {*std::min_element(I.x.begin(), I.x.end()), *std::min_element(I.y.begin(), I.y.end()), *std::min_element(I.z.begin(), I.z.end())};
    ll span = 1;
    for (int i = 0; i < n; ++i) span = std::max({span, I.x[i] - mn[0] + 1, I.y[i] - mn[1] + 1, I.z[i] - mn[2] + 1});
    const ll ORD = 1 << 21;
    auto spread = [](uint64_t v) {   // 21 bits -> every third bit
        v &= 0x1fffff;
        v = (v | v << 32) & 0x1f00000000ffffull; v = (v | v << 16) & 0x1f0000ff0000ffull; v = (v | v << 8) & 0x100f00f00f00f00full;
        v = (v | v << 4) & 0x10c30c30c30c30c3ull; v = (v | v << 2) & 0x1249249249249249ull; return v; };
    auto q = [&](ll c, ll m) { return (uint64_t)((long double)(c - m) / span * (ORD - 1)); };
    std::vector<uint64_t> key(n);
    for (int i = 0; i < n; ++i) key[i] = spread(q(I.x[i], mn[0])) << 2 | spread(q(I.y[i], mn[1])) << 1 | spread(q(I.z[i], mn[2]));
    std::vector<int> ord(n); std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return key[a] < key[b]; });
    return ord;
}

// local search on I from W: relabels the points in Morton order, runs, maps the result back
inline std::vector<Pair> runLocalSearch(const Instance& I, std::vector<Pair> W, const LSConfig& cf, LSResult& R) {
    auto tSetup = std::chrono::steady_clock::now();
    auto ord = mortonOrder(I); int n = I.n();
    std::vector<int> nw(n); Instance I2; I2.eps = I.eps;
    for (int i = 0; i < n; ++i) { nw[ord[i]] = i; I2.x.push_back(I.x[ord[i]]); I2.y.push_back(I.y[ord[i]]); I2.z.push_back(I.z[ord[i]]); }
    auto mapP = [](Pair& p, const std::vector<int>& m) { for (int& u : p.X) u = m[u]; for (int& u : p.Y) u = m[u];
        std::sort(p.X.begin(), p.X.end()); std::sort(p.Y.begin(), p.Y.end()); };
    for (auto& p : W) mapP(p, nw);
    { LocalSearch L(I2, std::move(W), cf.seed, cf.moves == "lns"); R.setup = seconds_since(tSetup); W = L.run(cf, R); }
    for (auto& p : W) mapP(p, ord);
    return W;
}

} // namespace wspd3

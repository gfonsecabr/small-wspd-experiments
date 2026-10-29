// Local search for minimum-size cover WSPDs (post-processing of any valid cover WSPD). Two moves, applied to targets T
// taken in random-permutation passes over the dumbbells:
//   reshape: replace T by the subdumbbell spanned by its private pairs, then add points one at a time as long as the
//            dumbbell stays eps-separated (free growth: the closest cross distance may decrease), choosing the point
//            that newly covers the most pairs private to other dumbbells (ties at random); then delete the dumbbells
//            that became redundant (only those that just lost a private pair can be). This is the heuristic of the paper.
//   lns    : remove T and its K nearest neighbours; re-cover the pairs left uncovered optimally (CPLEX set cover over
//            maximal separated pairs anchored at those pairs, plus the removed dumbbells); accept if smaller. Comparison.
// Both moves keep the cover valid and never increase the size.
// State: the dumbbells W (with alive flags), for every pair e the number cnt[e] of dumbbells covering it and the XOR
// own[e] of their ids (so the only cover of a private pair is own[e]), and, for LNS only, point -> dumbbell incidence
// lists. cnt and own are packed upper triangles (ElemIndex: n(n-1)/2 entries of 4 bytes each). runLocalSearch()
// relabels the points in Hilbert-curve order first (memory locality of the tables).
#pragma once
#include <random>
#include <fstream>
#include <numeric>
#include <climits>
#include "wspd.hpp"

namespace wspd {

struct LSConfig {
    double timeLimit = 10;          // seconds
    bool lns = false;               // the move: reshape (false) or LNS (true)
    int lnsK = 16;                  // LNS: number of neighbours removed with the target (lnsDefaults: 32 for n >= 5000)
    double lnsTL = 2.0;             // LNS: CPLEX time limit per call (1 thread)
    size_t lnsCap = 300;            // LNS: cap on maximal pairs enumerated per anchor
    double radius = 2.0;            // LNS: neighbourhood radius in units of eps * scale
    size_t lnsUmax = 1000;          // LNS: new candidate dumbbells are generated only if at most this many pairs are left
                                    // uncovered (lnsDefaults: 400 for n >= 5000)
    // The parameters of LNS for n points: the best values of the tests of 2026-10-12 (radius 1.5 to 5, 4 to 64 neighbours,
    // caps 50 to 4000 and none, on the local-search samples with n = 1000, 3000, 10000, 50000): radius 2 for every n,
    // 16 neighbours and cap 1000 for n <= 3000, 32 neighbours and cap 400 for n >= 10000 (the threshold 5000 is arbitrary).
    void lnsDefaults(int n) { lnsK = n >= 5000 ? 32 : 16; lnsUmax = n >= 5000 ? 400 : 1000; radius = 2.0; }
    unsigned seed = 1;
    std::string traceFile;          // if nonempty: "time,size" lines
    double traceEvery = 0.25;
    long maxRounds = 0;             // stop after this many rounds (0: no limit; tests of determinism)
};

struct LSResult {
    long rounds = 0, reshapeMoves = 0, reshapeKills = 0, lnsCalls = 0, lnsGain = 0, extraKills = 0;
    long startSize = 0, finalSize = 0; double time = 0; double setup = 0, finish = 0;
};

class LocalSearch {
public:
    const Instance& I; int n; ElemIndex EI; Spatial SP;
    std::vector<Pair> W; std::vector<char> alive; std::vector<uint32_t> cnt, own; std::vector<std::vector<int>> inc;
    bool useInc; long aliveCount = 0; std::mt19937_64 rng;
    std::vector<int> perm; size_t permPos = 0;
    std::vector<int> touched; bool track = false;
    std::unordered_map<size_t, std::vector<Pair>> anchorCache; size_t cacheSize = 0;

    LocalSearch(const Instance& I_, std::vector<Pair> W0, unsigned seed, bool needInc)
        : I(I_), n(I_.n()), EI(I_.n()), SP(I_), W(std::move(W0)), alive(W.size(), 1), cnt(EI.size(), 0),
          own(EI.size(), 0), inc(needInc ? I_.n() : 0), useInc(needInc), rng(seed) {
        for (size_t k = 0; k < W.size(); ++k) addIdx((int)k);
        aliveCount = (long)W.size();
    }

    // ------------------------------------------------------------ bookkeeping
    // with track set, every dumbbell that loses a private pair (count 1 -> 2) is recorded in touched
    void cover(int k, int s) {
        const Pair& z = W[k];
        for (int u : z.X) for (int v : z.Y) {
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
    int newPair(const Pair& z) { W.push_back(z); alive.push_back(1); int k = (int)W.size() - 1; addIdx(k); ++aliveCount; return k; }
    // after removeIdx(k). The point vectors of W[k] are kept: freeing them made reshape ~10% slower (heap fragmentation)
    void release(int k) { alive[k] = 0; --aliveCount; }
    void kill(int k) { removeIdx(k); release(k); }
    void replace(int k, const Pair& z) { removeIdx(k); W[k] = z; addIdx(k); }
    bool redundant(const Pair& z) const { for (int u : z.X) for (int v : z.Y) if (cnt[EI(u, v)] == 1) return false; return true; }
    ll dmin2(const std::vector<int>& X, const std::vector<int>& Y) const { ll m = LLONG_MAX; for (int u : X) for (int v : Y) m = std::min(m, I.d2(u, v)); return m; }
    std::pair<double,double> centroid(const std::vector<int>& S) const { double cx = 0, cy = 0; for (int u : S) { cx += I.x[u]; cy += I.y[u]; } return {cx / S.size(), cy / S.size()}; }
    static double dd(std::pair<double,double> a, std::pair<double,double> b) { return std::hypot(a.first - b.first, a.second - b.second); }

    // targets: random-permutation passes over the alive dumbbells (dumbbells created during a pass wait for the next)
    int pickTarget() {
        for (int pass = 0; pass < 2; ++pass) {
            while (permPos < perm.size()) { int k = perm[permPos++]; if (alive[k]) return k; }
            perm.clear(); for (size_t k = 0; k < W.size(); ++k) if (alive[k]) perm.push_back((int)k);
            std::shuffle(perm.begin(), perm.end(), rng); permPos = 0;
        }
        return -1;
    }

    // ------------------------------------------------------------ move: reshape
    // subdumbbell spanned by the private pairs (empty sides if none)
    Pair privSub(const Pair& Z) const {
        std::set<int> sx, sy;
        for (int u : Z.X) for (int v : Z.Y) if (cnt[EI(u, v)] == 1) { sx.insert(u); sy.insert(v); }
        return Pair{std::vector<int>(sx.begin(), sx.end()), std::vector<int>(sy.begin(), sy.end())};
    }
    // free growth of Z: add points while Z stays eps-separated (exact test), largest gain first. Candidates are the
    // points within eps * d(Z) of a point of each side (d(Z) only decreases). Per candidate: ms = max squared distance
    // to its side, mc = min squared distance to the other side. Diameters only grow and the closest distance only
    // shrinks, so a candidate that does not fit never fits later.
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

    // ------------------------------------------------------------ move: LNS (comparison)
    // neighbourhood of T: dumbbells with a point within radius * eps * d(T) of a point of T, scale within a factor
    // (1+2eps)^2, and side centroids within radius * eps * d(T) (in one orientation); nearest first
    struct Cand { int k; double dist; };
    // LNS never modifies a dumbbell, so the scale and the side centroids of a dumbbell are computed once (geo); the
    // dumbbells that share a nearby point with T are collected without duplicates by stamps
    struct Geo { double s = -1, ax, ay, bx, by; };
    std::vector<Geo> geo; std::vector<int> stamp, scan; int stampNow = 0;
    const Geo& geoOf(int k) {
        if ((size_t)k >= geo.size()) geo.resize(W.size());
        Geo& g = geo[k];
        if (g.s < 0) {
            g.s = std::sqrt((double)dmin2(W[k].X, W[k].Y));
            auto a = centroid(W[k].X), b = centroid(W[k].Y); g.ax = a.first; g.ay = a.second; g.bx = b.first; g.by = b.second;
        }
        return g;
    }
    std::vector<Cand> neighbours(int t, double radius) {
        const Pair& T = W[t];
        Geo gT = geoOf(t);   // a copy: geoOf may reallocate
        double sT = gT.s;
        std::pair<double,double> cA = {gT.ax, gT.ay}, cB = {gT.bx, gT.by};
        ll qD2 = (ll)std::ceil(sT * sT); long rnum = (long)std::lround(radius * 1000);
        if (stamp.size() < W.size()) stamp.resize(W.size(), 0);
        ++stampNow; scan.clear();
        for (auto* side : {&T.X, &T.Y}) for (int u : *side) for (int z : SP.ball(u, qD2, rnum, 1000))
            for (int k : inc[z]) if (k != t && alive[k] && stamp[k] != stampNow) { stamp[k] = stampNow; scan.push_back(k); }
        double e = I.eps.val(), lo = 1.0 / ((1 + 2 * e) * (1 + 2 * e)), hi = 1.0 / lo;
        std::vector<Cand> out;
        for (int k : scan) {
            const Geo& g = geoOf(k);
            if (g.s < lo * sT || g.s > hi * sT) continue;
            std::pair<double,double> cX = {g.ax, g.ay}, cY = {g.bx, g.by};
            double d1 = std::max(dd(cX, cA), dd(cY, cB)), d2 = std::max(dd(cX, cB), dd(cY, cA));
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
        if (U.size() <= cf.lnsUmax) {
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
        try { cr = setCoverCplex(n, U, cols, cf.lnsTL, 1); }
        catch (...) { for (int k : N) addIdx(k); return false; }
        if (!cr.chosen.empty() && cr.chosen.size() < N.size()) {
            for (int k : N) release(k);
            for (int c : cr.chosen) newPair(cols[c]);
            R.lnsGain += (long)(N.size() - cr.chosen.size());
            return true;
        }
        for (int k : N) addIdx(k);       // restore
        return false;
    }

    // ------------------------------------------------------------ driver
    std::vector<Pair> run(const LSConfig& cf, LSResult& R) {
        auto t0 = std::chrono::steady_clock::now();
        R.startSize = aliveCount;
        std::ofstream tr; if (!cf.traceFile.empty()) tr.open(cf.traceFile);
        double lastTrace = 0; long lastLogged = aliveCount;   // log every size change (at most once per ms) and a heartbeat
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
            if (cf.lns) lnsRound(t, cf, R); else reshapeRound(t, R);
        }
        R.time = seconds_since(t0);
        std::vector<Pair> out; for (size_t k = 0; k < W.size(); ++k) if (alive[k]) out.push_back(std::move(W[k]));   // ends the search
        R.finalSize = (long)out.size();
        if (tr) tr << R.time << "," << out.size() << "\n";
        return out;
    }
};

// Hilbert-curve order of the points (a permutation of 0..n-1)
inline std::vector<int> hilbertOrder(const Instance& I) {
    int n = I.n(); if (n == 0) return {};
    ll mx = *std::min_element(I.x.begin(), I.x.end()), my = *std::min_element(I.y.begin(), I.y.end()), span = 1;
    for (int i = 0; i < n; ++i) span = std::max({span, I.x[i] - mx + 1, I.y[i] - my + 1});
    const ll ORD = 1 << 20;
    auto hil = [&](ll px, ll py) {
        ll x = (ll)((long double)(px - mx) / span * (ORD - 1)), y = (ll)((long double)(py - my) / span * (ORD - 1)), d = 0;
        for (ll s = ORD / 2; s > 0; s /= 2) {
            ll rx = (x & s) > 0, ry = (y & s) > 0; d += s * s * ((3 * rx) ^ ry);
            if (ry == 0) { if (rx == 1) { x = s - 1 - x; y = s - 1 - y; } std::swap(x, y); }
        }
        return d; };
    std::vector<ll> key(n); for (int i = 0; i < n; ++i) key[i] = hil(I.x[i], I.y[i]);
    std::vector<int> ord(n); std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return key[a] < key[b]; });
    return ord;
}

// Local search on I from W (moved in: only one copy of the WSPD is kept): relabels the points in Hilbert order, runs,
// maps the result back.
inline std::vector<Pair> runLocalSearch(const Instance& I, std::vector<Pair> W, const LSConfig& cf, LSResult& R) {
    bool needInc = cf.lns;
    auto tSetup = std::chrono::steady_clock::now();
    auto ord = hilbertOrder(I); int n = I.n();
    std::vector<int> nw(n); Instance I2; I2.eps = I.eps;
    for (int i = 0; i < n; ++i) { nw[ord[i]] = i; I2.x.push_back(I.x[ord[i]]); I2.y.push_back(I.y[ord[i]]); }
    auto mapP = [](Pair& p, const std::vector<int>& m) { for (int& u : p.X) u = m[u]; for (int& u : p.Y) u = m[u];
        std::sort(p.X.begin(), p.X.end()); std::sort(p.Y.begin(), p.Y.end()); };
    for (auto& p : W) mapP(p, nw);
    { LocalSearch L(I2, std::move(W), cf.seed, needInc); R.setup = seconds_since(tSetup); W = L.run(cf, R); }
    auto tFin = std::chrono::steady_clock::now();
    for (auto& p : W) mapP(p, ord);
    R.finish = seconds_since(tFin);
    return W;
}

} // namespace wspd

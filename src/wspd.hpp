// Core routines: instances, exact separation tests, the sweep, Callahan-Kosaraju, the maximal dumbbells and exact OPT
// by set cover (CPLEX), the redundancy pass, and the verifier (points in the plane, integer coordinates).
// Geometry: pgl (kd-tree ShapeTree with labelled points). ILPs: CPLEX (Concert).
// Separation parameter eps = a/b (rational).  All separation tests are exact:
//   max(diam X, diam Y) <= eps * d(X,Y)   <=>   b^2 * diam^2 <= a^2 * d^2.
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <unordered_set>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "pgl.hpp"

namespace wspd {

using ll = long long;
using i128 = __int128;
using LPoint = pgl::Point<ll, int>;          // label = point index
using Rect = pgl::Rectangle<pgl::Point<ll>>;

struct Eps { ll a, b; double val() const { return double(a) / double(b); } };

struct Pair { std::vector<int> X, Y; };      // sides sorted increasingly

struct Instance {
    std::vector<ll> x, y;
    Eps eps;
    int n() const { return (int)x.size(); }
    ll d2(int i, int j) const { ll dx = x[i] - x[j], dy = y[i] - y[j]; return dx * dx + dy * dy; }
    // |ij| <= eps * sqrt(D2)
    bool withinEps(ll dij2, ll D2) const { return (i128)eps.b * eps.b * dij2 <= (i128)eps.a * eps.a * D2; }
};

inline double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------- separation
inline bool separated(const Instance& I, const std::vector<int>& X, const std::vector<int>& Y) {
    if (X.empty() || Y.empty()) return false;
    ll diam2 = 0, dmin2 = LLONG_MAX;
    for (size_t i = 0; i < X.size(); ++i) for (size_t j = i + 1; j < X.size(); ++j) diam2 = std::max(diam2, I.d2(X[i], X[j]));
    for (size_t i = 0; i < Y.size(); ++i) for (size_t j = i + 1; j < Y.size(); ++j) diam2 = std::max(diam2, I.d2(Y[i], Y[j]));
    for (int u : X) for (int v : Y) { if (u == v) return false; dmin2 = std::min(dmin2, I.d2(u, v)); }
    return I.withinEps(diam2, dmin2);
}

// index of the element {i,j}, i != j, in a packed upper triangle: row i holds j = i+1..n-1 contiguously, from
// i(2n-i-1)/2 (an integer: i or 2n-i-1 is even). Tables indexed by it have size() = n(n-1)/2 entries.
struct ElemIndex {
    int n;
    explicit ElemIndex(int n_) : n(n_) {}
    size_t operator()(int i, int j) const { if (i > j) std::swap(i, j); return (size_t)i * (2 * (size_t)n - i - 1) / 2 + (j - i - 1); }
    size_t size() const { return (size_t)n * (n - 1) / 2; }
};

// ---------------------------------------------------------------- spatial index
struct Spatial {
    const Instance& I;
    pgl::ShapeTree<LPoint> tree;
    static std::vector<LPoint> pts(const Instance& I) {
        std::vector<LPoint> v; v.reserve(I.n());
        for (int i = 0; i < I.n(); ++i) v.emplace_back(I.x[i], I.y[i], i);
        return v;
    }
    explicit Spatial(const Instance& I_) : I(I_), tree(pts(I_)) {}
    // indices of points z with |z c| <= eps * sqrt(D2) * factor   (factor = num/den, exact)
    std::vector<int> ball(int c, ll D2, ll num = 1, ll den = 1) const { std::vector<int> out; ball(c, D2, out, num, den); return out; }
    void ball(int c, ll D2, std::vector<int>& out, ll num = 1, ll den = 1) const {
        // radius^2 = eps^2 * D2 * num^2/den^2 ; bounding box half-side = ceil(radius)
        long double r = std::sqrt((long double)D2) * I.eps.val() * num / den;
        ll h = (ll)std::ceil(r) + 1;
        Rect box(pgl::Point<ll>(I.x[c] - h, I.y[c] - h), pgl::Point<ll>(I.x[c] + h, I.y[c] + h));
        out.clear();
        tree.visitContainedIn(box, [&](const LPoint& p) {
            int z = p.label();
            // |zc|^2 * b^2 * den^2 <= a^2 * num^2 * D2
            if ((i128)I.eps.b * I.eps.b * den * den * I.d2(z, c) <= (i128)I.eps.a * I.eps.a * num * num * D2) out.push_back(z);
            return false;
        });
    }
    // ball(c, D2) computed from a superset `sup` of it that was itself returned by ball(): filtering keeps the
    // kd-tree visiting order, so the result equals ball(c, D2) element by element
    void filterBall(const std::vector<int>& sup, int c, ll D2, std::vector<int>& out) const {
        out.clear();
        for (int z : sup) if ((i128)I.eps.b * I.eps.b * I.d2(z, c) <= (i128)I.eps.a * I.eps.a * D2) out.push_back(z);
    }
};

// ---------------------------------------------------------------- maximal separated pairs with a given anchor
// All separated pairs whose closest cross pair is (p,q) are exactly the cliques containing p,q in the graph on
//   Xc = {x : |xp| <= eps|pq|, |xq| >= |pq|},  Yc = {y : |yq| <= eps|pq|, |yp| >= |pq|}
// with edges x-x' iff |xx'| <= eps|pq|, y-y' likewise, x-y iff |xy| >= |pq|.  Enumerate maximal cliques (Bron-Kerbosch).
inline void anchoredMaximalPairs(const Instance& I, const Spatial& S, int p, int q, std::vector<Pair>& out,
                                 size_t cap = 200000) {
    ll D2 = I.d2(p, q);
    std::vector<int> V = {p, q};
    std::vector<char> side = {0, 1};
    for (int x : S.ball(p, D2)) if (x != p && x != q && I.d2(x, q) >= D2) { V.push_back(x); side.push_back(0); }
    for (int y : S.ball(q, D2)) if (y != p && y != q && I.d2(y, p) >= D2) { V.push_back(y); side.push_back(1); }
    int m = (int)V.size();
    auto adj = [&](int i, int j) {
        if (i == j) return false;
        if (side[i] == side[j]) return I.withinEps(I.d2(V[i], V[j]), D2);
        return I.d2(V[i], V[j]) >= D2;
    };
    // candidates adjacent to both p (index 0) and q (index 1)
    std::vector<int> P0, X0;
    for (int i = 2; i < m; ++i) if (adj(i, 0) && adj(i, 1)) P0.push_back(i);
    std::vector<int> R = {0, 1};
    std::function<void(std::vector<int>&, std::vector<int>&, std::vector<int>&)> bk =
        [&](std::vector<int>& Rr, std::vector<int>& Pp, std::vector<int>& Xx) {
            if (out.size() >= cap) return;
            if (Pp.empty() && Xx.empty()) {
                Pair pr;
                for (int i : Rr) (side[i] == 0 ? pr.X : pr.Y).push_back(V[i]);
                std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
                out.push_back(std::move(pr));
                return;
            }
            // pivot maximizing |P ∩ N(u)|
            int piv = -1, best = -1;
            for (auto* L : {&Pp, &Xx}) for (int u : *L) {
                int c = 0; for (int v : Pp) if (adj(u, v)) ++c;
                if (c > best) { best = c; piv = u; }
            }
            std::vector<int> cand;
            for (int v : Pp) if (piv < 0 || !adj(piv, v)) cand.push_back(v);
            for (int v : cand) {
                std::vector<int> R2 = Rr; R2.push_back(v);
                std::vector<int> P2, X2;
                for (int w : Pp) if (w != v && adj(v, w)) P2.push_back(w);
                for (int w : Xx) if (adj(v, w)) X2.push_back(w);
                bk(R2, P2, X2);
                Pp.erase(std::find(Pp.begin(), Pp.end(), v));
                Xx.push_back(v);
            }
        };
    bk(R, P0, X0);
}

struct PairHash {
    size_t operator()(const Pair& p) const {
        size_t h = 1469598103934665603ull;
        for (int v : p.X) h = (h ^ (size_t)v) * 1099511628211ull;
        h = (h ^ 0x9e3779b97f4a7c15ull) * 1099511628211ull;
        for (int v : p.Y) h = (h ^ (size_t)v) * 1099511628211ull;
        return h;
    }
};
struct PairEq { bool operator()(const Pair& a, const Pair& b) const { return a.X == b.X && a.Y == b.Y; } };

inline Pair canonical(Pair p) { if (p.Y < p.X) std::swap(p.X, p.Y); return p; }

// ---------------------------------------------------------------- set cover by CPLEX
// elements: list of (i,j); columns: pairs. Returns indices of chosen columns. timeLimit in seconds.
struct CoverResult { std::vector<int> chosen; double lb = 0; bool optimal = false; };

inline CoverResult setCoverCplex(int n, const std::vector<std::pair<int,int>>& elems, const std::vector<Pair>& cols,
                                 double timeLimit = 60, int threads = 1) {   // always call with threads = 1 (policy)
    CoverResult res;
    ElemIndex EI(n);
    std::unordered_map<size_t, int> eid; eid.reserve(elems.size() * 2);
    for (size_t k = 0; k < elems.size(); ++k) eid[EI(elems[k].first, elems[k].second)] = (int)k;
    std::vector<std::vector<int>> rows(elems.size());
    for (size_t c = 0; c < cols.size(); ++c)
        for (int u : cols[c].X) for (int v : cols[c].Y) {
            auto it = eid.find(EI(u, v));
            if (it != eid.end()) rows[it->second].push_back((int)c);
        }
    // exact brute force for tiny instances: optimum 1 or 2 with <= 64 elements (bitmasks)
    if (elems.size() <= 64) {
        std::vector<uint64_t> mask(cols.size(), 0); uint64_t full = elems.size() == 64 ? ~0ull : ((1ull << elems.size()) - 1);
        for (size_t r = 0; r < rows.size(); ++r) for (int c : rows[r]) mask[c] |= (1ull << r);
        for (size_t c = 0; c < cols.size(); ++c) if (mask[c] == full) { res.chosen = {(int)c}; res.lb = 1; res.optimal = true; return res; }
        if (cols.size() <= 400)
            for (size_t c = 0; c < cols.size(); ++c) for (size_t c2 = c + 1; c2 < cols.size(); ++c2)
                if ((mask[c] | mask[c2]) == full) { res.chosen = {(int)c, (int)c2}; res.lb = 2; res.optimal = true; return res; }
    }
    IloEnv env;   // fresh environment per call (a shared one accumulates constraint objects and slows down)
    try {
        IloModel model(env);
        IloBoolVarArray x(env, (IloInt)cols.size());
        model.add(IloMinimize(env, IloSum(x)));
        for (auto& r : rows) {
            if (r.empty()) throw std::runtime_error("uncoverable element");
            IloExpr e(env); for (int c : r) e += x[c]; model.add(e >= 1); e.end();
        }
        IloCplex cp(model);
        cp.setOut(env.getNullStream()); cp.setWarning(env.getNullStream());
        cp.setParam(IloCplex::Param::TimeLimit, timeLimit);
        cp.setParam(IloCplex::Param::Threads, threads);
        // the default relative gap (1e-4) lets CPLEX report "optimal" up to 0.01% above the optimum for n >= 400
        cp.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
        cp.solve();
        for (size_t c = 0; c < cols.size(); ++c) if (cp.getValue(x[(IloInt)c]) > 0.5) res.chosen.push_back((int)c);
        res.lb = cp.getBestObjValue();
        res.optimal = (cp.getStatus() == IloAlgorithm::Optimal);
    } catch (IloException& ex) { std::cerr << "CPLEX: " << ex << std::endl; env.end(); throw; }
    env.end();
    return res;
}

// ---------------------------------------------------------------- largest set of pairwise incompatible pairs
// Two pairs are compatible iff some eps-separated dumbbell covers both, iff some maximal one does. So a set of pairs is
// pairwise incompatible iff every column of cols (all maximal dumbbells) covers at most one of its pairs, and its size
// is a lower bound on OPT. Integer program: one 0/1 variable per pair, maximize their sum, at most one chosen pair per
// column. Its LP relaxation is the dual of the LP relaxation of the set cover of opt.
struct PackResult { std::vector<std::pair<int,int>> chosen; double ub = 0; bool optimal = false; };

inline PackResult maxIncompatibleCplex(int n, const std::vector<Pair>& cols, double timeLimit = 1800, int threads = 1) {
    PackResult res;
    ElemIndex EI(n);
    IloEnv env;
    try {
        IloModel model(env);
        IloBoolVarArray y(env, (IloInt)EI.size());
        model.add(IloMaximize(env, IloSum(y)));
        for (auto& c : cols) {
            if (c.X.size() * c.Y.size() < 2) continue;
            IloExpr e(env); for (int u : c.X) for (int v : c.Y) e += y[(IloInt)EI(u, v)];
            model.add(e <= 1); e.end();
        }
        IloCplex cp(model);
        cp.setOut(env.getNullStream()); cp.setWarning(env.getNullStream());
        cp.setParam(IloCplex::Param::TimeLimit, timeLimit);
        cp.setParam(IloCplex::Param::Threads, threads);
        cp.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
        cp.solve();
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
            if (cp.getValue(y[(IloInt)EI(i, j)]) > 0.5) res.chosen.push_back({i, j});
        res.ub = cp.getBestObjValue();
        res.optimal = (cp.getStatus() == IloAlgorithm::Optimal);
    } catch (IloException& ex) { std::cerr << "CPLEX: " << ex << std::endl; env.end(); throw; }
    env.end();
    return res;
}

// all maximal separated pairs (small n): union over all anchors, deduplicated
inline std::vector<Pair> allMaximalPairs(const Instance& I, const Spatial& S) {
    std::unordered_set<Pair, PairHash, PairEq> seen;
    std::vector<Pair> out, tmp;
    for (int p = 0; p < I.n(); ++p) for (int q = 0; q < I.n(); ++q) if (p != q) {
        tmp.clear(); anchoredMaximalPairs(I, S, p, q, tmp);
        for (auto& pr : tmp) { auto c = canonical(pr); if (seen.insert(c).second) out.push_back(c); }
    }
    return out;
}

// ---------------------------------------------------------------- verification
inline bool verify(const Instance& I, const std::vector<Pair>& W, std::string* why = nullptr) {
    int n = I.n(); ElemIndex EI(n);
    std::vector<bool> cov(EI.size(), false);   // one bit per element
    for (auto& p : W) {
        if (!separated(I, p.X, p.Y)) { if (why) *why = "pair not separated"; return false; }
        for (int u : p.X) for (int v : p.Y) cov[EI(u, v)] = true;
    }
    for (size_t e = 0; e < cov.size(); ++e)
        if (!cov[e]) { if (why) *why = "element uncovered"; return false; }
    return true;
}

// ---------------------------------------------------------------- polynomial candidate generation (greedy growth)
// Grow ONE separated pair whose closest cross pair is the anchor (a,b): start from {a},{b} and repeatedly add the
// vertex (of the anchor graph used in anchoredMaximalPairs) that is compatible with everything chosen so far and
// maximizes (#elements of R newly covered, #currently uncovered elements newly covered), until no vertex fits.
// Time O(m^2 * |pair|) for m candidate vertices; no enumeration.
// Ba, Bb (optional): precomputed S.ball(a, |ab|^2) and S.ball(b, |ab|^2).
template <class IsR, class IsUncov>
inline Pair greedyAnchoredPair(const Instance& I, const Spatial& S, int a, int b, IsR inR, IsUncov uncov,
                               const std::vector<int>* Ba = nullptr, const std::vector<int>* Bb = nullptr) {
    ll D2 = I.d2(a, b);
    // same-side test |vw| <= eps|ab| as an exact integer threshold: b^2 d <= a^2 D2  <=>  d <= floor(a^2 D2 / b^2)
    i128 lim128 = (i128)I.eps.a * I.eps.a * D2 / ((i128)I.eps.b * I.eps.b);
    ll lim = lim128 > (i128)LLONG_MAX ? LLONG_MAX : (ll)lim128;
    // scratch buffers, reused across calls (this function is called millions of times on small inputs)
    static thread_local struct { std::vector<int> V, C, cs[2], ballA, ballB; std::vector<char> side; std::vector<ll> vx, vy; std::vector<int64_t> key; } sc;
    auto& V = sc.V; auto& side = sc.side; auto& ballA = sc.ballA; auto& ballB = sc.ballB;
    V.assign({a, b}); side.assign({0, 1});
    if (!Ba) { S.ball(a, D2, ballA); Ba = &ballA; }
    if (!Bb) { S.ball(b, D2, ballB); Bb = &ballB; }
    for (int x : *Ba) if (x != a && x != b && I.d2(x, b) >= D2) { V.push_back(x); side.push_back(0); }
    for (int y : *Bb) if (y != a && y != b && I.d2(y, a) >= D2) { V.push_back(y); side.push_back(1); }
    int m = (int)V.size();
    auto& vx = sc.vx; auto& vy = sc.vy; vx.resize(m); vy.resize(m);   // local contiguous copies of the coordinates
    for (int i = 0; i < m; ++i) { vx[i] = I.x[V[i]]; vy[i] = I.y[V[i]]; }
    auto d2 = [&](int i, int j) { ll dx = vx[i] - vx[j], dy = vy[i] - vy[j]; return dx * dx + dy * dy; };
    auto adj = [&](int i, int j) { return side[i] == side[j] ? d2(i, j) <= lim : d2(i, j) >= D2; };
    // candidates, one list per side, each in increasing index order
    auto& C = sc.C; auto& cs = sc.cs; C.assign({0, 1}); cs[0].clear(); cs[1].clear();
    for (int i = 2; i < m; ++i) if (adj(i, 0) && adj(i, 1)) cs[(int)side[i]].push_back(i);
    // incremental gains gR[v], gU[v] = elements (v, c) with c in C on the other side that are in R / uncovered, packed
    // with the tie-break (smallest index first) as key[v] = gR[v] * 2^42 + gU[v] * 2^21 + (2^21 - 1 - v), so that the
    // choice is the unique maximum key, independent of the scan order
    constexpr int SH = 21;
    if (m >= (1 << SH)) throw std::runtime_error("greedyAnchoredPair: too many candidates");
    auto& key = sc.key; key.resize(m);
    for (int v = 0; v < m; ++v) key[v] = ((int64_t)1 << SH) - 1 - v;
    // inR and uncov must be symmetric; they are called with the vertex just added first, so that consecutive
    // calls share their first argument (memory locality for row-major tables)
    auto addGain = [&](int v, int c) { key[v] += ((int64_t)inR(V[c], V[v]) << (2 * SH)) + ((int64_t)uncov(V[c], V[v]) << SH); };
    int u = -1; int64_t best = -1;
    for (int s = 0; s < 2; ++s) for (int v : cs[s]) { addGain(v, 1 - s); if (key[v] > best) { best = key[v]; u = v; } }
    while (u >= 0) {
        C.push_back(u);
        int su = side[u], nu = -1; best = -1;
        // same side as u: keep the candidates within eps|ab| of u (in-place compaction, order kept)
        auto& L = cs[su]; size_t k = 0;
        for (int w : L) if (w != u && d2(u, w) <= lim) { L[k++] = w; if (key[w] > best) { best = key[w]; nu = w; } }
        L.resize(k);
        // other side: keep the candidates at distance >= |ab| from u, and update their gains
        auto& M = cs[1 - su]; k = 0;
        for (int w : M) if (d2(u, w) >= D2) { addGain(w, u); M[k++] = w; if (key[w] > best) { best = key[w]; nu = w; } }
        M.resize(k);
        u = nu;
    }
    // exact capacities: the pairs are moved into the output, which can hold millions of them
    Pair pr; size_t nx = 0; for (int i : C) nx += side[i] == 0;
    pr.X.reserve(nx); pr.Y.reserve(C.size() - nx);
    for (int i : C) (side[i] == 0 ? pr.X : pr.Y).push_back(V[i]);
    std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
    return pr;
}

// greedy set cover of elements R (pairs (x,y), x on the X side) by columns; returns chosen column indices
inline std::vector<int> greedySetCover(int n, const std::vector<std::pair<int,int>>& R, const std::vector<Pair>& cols) {
    if (R.size() <= 64) {   // fast path, same choices: one bitmask of R per column (sides are sorted: binary search)
        std::vector<uint64_t> mask(cols.size(), 0);
        auto in = [](const std::vector<int>& A, int v) { return std::binary_search(A.begin(), A.end(), v); };
        for (size_t c = 0; c < cols.size(); ++c) for (size_t k = 0; k < R.size(); ++k) {
            auto [x, y] = R[k];
            if ((in(cols[c].X, x) && in(cols[c].Y, y)) || (in(cols[c].X, y) && in(cols[c].Y, x))) mask[c] |= 1ull << k;
        }
        uint64_t left = R.size() == 64 ? ~0ull : (1ull << R.size()) - 1; std::vector<int> chosen;
        while (left) {
            int best = -1; int bg = 0;
            for (size_t c = 0; c < cols.size(); ++c) { int g = std::popcount(mask[c] & left); if (g > bg) { bg = g; best = (int)c; } }
            if (best < 0) throw std::runtime_error("greedySetCover: uncoverable element");
            chosen.push_back(best); left &= ~mask[best];
        }
        return chosen;
    }
    ElemIndex EI(n);
    std::unordered_map<size_t, int> eid; for (size_t k = 0; k < R.size(); ++k) eid[EI(R[k].first, R[k].second)] = (int)k;
    std::vector<std::vector<int>> colElems(cols.size());
    for (size_t c = 0; c < cols.size(); ++c)
        for (int u : cols[c].X) for (int v : cols[c].Y) { auto it = eid.find(EI(u, v)); if (it != eid.end()) colElems[c].push_back(it->second); }
    std::vector<char> done(R.size(), 0); size_t left = R.size(); std::vector<int> chosen;
    while (left > 0) {
        int best = -1; long bg = 0;
        for (size_t c = 0; c < cols.size(); ++c) { long g = 0; for (int e : colElems[c]) g += !done[e]; if (g > bg) { bg = g; best = (int)c; } }
        if (best < 0) throw std::runtime_error("greedySetCover: uncoverable element");
        chosen.push_back(best); for (int e : colElems[best]) if (!done[e]) { done[e] = 1; --left; }
    }
    return chosen;
}

// ---------------------------------------------------------------- lexicographic order
inline std::vector<int> lexOrder(const Instance& I) {
    std::vector<int> ord(I.n()); std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return I.x[a] != I.x[b] ? I.x[a] < I.x[b] : I.y[a] < I.y[b]; });
    return ord;
}

// ---------------------------------------------------------------- pieces of the fixed collection (Theorem sweep)
// piece index of point z relative to centre c with radius^2 = eps^2 D2: disk -> 8 pieces, right half-disk -> 5 pieces
inline int diskPiece(const Instance& I, int c, int z, ll D2) {
    // centre piece: |zc| <= 0.434 * eps * |pq|
    long double r = std::sqrt((long double)I.d2(c, z)), R = std::sqrt((long double)D2) * I.eps.val();
    if (r <= 0.434L * R) return 0;
    long double ang = std::atan2((long double)(I.y[z] - I.y[c]), (long double)(I.x[z] - I.x[c]));
    if (ang < 0) ang += 2 * M_PIl;
    int k = (int)(ang / (2 * M_PIl / 7)); if (k > 6) k = 6;
    return 1 + k;
}
inline int halfDiskPiece(const Instance& I, int c, int z, ll D2) {
    long double r = std::sqrt((long double)I.d2(c, z)), R = std::sqrt((long double)D2) * I.eps.val();
    if (r <= 0.4L * R) return 0;
    long double ang = std::atan2((long double)(I.y[z] - I.y[c]), (long double)(I.x[z] - I.x[c])); // in [-pi/2, pi/2]
    int k = (int)((ang + M_PIl / 2) / (M_PIl / 4)); k = std::clamp(k, 0, 3);
    return 1 + k;
}

struct SweepStats { long capHits = 0, steps = 0, fixedChosen = 0, candTotal = 0, rTotal = 0; double time = 0, tCand = 0, tCover = 0; };

// The sweep. Points p in lexicographic order; for each p, the points q after p by increasing distance; every pair {p,q}
// not yet covered is the anchor of a step, and the number of steps is a lower bound on the minimum size: the anchors
// are pairwise incompatible, since a step covers every uncovered pair after p that is compatible with its anchor (the
// set R_pq). If anchors is given, the anchors are appended to it, in the order of the steps (check_lb.cpp checks them).
// This is the greedy sweep: each step covers the must-cover set R_pq by grown dumbbells (greedy growth from (p,q), then
// from the closest still-uncovered pair of R_pq, repeatedly, if |R_pq| > LAZY_THRESHOLD, else a greedy set cover of the dumbbells grown from
// every pair of R_pq), and outputs the fixed collection of Theorem sweep (the 5 x 8 piece dumbbells) instead when it is
// smaller (eps <= 1/16 only). After GROWTH_CAP grown dumbbells, the fixed collection is used directly.
inline std::vector<Pair> sweep(const Instance& I, SweepStats& st, std::vector<std::pair<int,int>>* anchors = nullptr) {
    constexpr size_t LAZY_THRESHOLD = 16, GROWTH_CAP = 40;
    auto t0 = std::chrono::steady_clock::now();
    int n = I.n();
    Spatial S(I);
    auto ord = lexOrder(I);
    std::vector<int> rank(n); for (int k = 0; k < n; ++k) rank[ord[k]] = k;
    // one mark per element {u,v}: COV = covered by W; otherwise INR = in R and not yet covered in this step, DONE = in R
    // and covered by a pair grown in this step, 0 = neither. At the end of a step every element of R is covered (COV),
    // so no mark of an earlier step survives; the sweep aborts if one is not, since the number of steps would then not
    // be a lower bound. One byte per element (n^2 bytes: 10 GB for n = 100000). A single array keeps the lookups to one cache miss, and it is stored symmetrically (at u*n+v and
    // v*n+u) so that the lookups of greedyAnchoredPair, which share one endpoint, stay within one row.
    constexpr uint8_t INR = 1, DONE = 2, COV = 3;
    std::vector<uint8_t> mark((size_t)n * n, 0);
    auto MI = [n](int u, int v) { return (size_t)u * n + v; };
    auto setMark = [&](int u, int v, uint8_t m) { mark[MI(u, v)] = m; mark[MI(v, u)] = m; };
    auto isCov = [&](int u, int v) { return mark[MI(u, v)] == COV; };
    std::vector<Pair> W;
    auto coverPair = [&](const Pair& pr) { for (int u : pr.X) for (int v : pr.Y) setMark(u, v, COV); };
    bool fixedValid = (I.eps.a * 16 <= I.eps.b);
    std::vector<int> Bp, Bq, supP, supQ, ballA, ballB;   // per-step scratch
    std::vector<ll> distP(n);
    std::vector<std::pair<int,int>> R;
    std::array<std::vector<int>, 5> Apieces; std::array<std::vector<int>, 8> Cpieces;
    for (int k = 0; k < n; ++k) {
        int p = ord[k];
        std::vector<int> Q; for (int kk = k + 1; kk < n; ++kk) Q.push_back(ord[kk]);
        for (int z : Q) distP[z] = I.d2(p, z);
        std::sort(Q.begin(), Q.end(), [&](int a, int b) { return distP[a] < distP[b]; });
        for (int q : Q) {
            if (isCov(p, q)) continue;
            ++st.steps;
            if (anchors) anchors->push_back({p, q});
            ll D2 = I.d2(p, q);
            S.ball(p, D2, Bp); S.ball(q, D2, Bq);   // B(p, eps D), B(q, eps D), shared below
            // fixed collection: only the pieces and the collection size are computed here; the dumbbells are built by
            // makeFixed() when the collection is actually used
            size_t fixedSize = 0;
            if (fixedValid) {
                for (auto& a : Apieces) a.clear();
                for (auto& c : Cpieces) c.clear();
                for (int z : Bp) if (rank[z] >= rank[p]) Apieces[halfDiskPiece(I, p, z, D2)].push_back(z);
                for (int z : Bq) Cpieces[diskPiece(I, q, z, D2)].push_back(z);
                size_t na = 0, nc = 0;
                for (auto& a : Apieces) na += !a.empty();
                for (auto& c : Cpieces) nc += !c.empty();
                fixedSize = na * nc;
            }
            auto makeFixed = [&]() {
                std::vector<Pair> fam;
                for (auto& a : Apieces) for (auto& c : Cpieces) if (!a.empty() && !c.empty()) {
                    Pair pr{a, c}; std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
                    fam.push_back(pr);
                }
                return fam;
            };
            // must-cover set R: uncovered {x,y}, x in B(p,eps D), y in B(q, eps D), both >= p, compatible with {p,q}
            R.clear();
            for (int x : Bp) if (rank[x] >= rank[p]) for (int y : Bq) if (rank[y] >= rank[p] && x != y && !isCov(x, y)) {
                // separated(I, {p,x}, {q,y}), inlined without allocations
                if (y == p || x == q) continue;
                ll diam2 = std::max(I.d2(p, x), I.d2(q, y));
                ll dmin2 = std::min(std::min(D2, I.d2(p, y)), std::min(I.d2(x, q), I.d2(x, y)));
                if (I.withinEps(diam2, dmin2)) R.push_back({x, y});
            }
            st.rTotal += R.size(); auto tB = std::chrono::steady_clock::now();
            for (auto& e : R) setMark(e.first, e.second, INR);
            auto inR = [&](int x, int y) { uint8_t m = mark[MI(x, y)]; return (long)(m == INR || m == DONE); };
            auto uncov = [&](int x, int y) { return (long)(mark[MI(x, y)] != COV); };
            // growth from an anchor (x,y) of R or (p,q). For (x,y) in R, |xy| <= (1+2eps)|pq|, so the balls
            // B(x, eps|xy|) and B(y, eps|xy|) lie in B(p, (2+2eps) eps|pq|) and B(q, (2+2eps) eps|pq|): these
            // two are queried once per step and filtered for each anchor
            bool haveSup = false;
            auto growFrom = [&](int a, int b, auto inRfn) {
                if (a == p && b == q) return greedyAnchoredPair(I, S, a, b, inRfn, uncov, &Bp, &Bq);
                if (!haveSup) {
                    S.ball(p, D2, supP, 2 * (I.eps.a + I.eps.b), I.eps.b);
                    S.ball(q, D2, supQ, 2 * (I.eps.a + I.eps.b), I.eps.b);
                    haveSup = true;
                }
                ll Dab = I.d2(a, b);
                S.filterBall(supP, a, Dab, ballA); S.filterBall(supQ, b, Dab, ballB);
                return greedyAnchoredPair(I, S, a, b, inRfn, uncov, &ballA, &ballB);
            };
            std::unordered_set<Pair, PairHash, PairEq> seen; std::vector<Pair> cols;
            std::vector<int> ch;
            bool useLazy = R.size() > LAZY_THRESHOLD;
            if (useLazy) {
                // sequential cover: grow from (p,q), then from the closest still-uncovered element of R
                size_t left = R.size();
                auto isLeft = [&](int x, int y) { return mark[MI(x, y)] == INR; };
                auto inLeft = [&](int x, int y) { return (long)isLeft(x, y); };
                std::pair<int,int> anc = {p, q}; size_t idx = 0;
                while (true) {
                    auto pr = growFrom(anc.first, anc.second, inLeft);
                    for (int u : pr.X) for (int v : pr.Y) if (isLeft(u, v)) { setMark(u, v, DONE); --left; }
                    ch.push_back((int)cols.size()); cols.push_back(std::move(pr));
                    if (left == 0) break;
                    if (fixedSize > 0 && cols.size() >= GROWTH_CAP) {
                        // more than GROWTH_CAP >= |fixed collection| grown dumbbells would lose against the fixed
                        // collection anyway, so stop and use it (same output, less time)
                        ++st.capHits;
                        ch.clear(); cols = makeFixed(); for (size_t c = 0; c < cols.size(); ++c) ch.push_back((int)c);
                        break;
                    }
                    // next seed: the closest still-uncovered element of R (the first one in R among equals)
                    while (!isLeft(R[idx].first, R[idx].second)) ++idx;
                    ll best = LLONG_MAX;
                    for (size_t i = idx; i < R.size(); ++i) if (isLeft(R[i].first, R[i].second)) {
                        ll d = I.d2(R[i].first, R[i].second);
                        if (d < best) { best = d; anc = R[i]; }
                    }
                }
                st.tCand += seconds_since(tB); st.candTotal += cols.size();
            }
            auto tC = std::chrono::steady_clock::now();
            if (!useLazy) {
                std::vector<std::pair<int,int>> seeds = R; seeds.push_back({p, q});
                for (auto [a, b] : seeds) {
                    auto pr = growFrom(a, b, inR);
                    if (seen.insert(pr).second) cols.push_back(pr);
                }
                st.tCand += seconds_since(tB); st.candTotal += cols.size();
                tC = std::chrono::steady_clock::now();
                ch = greedySetCover(n, R, cols);
            }
            st.tCover += seconds_since(tC);
            std::vector<Pair> chosen;
            if (fixedSize == 0 || ch.size() <= fixedSize) { for (int c : ch) chosen.push_back(std::move(cols[c])); }
            else { ++st.fixedChosen; chosen = makeFixed(); }
            for (auto& pr : chosen) { coverPair(pr); W.push_back(std::move(pr)); }
            for (auto& e : R) if (!isCov(e.first, e.second)) {
                std::cerr << "sweep: the step of anchor {" << p << "," << q << "} leaves {" << e.first << "," << e.second
                          << "} of its must-cover set uncovered" << std::endl;
                std::abort();
            }
        }
    }
    st.time = seconds_since(t0);
    return W;
}

// ---------------------------------------------------------------- Callahan-Kosaraju WSPD (fair split tree)
struct CKNode { ll x0, y0, x1, y1; int left = -1, right = -1; std::vector<int> pts; };

// exact = false: separation tested on bounding boxes (box diagonal vs box distance) -- the classical, fast test.
// exact = true : separation tested on the actual point sets (hull diameter via pgl, exact minimum distance);
//                the box test is tried first, so exact accepts every pair that the box version accepts.
inline std::vector<Pair> callahanKosaraju(const Instance& I, bool exact) {
    std::vector<CKNode> T;
    std::function<int(std::vector<int>)> build = [&](std::vector<int> P) {
        CKNode nd; nd.x0 = nd.y0 = LLONG_MAX; nd.x1 = nd.y1 = LLONG_MIN;
        for (int i : P) { nd.x0 = std::min(nd.x0, I.x[i]); nd.x1 = std::max(nd.x1, I.x[i]); nd.y0 = std::min(nd.y0, I.y[i]); nd.y1 = std::max(nd.y1, I.y[i]); }
        nd.pts = P; int id = (int)T.size(); T.push_back(nd);
        if (P.size() > 1) {
            bool byx = (nd.x1 - nd.x0) >= (nd.y1 - nd.y0);
            ll lo = byx ? nd.x0 : nd.y0, hi = byx ? nd.x1 : nd.y1;
            // split at midpoint (fair split); points on the line go left
            std::vector<int> L, Rr;
            for (int i : P) { ll c = byx ? I.x[i] : I.y[i]; (2 * c <= lo + hi ? L : Rr).push_back(i); }
            if (Rr.empty()) { Rr.push_back(L.back()); L.pop_back(); }
            int l = build(L); int r = build(Rr); T[id].left = l; T[id].right = r;
        }
        return id;
    };
    std::vector<int> all(I.n()); std::iota(all.begin(), all.end(), 0);
    int root = build(all);
    auto diag2 = [&](const CKNode& a) { ll dx = a.x1 - a.x0, dy = a.y1 - a.y0; return dx * dx + dy * dy; };
    auto boxd2 = [&](const CKNode& a, const CKNode& b) {
        ll dx = std::max({(ll)0, a.x0 - b.x1, b.x0 - a.x1}), dy = std::max({(ll)0, a.y0 - b.y1, b.y0 - a.y1});
        return dx * dx + dy * dy;
    };
    // exact point-set diameters (squared), computed lazily from the convex hull
    std::vector<ll> pdiam2(T.size(), -1);
    auto pd2 = [&](int u) {
        if (pdiam2[u] >= 0) return pdiam2[u];
        const auto& P = T[u].pts; ll d = 0;
        if (P.size() == 2) d = I.d2(P[0], P[1]);
        else if (P.size() > 2) {
            std::vector<pgl::Point<ll>> pts; for (int i : P) pts.emplace_back(I.x[i], I.y[i]);
            pgl::Convex<pgl::Point<ll>> hull(pts);
            if (hull.size() >= 2) d = (ll)hull.diameter().squaredLength();
            else { for (size_t i = 0; i < P.size(); ++i) for (size_t j = i + 1; j < P.size(); ++j) d = std::max(d, I.d2(P[i], P[j])); }
        }
        return pdiam2[u] = d;
    };
    auto exactSeparated = [&](int u, int v) {
        ll r2 = std::max(pd2(u), pd2(v));
        // quick reject: d(Pu,Pv) <= distance of any two representatives
        if (!I.withinEps(r2, I.d2(T[u].pts[0], T[v].pts[0]))) return false;
        ll dm = LLONG_MAX;
        for (int x : T[u].pts) for (int y : T[v].pts) { dm = std::min(dm, I.d2(x, y)); if (!I.withinEps(r2, dm)) return false; }
        return I.withinEps(r2, dm);
    };
    std::vector<Pair> W;
    std::function<void(int,int)> findPairs = [&](int u, int v) {
        const CKNode &a = T[u], &b = T[v];
        ll r2 = std::max(diag2(a), diag2(b)), d2 = boxd2(a, b);
        if ((d2 > 0 && I.withinEps(r2, d2)) || (exact && exactSeparated(u, v))) {
            Pair pr{a.pts, b.pts}; std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
            W.push_back(pr); return;
        }
        if (diag2(a) >= diag2(b) && a.left >= 0) { findPairs(a.left, v); findPairs(a.right, v); }
        else if (b.left >= 0) { findPairs(u, b.left); findPairs(u, b.right); }
        else { findPairs(a.left, v); findPairs(a.right, v); }
    };
    std::function<void(int)> rec = [&](int u) {
        if (T[u].left < 0) return;
        findPairs(T[u].left, T[u].right); rec(T[u].left); rec(T[u].right);
    };
    rec(root);
    return W;
}

// greedy redundancy removal (keeps a valid cover): drop pairs whose elements are all covered by others.
// The counts saturate at 65535 (2 bytes per element). A stored count never exceeds the true count, and a pair is
// dropped only if all its stored counts are >= 2, so the cover stays valid; saturation can only keep a redundant pair.
inline std::vector<Pair> pruneRedundant(const Instance& I, std::vector<Pair> W) {
    int n = I.n(); ElemIndex EI(n);
    std::vector<uint16_t> cnt(EI.size(), 0);
    for (auto& p : W) for (int u : p.X) for (int v : p.Y) { auto& c = cnt[EI(u, v)]; if (c != UINT16_MAX) ++c; }
    std::vector<size_t> idx(W.size()); std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return W[a].X.size() * W[a].Y.size() < W[b].X.size() * W[b].Y.size(); });
    std::vector<char> keep(W.size(), 1);
    for (size_t k : idx) {
        bool red = true;
        for (int u : W[k].X) { for (int v : W[k].Y) if (cnt[EI(u, v)] < 2) { red = false; break; } if (!red) break; }
        if (red) { keep[k] = 0; for (int u : W[k].X) for (int v : W[k].Y) --cnt[EI(u, v)]; }
    }
    std::vector<Pair> out; for (size_t k = 0; k < W.size(); ++k) if (keep[k]) out.push_back(std::move(W[k]));
    return out;
}

} // namespace wspd

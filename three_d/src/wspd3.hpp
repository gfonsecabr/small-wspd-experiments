// Core routines for the experiments in three dimensions (separate from the planar code; nothing is shared).
// Points have integer coordinates, eps = a/b is rational, and all separation tests are exact:
//   max(diam X, diam Y) <= eps * d(X,Y)   <=>   b^2 * diam^2 <= a^2 * d^2   (128-bit integers).
// Same algorithms as the planar wspd.hpp, in R^3: the greedy sweep (nearest seed, fallback to the fixed collection of the grid
// pieces of the paper: 64 cubes in the ball, 32 in the half-ball), Callahan-Kosaraju on the fair split tree (box and
// exact test), exact OPT (set cover over all maximal dumbbells, CPLEX), redundancy removal and a verifier.
// Geometry: a static kd-tree written here (pgl is planar). ILPs: CPLEX (Concert).
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <ilcplex/ilocplex.h>

namespace wspd3 {

using ll = long long;
using i128 = __int128;

struct Eps { ll a, b; double val() const { return double(a) / double(b); } };

struct Pair { std::vector<int> X, Y; };      // sides sorted increasingly

struct Instance {
    std::vector<ll> x, y, z;
    Eps eps;
    int n() const { return (int)x.size(); }
    ll d2(int i, int j) const { ll dx = x[i] - x[j], dy = y[i] - y[j], dz = z[i] - z[j]; return dx * dx + dy * dy + dz * dz; }
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

// index of the element {i,j}, i != j, in a packed upper triangle (n(n-1)/2 entries)
struct ElemIndex {
    int n;
    explicit ElemIndex(int n_) : n(n_) {}
    size_t operator()(int i, int j) const { if (i > j) std::swap(i, j); return (size_t)i * (2 * (size_t)n - i - 1) / 2 + (j - i - 1); }
    size_t size() const { return (size_t)n * (n - 1) / 2; }
};

// ---------------------------------------------------------------- spatial index: static kd-tree
// Leaves of at most 8 points, children visited in a fixed order, so the points of a query are reported in an order
// that does not depend on the query (filterBall relies on this).
struct Spatial {
    const Instance& I;
    struct Node { ll lo[3], hi[3]; int left = -1, right = -1, begin = 0, end = 0; };
    std::vector<Node> T; std::vector<int> idx;
    ll coord(int i, int a) const { return a == 0 ? I.x[i] : a == 1 ? I.y[i] : I.z[i]; }
    int build(int b, int e) {
        Node nd; nd.begin = b; nd.end = e;
        for (int a = 0; a < 3; ++a) { nd.lo[a] = LLONG_MAX; nd.hi[a] = LLONG_MIN; }
        for (int k = b; k < e; ++k) for (int a = 0; a < 3; ++a) { ll c = coord(idx[k], a); nd.lo[a] = std::min(nd.lo[a], c); nd.hi[a] = std::max(nd.hi[a], c); }
        int id = (int)T.size(); T.push_back(nd);
        if (e - b > 8) {
            int ax = 0; for (int a = 1; a < 3; ++a) if (nd.hi[a] - nd.lo[a] > nd.hi[ax] - nd.lo[ax]) ax = a;
            int m = (b + e) / 2;
            std::nth_element(idx.begin() + b, idx.begin() + m, idx.begin() + e, [&](int u, int v) { ll cu = coord(u, ax), cv = coord(v, ax); return cu != cv ? cu < cv : u < v; });
            int l = build(b, m), r = build(m, e); T[id].left = l; T[id].right = r;
        }
        return id;
    }
    explicit Spatial(const Instance& I_) : I(I_) {
        idx.resize(I.n()); std::iota(idx.begin(), idx.end(), 0);
        if (I.n() > 0) build(0, I.n());
    }
    // indices of points z with |z c| <= eps * sqrt(D2) * factor   (factor = num/den, exact)
    std::vector<int> ball(int c, ll D2, ll num = 1, ll den = 1) const { std::vector<int> out; ball(c, D2, out, num, den); return out; }
    void ball(int c, ll D2, std::vector<int>& out, ll num = 1, ll den = 1) const {
        out.clear(); if (T.empty()) return;
        long double r = std::sqrt((long double)D2) * I.eps.val() * num / den;
        ll h = (ll)std::ceil(r) + 1;
        ll cc[3] = {I.x[c], I.y[c], I.z[c]};
        i128 lhsF = (i128)I.eps.b * I.eps.b * den * den, rhs = (i128)I.eps.a * I.eps.a * num * num * D2;
        std::function<void(int)> rec = [&](int u) {
            const Node& nd = T[u];
            for (int a = 0; a < 3; ++a) if (cc[a] + h < nd.lo[a] || cc[a] - h > nd.hi[a]) return;
            if (nd.left < 0) {
                for (int k = nd.begin; k < nd.end; ++k) { int zz = idx[k]; if (lhsF * I.d2(zz, c) <= rhs) out.push_back(zz); }
                return;
            }
            rec(nd.left); rec(nd.right);
        };
        rec(0);
    }
    // ball(c, D2) computed from a superset `sup` returned by ball(): same points in the same order
    void filterBall(const std::vector<int>& sup, int c, ll D2, std::vector<int>& out) const {
        out.clear();
        for (int zz : sup) if ((i128)I.eps.b * I.eps.b * I.d2(zz, c) <= (i128)I.eps.a * I.eps.a * D2) out.push_back(zz);
    }
};

// ---------------------------------------------------------------- maximal separated pairs with a given anchor
// The eps-separated dumbbells whose closest cross pair is (p,q) are the cliques containing p,q in the graph on
//   Xc = {x : |xp| <= eps|pq|, |xq| >= |pq|},  Yc = {y : |yq| <= eps|pq|, |yp| >= |pq|}
// with edges x-x' iff |xx'| <= eps|pq|, y-y' likewise, x-y iff |xy| >= |pq|. Maximal cliques by Bron-Kerbosch.
inline void anchoredMaximalPairs(const Instance& I, const Spatial& S, int p, int q, std::vector<Pair>& out,
                                 size_t cap = 200000, bool* capHit = nullptr) {
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
    std::vector<int> P0, X0;
    for (int i = 2; i < m; ++i) if (adj(i, 0) && adj(i, 1)) P0.push_back(i);
    std::vector<int> R = {0, 1};
    std::function<void(std::vector<int>&, std::vector<int>&, std::vector<int>&)> bk =
        [&](std::vector<int>& Rr, std::vector<int>& Pp, std::vector<int>& Xx) {
            if (out.size() >= cap) { if (capHit) *capHit = true; return; }
            if (Pp.empty() && Xx.empty()) {
                Pair pr;
                for (int i : Rr) (side[i] == 0 ? pr.X : pr.Y).push_back(V[i]);
                std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
                out.push_back(std::move(pr));
                return;
            }
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

// ---------------------------------------------------------------- set cover by CPLEX (MIP gap 0, one thread)
struct CoverResult { std::vector<int> chosen; double lb = 0; bool optimal = false; };

inline CoverResult setCoverCplex(int n, const std::vector<std::pair<int,int>>& elems, const std::vector<Pair>& cols, double timeLimit) {
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
    // exact brute force for tiny instances (as in the planar wspd.hpp): optimum 1 or 2 with <= 64 elements (LNS calls)
    if (elems.size() <= 64) {
        std::vector<uint64_t> mask(cols.size(), 0); uint64_t full = elems.size() == 64 ? ~0ull : ((1ull << elems.size()) - 1);
        for (size_t r = 0; r < rows.size(); ++r) for (int c : rows[r]) mask[c] |= (1ull << r);
        for (size_t c = 0; c < cols.size(); ++c) if (mask[c] == full) { res.chosen = {(int)c}; res.lb = 1; res.optimal = true; return res; }
        if (cols.size() <= 400)
            for (size_t c = 0; c < cols.size(); ++c) for (size_t c2 = c + 1; c2 < cols.size(); ++c2)
                if ((mask[c] | mask[c2]) == full) { res.chosen = {(int)c, (int)c2}; res.lb = 2; res.optimal = true; return res; }
    }
    IloEnv env;
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
        cp.setParam(IloCplex::Param::Threads, 1);
        cp.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
        cp.solve();
        for (size_t c = 0; c < cols.size(); ++c) if (cp.getValue(x[(IloInt)c]) > 0.5) res.chosen.push_back((int)c);
        res.lb = cp.getBestObjValue();
        res.optimal = (cp.getStatus() == IloAlgorithm::Optimal);
    } catch (IloException& ex) { std::cerr << "CPLEX: " << ex << std::endl; env.end(); throw; }
    env.end();
    return res;
}

// The largest set of pairwise incompatible pairs (alpha, the independence number of the compatibility graph), by CPLEX:
// one 0/1 variable per pair, at most one chosen pair in every maximal dumbbell (two pairs are compatible if and only if
// some maximal dumbbell covers both). As maxIncompatibleCplex of the planar wspd.hpp.
struct PackResult { std::vector<std::pair<int,int>> chosen; double ub = 0; bool optimal = false; };

inline PackResult maxIncompatibleCplex(int n, const std::vector<Pair>& cols, double timeLimit) {
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
        cp.setParam(IloCplex::Param::Threads, 1);
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

// all maximal separated pairs: union over all anchors, deduplicated. capHit is set if an anchor reached the cap (then
// the set cover over these columns is not the exact optimum)
inline std::vector<Pair> allMaximalPairs(const Instance& I, const Spatial& S, bool* capHit) {
    std::unordered_set<Pair, PairHash, PairEq> seen;
    std::vector<Pair> out, tmp;
    for (int p = 0; p < I.n(); ++p) for (int q = 0; q < I.n(); ++q) if (p != q) {
        tmp.clear(); anchoredMaximalPairs(I, S, p, q, tmp, 200000, capHit);
        for (auto& pr : tmp) { auto c = canonical(pr); if (seen.insert(c).second) out.push_back(c); }
    }
    return out;
}

// ---------------------------------------------------------------- verification
inline bool verify(const Instance& I, const std::vector<Pair>& W, std::string* why = nullptr) {
    int n = I.n(); ElemIndex EI(n);
    std::vector<bool> cov(EI.size(), false);
    for (auto& p : W) {
        if (!separated(I, p.X, p.Y)) { if (why) *why = "pair not separated"; return false; }
        for (int u : p.X) for (int v : p.Y) cov[EI(u, v)] = true;
    }
    for (size_t e = 0; e < cov.size(); ++e)
        if (!cov[e]) { if (why) *why = "element uncovered"; return false; }
    return true;
}

// ---------------------------------------------------------------- greedy growth of one dumbbell from an anchor
// As in the planar wspd.hpp: start from {a},{b} and repeatedly add the vertex of the anchor graph that is compatible with
// everything chosen and maximizes (#elements of R newly covered, #uncovered elements newly covered), smallest index on ties.
template <class IsR, class IsUncov>
inline Pair greedyAnchoredPair(const Instance& I, const Spatial& S, int a, int b, IsR inR, IsUncov uncov,
                               const std::vector<int>* Ba = nullptr, const std::vector<int>* Bb = nullptr) {
    ll D2 = I.d2(a, b);
    i128 lim128 = (i128)I.eps.a * I.eps.a * D2 / ((i128)I.eps.b * I.eps.b);
    ll lim = lim128 > (i128)LLONG_MAX ? LLONG_MAX : (ll)lim128;
    static thread_local struct { std::vector<int> V, C, cs[2], ballA, ballB; std::vector<char> side; std::vector<ll> vx, vy, vz; std::vector<int64_t> key; } sc;
    auto& V = sc.V; auto& side = sc.side; auto& ballA = sc.ballA; auto& ballB = sc.ballB;
    V.assign({a, b}); side.assign({0, 1});
    if (!Ba) { S.ball(a, D2, ballA); Ba = &ballA; }
    if (!Bb) { S.ball(b, D2, ballB); Bb = &ballB; }
    for (int x : *Ba) if (x != a && x != b && I.d2(x, b) >= D2) { V.push_back(x); side.push_back(0); }
    for (int y : *Bb) if (y != a && y != b && I.d2(y, a) >= D2) { V.push_back(y); side.push_back(1); }
    int m = (int)V.size();
    auto& vx = sc.vx; auto& vy = sc.vy; auto& vz = sc.vz; vx.resize(m); vy.resize(m); vz.resize(m);
    for (int i = 0; i < m; ++i) { vx[i] = I.x[V[i]]; vy[i] = I.y[V[i]]; vz[i] = I.z[V[i]]; }
    auto d2 = [&](int i, int j) { ll dx = vx[i] - vx[j], dy = vy[i] - vy[j], dz = vz[i] - vz[j]; return dx * dx + dy * dy + dz * dz; };
    auto adj = [&](int i, int j) { return side[i] == side[j] ? d2(i, j) <= lim : d2(i, j) >= D2; };
    auto& C = sc.C; auto& cs = sc.cs; C.assign({0, 1}); cs[0].clear(); cs[1].clear();
    for (int i = 2; i < m; ++i) if (adj(i, 0) && adj(i, 1)) cs[(int)side[i]].push_back(i);
    constexpr int SH = 21;
    if (m >= (1 << SH)) throw std::runtime_error("greedyAnchoredPair: too many candidates");
    auto& key = sc.key; key.resize(m);
    for (int v = 0; v < m; ++v) key[v] = ((int64_t)1 << SH) - 1 - v;
    auto addGain = [&](int v, int c) { key[v] += ((int64_t)inR(V[c], V[v]) << (2 * SH)) + ((int64_t)uncov(V[c], V[v]) << SH); };
    int u = -1; int64_t best = -1;
    for (int s = 0; s < 2; ++s) for (int v : cs[s]) { addGain(v, 1 - s); if (key[v] > best) { best = key[v]; u = v; } }
    while (u >= 0) {
        C.push_back(u);
        int su = side[u], nu = -1; best = -1;
        auto& L = cs[su]; size_t k = 0;
        for (int w : L) if (w != u && d2(u, w) <= lim) { L[k++] = w; if (key[w] > best) { best = key[w]; nu = w; } }
        L.resize(k);
        auto& M = cs[1 - su]; k = 0;
        for (int w : M) if (d2(u, w) >= D2) { addGain(w, u); M[k++] = w; if (key[w] > best) { best = key[w]; nu = w; } }
        M.resize(k);
        u = nu;
    }
    Pair pr; size_t nx = 0; for (int i : C) nx += side[i] == 0;
    pr.X.reserve(nx); pr.Y.reserve(C.size() - nx);
    for (int i : C) (side[i] == 0 ? pr.X : pr.Y).push_back(V[i]);
    std::sort(pr.X.begin(), pr.X.end()); std::sort(pr.Y.begin(), pr.Y.end());
    return pr;
}

// greedy set cover of the elements R (pairs (x,y), x on the X side) by columns; returns the chosen column indices
inline std::vector<int> greedySetCover(int n, const std::vector<std::pair<int,int>>& R, const std::vector<Pair>& cols) {
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
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        if (I.x[a] != I.x[b]) return I.x[a] < I.x[b];
        if (I.y[a] != I.y[b]) return I.y[a] < I.y[b];
        return I.z[a] < I.z[b]; });
    return ord;
}

// ---------------------------------------------------------------- grid pieces (theorem in R^3)
// The bounding cube of the ball B(c, rho), rho = eps sqrt(D2), is cut into 4^3 cubes of side rho/2 (diameter
// 0.866 rho <= (1-2eps) rho for eps <= 1/16); the half-ball (first coordinate >= that of c) meets the 2 x 4 x 4 = 32 cubes
// whose first index is 2 or 3. A coordinate offset of exactly rho goes to the last cube. Floating point, as in the plane.
inline int cubePiece(const Instance& I, int c, int zz, ll D2) {
    long double s = std::sqrt((long double)D2) * I.eps.val() / 2;   // cube side
    auto cell = [&](ll off) { int k = (int)std::floor((long double)off / s) + 2; return std::clamp(k, 0, 3); };
    return cell(I.x[zz] - I.x[c]) * 16 + cell(I.y[zz] - I.y[c]) * 4 + cell(I.z[zz] - I.z[c]);
}

struct SweepStats { long capHits = 0, steps = 0, fixedChosen = 0, adaptiveChosen = 0, candTotal = 0, rTotal = 0; double time = 0, tR = 0, tCand = 0, tCover = 0; };

// The greedy sweep of the paper (the planar wspd.hpp), in R^3. In the order <, for
// every point p and every q > p by increasing distance, a step with anchor {p,q} starts if {p,q} is uncovered. It covers
// the must-cover set R_pq (the uncovered pairs {x,y}, x, y >= p, with {{p,x},{q,y}} eps-separated) by dumbbells grown
// from (p,q) and then from the closest pair of R_pq still uncovered (if |R_pq| > 16; otherwise grown from every pair of
// R_pq and chosen by greedy set cover), and keeps them if they are not more numerous than the nonempty products of the
// fixed collection (eps <= 1/16), else the fixed collection. The anchors are pairwise incompatible, so the number of
// steps is a lower bound; they are appended to anchorsOut if given (check_lb3.cpp checks them).
inline std::vector<Pair> sweep(const Instance& I, SweepStats& st, size_t lazyThreshold = 16, std::vector<std::pair<int,int>>* anchorsOut = nullptr) {
    auto t0 = std::chrono::steady_clock::now();
    int n = I.n();
    Spatial S(I);
    auto ord = lexOrder(I);
    std::vector<int> rank(n); for (int k = 0; k < n; ++k) rank[ord[k]] = k;
    // one mark per element, stored symmetrically: COV = covered; INR = in R, not covered in this step; DONE = in R and
    // covered by a dumbbell grown in this step
    constexpr uint8_t INR = 1, DONE = 2, COV = 3;
    std::vector<uint8_t> mark((size_t)n * n, 0);
    auto MI = [n](int u, int v) { return (size_t)u * n + v; };
    auto setMark = [&](int u, int v, uint8_t m) { mark[MI(u, v)] = m; mark[MI(v, u)] = m; };
    auto isCov = [&](int u, int v) { return mark[MI(u, v)] == COV; };
    std::vector<Pair> W;
    auto coverPair = [&](const Pair& pr) { for (int u : pr.X) for (int v : pr.Y) setMark(u, v, COV); };
    bool fixedValid = (I.eps.a * 16 <= I.eps.b);
    std::vector<int> Bp, Bq, supP, supQ, ballA, ballB;
    std::vector<ll> distP(n);
    std::vector<std::pair<int,int>> R;
    std::array<std::vector<int>, 64> Apieces, Cpieces;   // half-ball pieces use the indices 32..63
    for (int k = 0; k < n; ++k) {
        int p = ord[k];
        std::vector<int> Q; for (int kk = k + 1; kk < n; ++kk) Q.push_back(ord[kk]);
        for (int zz : Q) distP[zz] = I.d2(p, zz);
        std::sort(Q.begin(), Q.end(), [&](int a, int b) { return distP[a] < distP[b]; });
        for (int q : Q) {
            if (isCov(p, q)) continue;
            ++st.steps;
            if (anchorsOut) anchorsOut->push_back({p, q});
            ll D2 = I.d2(p, q);
            auto uncoveredAbort = [&](int x, int y) {
                std::cerr << "sweep: the step of anchor {" << p << "," << q << "} leaves {" << x << "," << y << "} of its must-cover set uncovered" << std::endl;
                std::abort();
            };
            S.ball(p, D2, Bp); S.ball(q, D2, Bq);
            size_t fixedSize = 0;
            if (fixedValid) {
                for (auto& a : Apieces) a.clear();
                for (auto& c : Cpieces) c.clear();
                for (int zz : Bp) if (rank[zz] >= rank[p]) Apieces[cubePiece(I, p, zz, D2)].push_back(zz);
                for (int zz : Bq) Cpieces[cubePiece(I, q, zz, D2)].push_back(zz);
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
            auto tA = std::chrono::steady_clock::now();
            R.clear();
            for (int x : Bp) if (rank[x] >= rank[p]) for (int y : Bq) if (rank[y] >= rank[p] && x != y && !isCov(x, y)) {
                if (y == p || x == q) continue;
                ll diam2 = std::max(I.d2(p, x), I.d2(q, y));
                ll dmin2 = std::min(std::min(D2, I.d2(p, y)), std::min(I.d2(x, q), I.d2(x, y)));
                if (I.withinEps(diam2, dmin2)) R.push_back({x, y});
            }
            st.tR += seconds_since(tA); st.rTotal += R.size(); auto tB = std::chrono::steady_clock::now();
            for (auto& e : R) setMark(e.first, e.second, INR);
            auto inR = [&](int x, int y) { uint8_t m = mark[MI(x, y)]; return (long)(m == INR || m == DONE); };
            auto uncov = [&](int x, int y) { return (long)(mark[MI(x, y)] != COV); };
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
            std::vector<Pair> cols; std::vector<int> ch;
            bool useLazy = R.size() > lazyThreshold;
            if (useLazy) {
                size_t left = R.size();
                auto isLeft = [&](int x, int y) { return mark[MI(x, y)] == INR; };
                auto inLeft = [&](int x, int y) { return (long)isLeft(x, y); };
                std::pair<int,int> anc = {p, q}; size_t idx = 0;
                while (true) {
                    auto pr = growFrom(anc.first, anc.second, inLeft);
                    for (int u : pr.X) for (int v : pr.Y) if (isLeft(u, v)) { setMark(u, v, DONE); --left; }
                    ch.push_back((int)cols.size()); cols.push_back(std::move(pr));
                    if (left == 0) break;
                    // more grown dumbbells than the fixed collection: the fixed collection wins anyway (same output)
                    if (fixedSize > 0 && cols.size() >= fixedSize) {
                        ++st.capHits;
                        ch.clear(); cols = makeFixed(); for (size_t c = 0; c < cols.size(); ++c) ch.push_back((int)c);
                        break;
                    }
                    // next seed: the closest pair of R still uncovered (the first one in R among equals)
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
                std::unordered_set<Pair, PairHash, PairEq> seen;
                std::vector<std::pair<int,int>> anchors = R; anchors.push_back({p, q});
                for (auto [a, b] : anchors) { auto pr = growFrom(a, b, inR); if (seen.insert(pr).second) cols.push_back(pr); }
                st.tCand += seconds_since(tB); st.candTotal += cols.size();
                tC = std::chrono::steady_clock::now();
                ch = greedySetCover(n, R, cols);
            }
            st.tCover += seconds_since(tC);
            std::vector<Pair> chosen;
            if (fixedSize == 0 || ch.size() <= fixedSize) { for (int c : ch) chosen.push_back(std::move(cols[c])); ++st.adaptiveChosen; }
            else { ++st.fixedChosen; chosen = makeFixed(); }
            for (auto& pr : chosen) { coverPair(pr); W.push_back(std::move(pr)); }
            for (auto& e : R) if (!isCov(e.first, e.second)) uncoveredAbort(e.first, e.second);
            if (!isCov(p, q)) uncoveredAbort(p, q);
        }
    }
    st.time = seconds_since(t0);
    return W;
}

// ---------------------------------------------------------------- Callahan-Kosaraju WSPD (fair split tree)
// A node is split at the middle of the longest side of its bounding box (points on the plane go left). Two nodes are
// accepted if their bounding boxes are eps-separated (box test: larger box diagonal <= eps * box distance) or, with
// exact = true, also if their point sets are (exact diameters and closest distance). Otherwise the node whose box has
// the larger diagonal is split.
struct CKNode { ll lo[3], hi[3]; int left = -1, right = -1; std::vector<int> pts; };

inline std::vector<Pair> callahanKosaraju(const Instance& I, bool exact) {
    std::vector<CKNode> T;
    auto coord = [&](int i, int a) { return a == 0 ? I.x[i] : a == 1 ? I.y[i] : I.z[i]; };
    std::function<int(std::vector<int>)> build = [&](std::vector<int> P) {
        CKNode nd; for (int a = 0; a < 3; ++a) { nd.lo[a] = LLONG_MAX; nd.hi[a] = LLONG_MIN; }
        for (int i : P) for (int a = 0; a < 3; ++a) { nd.lo[a] = std::min(nd.lo[a], coord(i, a)); nd.hi[a] = std::max(nd.hi[a], coord(i, a)); }
        int ax = 0; for (int a = 1; a < 3; ++a) if (nd.hi[a] - nd.lo[a] > nd.hi[ax] - nd.lo[ax]) ax = a;
        ll lo = nd.lo[ax], hi = nd.hi[ax];
        nd.pts = std::move(P); int id = (int)T.size(); T.push_back(std::move(nd));
        if (T[id].pts.size() > 1) {
            std::vector<int> L, Rr;
            for (int i : T[id].pts) (2 * coord(i, ax) <= lo + hi ? L : Rr).push_back(i);
            if (Rr.empty()) { Rr.push_back(L.back()); L.pop_back(); }
            int l = build(std::move(L)); int r = build(std::move(Rr)); T[id].left = l; T[id].right = r;
        }
        return id;
    };
    std::vector<int> all(I.n()); std::iota(all.begin(), all.end(), 0);
    int root = build(all);
    auto diag2 = [&](const CKNode& a) { ll s = 0; for (int k = 0; k < 3; ++k) { ll d = a.hi[k] - a.lo[k]; s += d * d; } return s; };
    auto side2 = [&](const CKNode& a) { ll s = 0; for (int k = 0; k < 3; ++k) { ll d = a.hi[k] - a.lo[k]; s = std::max(s, d * d); } return s; };
    auto boxd2 = [&](const CKNode& a, const CKNode& b) {
        ll s = 0; for (int k = 0; k < 3; ++k) { ll d = std::max({(ll)0, a.lo[k] - b.hi[k], b.lo[k] - a.hi[k]}); s += d * d; } return s; };
    // exact diameters (squared), computed lazily by brute force
    std::vector<ll> pdiam2(T.size(), -1);
    auto pd2 = [&](int u) {
        if (pdiam2[u] >= 0) return pdiam2[u];
        const auto& P = T[u].pts; ll d = 0;
        for (size_t i = 0; i < P.size(); ++i) for (size_t j = i + 1; j < P.size(); ++j) d = std::max(d, I.d2(P[i], P[j]));
        return pdiam2[u] = d;
    };
    auto exactSeparated = [&](int u, int v) {
        ll rep = I.d2(T[u].pts[0], T[v].pts[0]);
        // quick reject: the longest box side is at most the diameter, and the closest distance at most rep
        if (!I.withinEps(std::max(side2(T[u]), side2(T[v])), rep)) return false;
        ll r2 = std::max(pd2(u), pd2(v));
        if (!I.withinEps(r2, rep)) return false;
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
            W.push_back(std::move(pr)); return;
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

// greedy redundancy removal (keeps a valid cover), smallest dumbbells first; counts saturate at 65535 (safe)
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

} // namespace wspd3

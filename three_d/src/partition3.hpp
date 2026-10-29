// Partition WSPDs in R^3 (the planar partition.hpp; every pair covered exactly once). toPartition, the merge pass, and
// verifyPartition are those of the planar partition.hpp (they do not depend on the dimension); quadtreePartition is the octree
// version (cubes of diameter sqrt3 times the side, diameters by brute force).
//   toPartition       converts any cover WSPD into a partition WSPD: the dumbbells are processed in some order,
//                     and the pairs of a dumbbell that are not yet covered are split into dumbbells whose sides are
//                     subsets of its sides (separated by heredity). Optional merge pass afterwards.
//   verifyPartition   exact separation of every dumbbell, every pair covered exactly once.
//   quadtreePartition the maximal same-level octree cell dumbbells (the construction of thm:partition in R^3), with
//                     the separation tested on the points (exact) or on the cells (box), deterministic or randomly
//                     shifted and scaled.
#pragma once
#include <cassert>
#include <climits>
#include <map>
#include <queue>
#include <tuple>
#include <unordered_map>
#include "wspd3.hpp"

namespace wspd3 {

// ---------------------------------------------------------------- bitset over the elements (packed triangle)
struct ElemBits {
    std::vector<uint64_t> w;
    explicit ElemBits(size_t m) : w((m + 63) / 64, 0) {}
    bool get(size_t e) const { return (w[e >> 6] >> (e & 63)) & 1; }
    void set(size_t e) { w[e >> 6] |= uint64_t(1) << (e & 63); }
};

struct PartConfig {
    std::string order = "big";     // big (largest |X||Y| first) | small | input | random | fewest (fewest pieces first)
    std::string split = "min";     // rows | cols | min (fewer of rows and cols) | greedy (greedy bicliques, min of both)
    bool merge = true;             // merge dumbbells that share a side when the union stays separated
    unsigned seed = 1;
    std::string name() const { return order + "/" + split + "/" + (merge ? "merge" : "nomerge"); }
    static PartConfig parse(const std::string& s) {   // "order/split/merge|nomerge"; missing fields keep the default
        PartConfig c; std::vector<std::string> f; std::string cur;
        for (char ch : s) { if (ch == '/') { f.push_back(cur); cur.clear(); } else cur += ch; }
        f.push_back(cur);
        if (f.size() > 0 && !f[0].empty()) c.order = f[0];
        if (f.size() > 1 && !f[1].empty()) c.split = f[1];
        if (f.size() > 2 && !f[2].empty()) c.merge = f[2] != "nomerge";
        return c;
    }
};

struct PartStats {
    size_t coverSize = 0, dropped = 0, split = 0, maxPieces = 0, beforeMerge = 0, merges = 0;
    double overlap = 0, time = 0, tMerge = 0;   // overlap = sum |X||Y| / C(n,2)
};

namespace detail {

// Uncovered pairs of {X,Y} as rows: for row r (a point of A), the positions j in B with {A[r],B[j]} uncovered.
inline std::vector<std::vector<int>> uncoveredRows(const std::vector<int>& A, const std::vector<int>& B, const ElemBits& cov, const ElemIndex& EI) {
    std::vector<std::vector<int>> R(A.size());
    for (size_t r = 0; r < A.size(); ++r) for (size_t j = 0; j < B.size(); ++j) if (!cov.get(EI(A[r], B[j]))) R[r].push_back((int)j);
    return R;
}

struct VecHash { size_t operator()(const std::vector<int>& v) const { size_t h = v.size(); for (int x : v) h = h * 1000003u ^ (size_t)x; return h; } };

// Rows with identical neighbourhoods form one biclique.
inline std::vector<Pair> rowClasses(const std::vector<int>& A, const std::vector<int>& B, const std::vector<std::vector<int>>& R) {
    std::unordered_map<std::vector<int>, int, VecHash> id; std::vector<Pair> out;
    for (size_t r = 0; r < A.size(); ++r) {
        if (R[r].empty()) continue;
        auto [it, fresh] = id.try_emplace(R[r], (int)out.size());
        if (fresh) { Pair p; for (int j : R[r]) p.Y.push_back(B[j]); out.push_back(std::move(p)); }
        out[it->second].X.push_back(A[r]);
    }
    return out;
}

// Greedy bicliques: among the candidates C = N(r), take the one maximizing |{r' : N(r') ⊇ C}| * |C|, remove it, repeat.
// O(rows^3 |B| / 64); falls back to row classes for more than `cap` nonempty rows.
inline std::vector<Pair> greedyBicliques(const std::vector<int>& A, const std::vector<int>& B, const std::vector<std::vector<int>>& R, size_t cap = 400) {
    std::vector<int> rows; for (size_t r = 0; r < A.size(); ++r) if (!R[r].empty()) rows.push_back((int)r);
    if (rows.size() > cap) return rowClasses(A, B, R);
    size_t W = (B.size() + 63) / 64;
    std::vector<std::vector<uint64_t>> N(rows.size(), std::vector<uint64_t>(W, 0));
    for (size_t k = 0; k < rows.size(); ++k) for (int j : R[rows[k]]) N[k][j >> 6] |= uint64_t(1) << (j & 63);
    auto empty = [&](size_t k) { for (auto w : N[k]) if (w) return false; return true; };
    auto contains = [&](size_t k, size_t c) { for (size_t t = 0; t < W; ++t) if ((N[c][t] & ~N[k][t]) != 0) return false; return true; };   // N(k) ⊇ N(c)
    auto pop = [&](size_t k) { size_t s = 0; for (auto w : N[k]) s += std::popcount(w); return s; };
    std::vector<Pair> out;
    while (true) {
        size_t best = SIZE_MAX, bestVal = 0;
        for (size_t c = 0; c < rows.size(); ++c) {
            if (empty(c)) continue;
            size_t cnt = 0; for (size_t k = 0; k < rows.size(); ++k) if (!empty(k) && contains(k, c)) ++cnt;
            size_t v = cnt * pop(c); if (v > bestVal) { bestVal = v; best = c; }
        }
        if (best == SIZE_MAX) break;
        std::vector<uint64_t> C = N[best]; Pair p;
        for (size_t j = 0; j < B.size(); ++j) if ((C[j >> 6] >> (j & 63)) & 1) p.Y.push_back(B[j]);
        for (size_t k = 0; k < rows.size(); ++k) {
            bool sup = !empty(k); for (size_t t = 0; t < W && sup; ++t) if ((C[t] & ~N[k][t]) != 0) sup = false;
            if (!sup) continue;
            p.X.push_back(A[rows[k]]); for (size_t t = 0; t < W; ++t) N[k][t] &= ~C[t];
        }
        out.push_back(std::move(p));
    }
    return out;
}

inline std::vector<int> unionSorted(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> u; u.reserve(a.size() + b.size()); std::merge(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(u)); return u;
}

// Merge pass: dumbbells {S,O1}, {S,O2} with the same side S cover disjoint pairs, so O1 and O2 are disjoint and
// {S, O1 ∪ O2} covers exactly their pairs; replace them by it if it is separated. Repeated until no merge (at most
// 20 passes). Grids of dumbbells {Xi,Yj} merge in several steps, since every intermediate union is a subdumbbell of
// the final one. Sides are grouped by a 64-bit hash (sorted, so the order is deterministic; equality is checked).
// Per dumbbell, the squared diameters of both sides and the squared distance are cached, so a merge test scans only
// the cross distances of O1 and O2. A group tries at most groupCap dumbbells per pass.
inline size_t mergePass(const Instance& I, std::vector<Pair>& W, size_t groupCap = 256) {
    auto diam2 = [&](const std::vector<int>& A) { ll d = 0; for (size_t i = 0; i < A.size(); ++i) for (size_t j = i + 1; j < A.size(); ++j) d = std::max(d, I.d2(A[i], A[j])); return d; };
    struct Info { ll dX, dY, dm; };
    std::vector<Info> info(W.size());
    for (size_t i = 0; i < W.size(); ++i) {
        ll dm = LLONG_MAX; for (int u : W[i].X) for (int v : W[i].Y) dm = std::min(dm, I.d2(u, v));
        info[i] = {diam2(W[i].X), diam2(W[i].Y), dm};
    }
    auto hashOf = [](const std::vector<int>& v) { uint64_t h = 0x9e3779b97f4a7c15ull ^ v.size(); for (int x : v) { h ^= (uint64_t)(uint32_t)x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2); h *= 0xff51afd7ed558ccdull; } return h; };
    struct Entry { uint64_t h; uint32_t i; uint8_t which; };   // which = 0: S = X, 1: S = Y
    size_t merges = 0;
    std::vector<char> alive(W.size(), 1);
    for (int pass = 0; pass < 20; ++pass) {
        std::vector<Entry> E;
        for (size_t i = 0; i < W.size(); ++i) if (alive[i]) { E.push_back({hashOf(W[i].X), (uint32_t)i, 0}); E.push_back({hashOf(W[i].Y), (uint32_t)i, 1}); }
        std::sort(E.begin(), E.end(), [](const Entry& a, const Entry& b) { return a.h != b.h ? a.h < b.h : a.i != b.i ? a.i < b.i : a.which < b.which; });
        std::vector<char> touched(W.size(), 0); size_t before = merges;
        auto side = [&](const Entry& e) -> std::vector<int>& { return e.which ? W[e.i].Y : W[e.i].X; };
        auto other = [&](const Entry& e) -> std::vector<int>& { return e.which ? W[e.i].X : W[e.i].Y; };
        auto otherDiam = [&](const Entry& e) -> ll& { return e.which ? info[e.i].dX : info[e.i].dY; };
        for (size_t g0 = 0, g1; g0 < E.size(); g0 = g1) {
            for (g1 = g0 + 1; g1 < E.size() && E[g1].h == E[g0].h; ++g1) {}
            if (g1 - g0 < 2) continue;
            std::vector<char> used(g1 - g0, 0);
            for (size_t r = g0; r < g1; ++r) {   // one class of equal sides per representative (hash collisions)
                if (used[r - g0]) continue;
                std::vector<Entry> cand;
                for (size_t k = r; k < g1; ++k) if (!used[k - g0] && side(E[k]) == side(E[r])) {
                    used[k - g0] = 1;
                    if (!touched[E[k].i] && alive[E[k].i] && cand.size() < groupCap) cand.push_back(E[k]);
                }
                if (cand.size() < 2) continue;
                std::stable_sort(cand.begin(), cand.end(), [&](const Entry& a, const Entry& b) { return other(a).size() > other(b).size(); });
                ll dS = cand[0].which ? info[cand[0].i].dY : info[cand[0].i].dX;
                std::vector<Entry> heads;
                for (auto& e : cand) {
                    bool done = false;
                    for (auto& h : heads) {
                        ll dm = std::min(info[h.i].dm, info[e.i].dm);
                        ll dU = std::max(otherDiam(h), otherDiam(e));
                        if (!I.withinEps(std::max(dS, dU), dm)) continue;
                        bool ok = true;
                        for (int u : other(h)) { for (int v : other(e)) { dU = std::max(dU, I.d2(u, v)); if (!I.withinEps(dU, dm)) { ok = false; break; } } if (!ok) break; }
                        if (!ok) continue;
                        other(h) = unionSorted(other(h), other(e)); otherDiam(h) = dU; info[h.i].dm = dm;
                        alive[e.i] = 0; touched[h.i] = touched[e.i] = 1; ++merges; done = true; break;
                    }
                    if (!done) heads.push_back(e);
                }
            }
        }
        if (merges == before) break;
    }
    std::vector<Pair> out; for (size_t i = 0; i < W.size(); ++i) if (alive[i]) out.push_back(std::move(W[i]));
    W = std::move(out);
    return merges;
}

} // namespace detail

// W is taken by value and consumed dumbbell by dumbbell (at most one copy of the pairs in memory).
inline std::vector<Pair> toPartition(const Instance& I, std::vector<Pair> W, const PartConfig& cf, PartStats& st) {
    auto t0 = std::chrono::steady_clock::now();
    int n = I.n(); ElemIndex EI(n); ElemBits cov(EI.size());
    st.coverSize = W.size();
    double tot = 0; for (auto& p : W) tot += double(p.X.size()) * p.Y.size();
    st.overlap = tot / (double(n) * (n - 1) / 2);
    auto sz = [&](size_t i) { return W[i].X.size() * W[i].Y.size(); };
    // the uncovered pairs of p, split into dumbbells (empty if p covers no new pair)
    auto splitOf = [&](const Pair& p) {
        std::vector<Pair> pieces;
        auto RX = detail::uncoveredRows(p.X, p.Y, cov, EI);
        size_t unc = 0; for (auto& r : RX) unc += r.size();
        if (unc == 0) return pieces;
        if (unc == p.X.size() * p.Y.size()) { pieces.push_back(p); return pieces; }
        auto RY = detail::uncoveredRows(p.Y, p.X, cov, EI);
        auto swapSides = [](std::vector<Pair> v) { for (auto& q : v) std::swap(q.X, q.Y); return v; };
        if (cf.split == "rows") pieces = detail::rowClasses(p.X, p.Y, RX);
        else if (cf.split == "cols") pieces = swapSides(detail::rowClasses(p.Y, p.X, RY));
        else if (cf.split == "min" || cf.split == "greedy") {
            bool g = cf.split == "greedy";
            auto a = g ? detail::greedyBicliques(p.X, p.Y, RX) : detail::rowClasses(p.X, p.Y, RX);
            auto b = swapSides(g ? detail::greedyBicliques(p.Y, p.X, RY) : detail::rowClasses(p.Y, p.X, RY));
            pieces = b.size() < a.size() ? std::move(b) : std::move(a);
        } else { std::cerr << "unknown split " << cf.split << std::endl; exit(1); }
        return pieces;
    };
    std::vector<Pair> out;
    auto commit = [&](size_t i, std::vector<Pair> pieces) {
        if (pieces.empty()) { ++st.dropped; }
        else {
            if (pieces.size() > 1 || pieces[0].X.size() * pieces[0].Y.size() < sz(i)) ++st.split;
            st.maxPieces = std::max(st.maxPieces, pieces.size());
            for (auto& q : pieces) for (int u : q.X) for (int v : q.Y) cov.set(EI(u, v));
            for (auto& q : pieces) out.push_back(std::move(q));
        }
        W[i] = Pair{};
    };
    if (cf.order == "fewest") {
        // lazy priority queue on (number of pieces, -|X||Y|): a popped dumbbell whose number of pieces grew since it
        // was pushed goes back with the new number. Initially every dumbbell has one piece, so ties are broken as in "big".
        using Key = std::tuple<size_t, size_t, size_t>;   // pieces, -size (as SIZE_MAX - size), index
        std::priority_queue<Key, std::vector<Key>, std::greater<Key>> pq;
        for (size_t i = 0; i < W.size(); ++i) pq.push({1, SIZE_MAX - sz(i), i});
        while (!pq.empty()) {
            auto [c, ns, i] = pq.top(); pq.pop();
            auto pieces = splitOf(W[i]);
            if (pieces.size() > c) { pq.push({pieces.size(), ns, i}); continue; }
            commit(i, std::move(pieces));
        }
    } else {
        std::vector<size_t> ord(W.size()); std::iota(ord.begin(), ord.end(), 0);
        if (cf.order == "big") std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return sz(a) > sz(b); });
        else if (cf.order == "small") std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return sz(a) < sz(b); });
        else if (cf.order == "random") { std::mt19937_64 rng(cf.seed); std::shuffle(ord.begin(), ord.end(), rng); }
        else if (cf.order != "input") { std::cerr << "unknown order " << cf.order << std::endl; exit(1); }
        for (size_t i : ord) commit(i, splitOf(W[i]));
    }
    st.beforeMerge = out.size();
    if (cf.merge) { auto tm = std::chrono::steady_clock::now(); st.merges = detail::mergePass(I, out); st.tMerge = seconds_since(tm); }
    st.time = seconds_since(t0);
    return out;
}

// Exact separation of every dumbbell and every pair covered exactly once.
inline bool verifyPartition(const Instance& I, const std::vector<Pair>& W, std::string* why = nullptr) {
    int n = I.n(); ElemIndex EI(n); ElemBits cov(EI.size());
    for (auto& p : W) {
        if (!separated(I, p.X, p.Y)) { if (why) *why = "pair not separated"; return false; }
        for (int u : p.X) for (int v : p.Y) { size_t e = EI(u, v); if (cov.get(e)) { if (why) *why = "element covered twice"; return false; } cov.set(e); }
    }
    for (size_t e = 0; e < EI.size(); ++e) if (!cov.get(e)) { if (why) *why = "element uncovered"; return false; }
    return true;
}

// ---------------------------------------------------------------- octree partition WSPD (thm:partition in R^3)
// As quadtreePartition of the planar partition.hpp, with cubes: root cube o + [0,S)^3, P in m + [0,L)^3 (L = largest side of the
// bounding box + 1). A point gets the integer coordinates A = floor((x - o)/S * 2^K), K = 62, so the cells of all levels
// nest exactly. Visit(c,c') outputs {P_c,P_c'} if separated, else recurses into the cells of the first level at which
// P_c or P_c' meets two cells.
//   exact = true : separation tested on the points;  exact = false: on the cubes (diameter sqrt3 s vs. cube distance).
struct QTStats { long calls = 0, levels = 0; double S = 0, L = 0; };

inline std::vector<Pair> quadtreePartition(const Instance& I, bool exact, bool randomized, unsigned seed, QTStats& st) {
    const int K = 62; int n = I.n();
    const std::vector<ll>* C[3] = {&I.x, &I.y, &I.z};
    ll m[3], M[3];
    for (int a = 0; a < 3; ++a) { m[a] = *std::min_element(C[a]->begin(), C[a]->end()); M[a] = *std::max_element(C[a]->begin(), C[a]->end()); }
    ll Li = std::max({M[0] - m[0], M[1] - m[1], M[2] - m[2]}) + 1;
    long double L = (long double)Li, S = L, o[3] = {(long double)m[0], (long double)m[1], (long double)m[2]};
    if (randomized) {
        std::mt19937_64 rng(seed); std::uniform_real_distribution<long double> U(0, 1);
        S = std::pow((long double)2, 1 + U(rng)) * L;
        for (int a = 0; a < 3; ++a) o[a] = m[a] - U(rng) * S / 2;
    }
    st.S = (double)S; st.L = (double)L;
    std::vector<std::array<uint64_t, 3>> A(n);
    long double scale = std::ldexp((long double)1, K) / S;
    for (int i = 0; i < n; ++i) for (int a = 0; a < 3; ++a) {
        ll c = (*C[a])[i];
        A[i][a] = randomized ? (uint64_t)std::floor((c - o[a]) * scale) : (uint64_t)(((i128)(c - m[a]) << K) / Li);
        assert(A[i][a] < (uint64_t(1) << K));
    }
    auto splitLevel = [&](const std::vector<int>& P) {
        uint64_t d = 0; for (int i : P) for (int a = 0; a < 3; ++a) d |= A[i][a] ^ A[P[0]][a];
        if (!d) return K + 1;
        return K - (63 - std::countl_zero(d));
    };
    using Cell = std::array<uint64_t, 3>;
    auto cellOf = [&](int i, int j) { return Cell{A[i][0] >> (K - j), A[i][1] >> (K - j), A[i][2] >> (K - j)}; };
    auto cellsAt = [&](const std::vector<int>& P, int j) {
        std::map<Cell, std::vector<int>> mp;
        for (int i : P) mp[cellOf(i, j)].push_back(i);
        std::vector<std::vector<int>> out; for (auto& kv : mp) out.push_back(std::move(kv.second));
        return out;
    };
    auto diam2 = [&](const std::vector<int>& P) -> ll {
        ll d = 0; for (size_t a = 0; a < P.size(); ++a) for (size_t b = a + 1; b < P.size(); ++b) d = std::max(d, I.d2(P[a], P[b]));
        return d;
    };
    auto exactSep = [&](const std::vector<int>& X, const std::vector<int>& Y) {
        ll r2 = std::max(diam2(X), diam2(Y));
        if (!I.withinEps(r2, I.d2(X[0], Y[0]))) return false;
        ll dm = LLONG_MAX;
        for (int x : X) for (int y : Y) { dm = std::min(dm, I.d2(x, y)); if (!I.withinEps(r2, dm)) return false; }
        return true;
    };
    // box test on the cubes of a level, exact: in units of the cell side, the cubes have diameter sqrt3 and the axis gaps
    // are integers, so the test is 3 b^2 <= a^2 (gx^2 + gy^2 + gz^2). Gaps >= 2^40 always pass.
    auto boxSep = [&](const Cell& c, const Cell& d) {
        auto gap = [](uint64_t a, uint64_t b) -> uint64_t { return a > b ? a - b - 1 : b > a ? b - a - 1 : 0; };
        uint64_t g[3]; for (int a = 0; a < 3; ++a) g[a] = gap(c[a], d[a]);
        if (std::max({g[0], g[1], g[2]}) >= (uint64_t(1) << 40)) return true;
        return (i128)3 * I.eps.b * I.eps.b <= (i128)I.eps.a * I.eps.a * ((i128)g[0] * g[0] + (i128)g[1] * g[1] + (i128)g[2] * g[2]);
    };
    std::vector<Pair> W;
    std::function<void(std::vector<int>, std::vector<int>, int)> visit = [&](std::vector<int> X, std::vector<int> Y, int j) {
        ++st.calls; st.levels = std::max<long>(st.levels, j);
        int k = std::min(splitLevel(X), splitLevel(Y));
        bool sep = false;
        if (exact) sep = exactSep(X, Y);
        else for (int l = j; l < std::min(k, K + 1) && !sep; ++l) sep = boxSep(cellOf(X[0], l), cellOf(Y[0], l));
        if (sep) { std::sort(X.begin(), X.end()); std::sort(Y.begin(), Y.end()); W.push_back(Pair{std::move(X), std::move(Y)}); return; }
        assert(k <= K);
        auto CX = cellsAt(X, k), CY = cellsAt(Y, k);
        for (auto& a : CX) for (auto& b : CY) visit(a, b, k);
    };
    std::function<void(std::vector<int>)> node = [&](std::vector<int> P) {
        int k = splitLevel(P); if (k > K) return;
        auto Cs = cellsAt(P, k);
        for (size_t a = 0; a < Cs.size(); ++a) for (size_t b = a + 1; b < Cs.size(); ++b) visit(Cs[a], Cs[b], k);
        for (auto& c : Cs) node(std::move(c));
    };
    std::vector<int> all(n); std::iota(all.begin(), all.end(), 0);
    node(all);
    return W;
}

} // namespace wspd3

// Partition WSPDs (every pair covered exactly once).
//   toPartition       converts any cover WSPD into a partition WSPD: the dumbbells are processed in the order of fewest
//                     pieces, and the pairs of a dumbbell that are not yet covered are split into dumbbells whose sides
//                     are subsets of its sides (separated by heredity); then the merge pass.
//   verifyPartition   exact separation of every dumbbell, every pair covered exactly once.
//   quadtreePartition the maximal same-level quadtree cell dumbbells of Theorem partition, with the separation tested
//                     on the points, deterministic or randomly shifted and scaled.
#pragma once
#include <cassert>
#include <climits>
#include <map>
#include <queue>
#include <tuple>
#include <unordered_map>
#include "wspd.hpp"

namespace wspd {

// ---------------------------------------------------------------- bitset over the elements (packed triangle)
struct ElemBits {
    std::vector<uint64_t> w;
    explicit ElemBits(size_t m) : w((m + 63) / 64, 0) {}
    bool get(size_t e) const { return (w[e >> 6] >> (e & 63)) & 1; }
    void set(size_t e) { w[e >> 6] |= uint64_t(1) << (e & 63); }
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
// Order "fewest pieces": a lazy priority queue on (number of pieces of the uncovered pairs, -|X||Y|); split: rows (or
// columns, whichever gives fewer pieces) with equal uncovered neighbourhoods form one piece.
inline std::vector<Pair> toPartition(const Instance& I, std::vector<Pair> W, PartStats& st) {
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
        auto a = detail::rowClasses(p.X, p.Y, RX);
        auto b = swapSides(detail::rowClasses(p.Y, p.X, RY));
        pieces = b.size() < a.size() ? std::move(b) : std::move(a);
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
    {
        // a popped dumbbell whose number of pieces grew since it was pushed goes back with the new number. Initially
        // every dumbbell has one piece, so ties are broken by size (largest first), then by index.
        using Key = std::tuple<size_t, size_t, size_t>;   // pieces, -size (as SIZE_MAX - size), index
        std::priority_queue<Key, std::vector<Key>, std::greater<Key>> pq;
        for (size_t i = 0; i < W.size(); ++i) pq.push({1, SIZE_MAX - sz(i), i});
        while (!pq.empty()) {
            auto [c, ns, i] = pq.top(); pq.pop();
            auto pieces = splitOf(W[i]);
            if (pieces.size() > c) { pq.push({pieces.size(), ns, i}); continue; }
            commit(i, std::move(pieces));
        }
    }
    st.beforeMerge = out.size();
    auto tm = std::chrono::steady_clock::now(); st.merges = detail::mergePass(I, out); st.tMerge = seconds_since(tm);
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

// ---------------------------------------------------------------- quadtree partition WSPD (thm:partition)
// Root square R0 = o + [0,S)^2 with P ⊂ m + [0,L)^2 (m = lower-left corner of the bounding box, L = its larger side
// + 1, so P ⊂ [0,L)^2 after translation). Deterministic: o = m, S = L. Randomized: U uniform in [0,1), S = 2^{1+U} L,
// o = m + uniform in (-S/2,0]^2. A point gets the integer coordinates A = floor((x - o)/S * 2^K), K = 62, computed once;
// its cell of level i has index A >> (K - i), so the cells of all levels are nested exactly.
// A cell dumbbell {c,c'} (same level, nonempty) is maximal if no pair of ancestors is a cell dumbbell. The recursion
// follows the proof: Visit(c,c') outputs {P_c,P_c'} if separated, else recurses into the cells of the first level at
// which P_c or P_c' meets two cells. Separation is tested on the points.
struct QTStats { long calls = 0, levels = 0; double S = 0, L = 0; };

inline std::vector<Pair> quadtreePartition(const Instance& I, bool randomized, unsigned seed, QTStats& st) {
    const int K = 62; int n = I.n();
    ll mx = *std::min_element(I.x.begin(), I.x.end()), my = *std::min_element(I.y.begin(), I.y.end());
    ll Mx = *std::max_element(I.x.begin(), I.x.end()), My = *std::max_element(I.y.begin(), I.y.end());
    long double L = (long double)std::max(Mx - mx, My - my) + 1, S = L, ox = mx, oy = my;
    if (randomized) {
        std::mt19937_64 rng(seed); std::uniform_real_distribution<long double> U(0, 1);
        S = std::pow((long double)2, 1 + U(rng)) * L;
        ox = mx - U(rng) * S / 2; oy = my - U(rng) * S / 2;   // (-S/2,0]: the endpoint has measure zero
    }
    st.S = (double)S; st.L = (double)L;
    std::vector<uint64_t> AX(n), AY(n);
    long double scale = std::ldexp((long double)1, K) / S;
    ll Li = std::max(Mx - mx, My - my) + 1;
    for (int i = 0; i < n; ++i) {
        if (!randomized) {   // exact: floor((x - m) 2^K / L) in integers, so points on cell boundaries go to the right cell
            AX[i] = (uint64_t)(((i128)(I.x[i] - mx) << K) / Li); AY[i] = (uint64_t)(((i128)(I.y[i] - my) << K) / Li);
        } else {             // the boundaries are at random real positions (measure zero)
            long double tx = (I.x[i] - ox) * scale, ty = (I.y[i] - oy) * scale;
            AX[i] = (uint64_t)std::floor(tx); AY[i] = (uint64_t)std::floor(ty);
        }
        assert(AX[i] < (uint64_t(1) << K) && AY[i] < (uint64_t(1) << K));
    }
    // first level at which the points of P meet two cells (K + 1 if |P| = 1)
    auto splitLevel = [&](const std::vector<int>& P) {
        uint64_t dx = 0, dy = 0; for (int i : P) { dx |= AX[i] ^ AX[P[0]]; dy |= AY[i] ^ AY[P[0]]; }
        uint64_t d = dx | dy; if (!d) return K + 1;
        return K - (63 - std::countl_zero(d));
    };
    // nonempty cells of level j inside P (P lies in one cell of level j-1 or above)
    auto cellsAt = [&](const std::vector<int>& P, int j) {
        std::map<std::pair<uint64_t,uint64_t>, std::vector<int>> m;
        for (int i : P) m[{AX[i] >> (K - j), AY[i] >> (K - j)}].push_back(i);
        std::vector<std::vector<int>> out; for (auto& kv : m) out.push_back(std::move(kv.second));
        return out;
    };
    // squared diameter of a point set via the convex hull
    auto diam2 = [&](const std::vector<int>& P) -> ll {
        if (P.size() < 2) return 0;
        if (P.size() == 2) return I.d2(P[0], P[1]);
        std::vector<pgl::Point<ll>> pts; for (int i : P) pts.emplace_back(I.x[i], I.y[i]);
        pgl::Convex<pgl::Point<ll>> hull(pts);
        if (hull.size() >= 2) return (ll)hull.diameter().squaredLength();
        ll d = 0; for (size_t a = 0; a < P.size(); ++a) for (size_t b = a + 1; b < P.size(); ++b) d = std::max(d, I.d2(P[a], P[b]));
        return d;
    };
    auto exactSep = [&](const std::vector<int>& A, const std::vector<int>& B) {
        ll r2 = std::max(diam2(A), diam2(B));
        if (!I.withinEps(r2, I.d2(A[0], B[0]))) return false;
        ll dm = LLONG_MAX;
        for (int x : A) for (int y : B) { dm = std::min(dm, I.d2(x, y)); if (!I.withinEps(r2, dm)) return false; }
        return true;
    };
    std::vector<Pair> W;
    std::function<void(std::vector<int>, std::vector<int>, int)> visit = [&](std::vector<int> A, std::vector<int> B, int j) {
        ++st.calls; st.levels = std::max<long>(st.levels, j);
        int k = std::min(splitLevel(A), splitLevel(B));
        // up to level k-1, A and B stay in one cell each, so the test sees the same point sets at every such level
        if (exactSep(A, B)) { std::sort(A.begin(), A.end()); std::sort(B.begin(), B.end()); W.push_back(Pair{std::move(A), std::move(B)}); return; }
        assert(k <= K);
        auto CA = cellsAt(A, k), CB = cellsAt(B, k);
        for (auto& a : CA) for (auto& b : CB) visit(a, b, k);
    };
    // every node of the compressed quadtree: a cell whose points meet >= 2 cells of the next split level
    std::function<void(std::vector<int>)> node = [&](std::vector<int> P) {
        int k = splitLevel(P); if (k > K) return;
        auto C = cellsAt(P, k);
        for (size_t a = 0; a < C.size(); ++a) for (size_t b = a + 1; b < C.size(); ++b) visit(C[a], C[b], k);
        for (auto& c : C) node(std::move(c));
    };
    std::vector<int> all(n); std::iota(all.begin(), all.end(), 0);
    node(all);
    return W;
}

} // namespace wspd

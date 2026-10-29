// lb_swap: enlarges the lower bound of the sweep by local search (notes.md §23). Like check_lb, this file shares
// nothing with wspd.hpp and needs neither pgl nor CPLEX:   g++ -std=c++17 -O2 -o lb_swap lb_swap.cpp
//
// usage: lb_swap POINTS ANCHORS a b OUT [PLATEAU]   (eps = a/b < 1/2; files as for check_lb, written by exp with WSPD_ANCHORS)
//
// Any set of pairwise incompatible pairs is a lower bound on the size of an eps-WSPD (see check_lb.cpp). Starting from
// the anchors, the program (1) adds every pair that is incompatible with all pairs of the set, and (2) repeats 1-for-2
// swaps: a pair e of the set is replaced by two pairs f, g that are incompatible with each other and compatible with
// no pair of the set other than e; the pairs that become free by a swap are added too. It stops when a pass over the
// set finds no swap. With PLATEAU > 0, that many more passes follow in which a pair e without a swap is replaced, with
// probability 1/2, by a random pair compatible with no other pair of the set (the size stays the same, and later
// passes may find new swaps). The result is written to OUT in the format of ANCHORS; it is certified by check_lb, not here.
// With PLATEAU > 0, the set before the plateau passes is written to OUT.swaps, so that both sets can be certified.
//
// For every pair w, t[w] is the number of pairs of the set that are compatible with w. The pairs compatible with
// {x,y} are the pairs {u,v} with {{x,u},{y,v}} eps-separated, so |xu|, |yv| <= eps|xy|: they are enumerated from the
// lists of the points sorted by distance from x and from y (4 n^2 bytes). All tests are exact, as in check_lb.
// Output: one line  n,anchors,free,swaps,swapFree,final,passes,converged,capHits,seconds   (final = size after the
// passes without plateau moves; with PLATEAU > 0, five more fields: the size at the end, and passes, swaps, swapFree,
// seconds when the passes without plateau moves ended; the first swaps and swapFree count all passes)
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

using ll = long long;
using i128 = __int128;

static std::vector<ll> X, Y;
static ll A, B;
static int n;
static std::vector<std::vector<int>> nearL;   // nearL[c]: all points sorted by distance from c (c itself first)
static std::vector<uint8_t> t;                // t[min*n+max]
static std::vector<bool> inI;

static ll d2(int i, int j) { ll dx = X[i] - X[j], dy = Y[i] - Y[j]; return dx * dx + dy * dy; }
static bool withinEps(ll duv2, ll D2) { return (i128)B * B * duv2 <= (i128)A * A * D2; }
static bool separated(int p, int p2, int q, int q2) {
    ll diam2 = std::max(d2(p, p2), d2(q, q2));
    ll dist2 = std::min(std::min(d2(p, q), d2(p, q2)), std::min(d2(p2, q), d2(p2, q2)));
    return withinEps(diam2, dist2);
}
static bool compatible(int p, int q, int p2, int q2) { return separated(p, p2, q, q2) || separated(p, q2, q, p2); }
static size_t idx(int i, int j) { return i < j ? (size_t)i * n + j : (size_t)j * n + i; }

// calls f(u, v) for every pair {u,v} != {x,y} compatible with {x,y}, with u on the side of x; each pair once (eps < 1/2)
template <class F> static void forCompatible(int x, int y, F f) {
    ll L2 = d2(x, y);
    for (int u : nearL[x]) {
        if (!withinEps(d2(x, u), L2)) break;
        for (int v : nearL[y]) {
            if (!withinEps(d2(y, v), L2)) break;
            if (u == x && v == y) continue;
            if (separated(x, u, y, v)) f(u, v);
        }
    }
}
static void add(int x, int y) {
    inI[idx(x, y)] = true;
    forCompatible(x, y, [&](int u, int v) { uint8_t& c = t[idx(u, v)]; if (c == 255) { std::puts("FAIL count overflow"); std::exit(1); } ++c; });
}
static void removePair(int x, int y) {
    inI[idx(x, y)] = false;
    forCompatible(x, y, [&](int u, int v) { --t[idx(u, v)]; });
}

int main(int argc, char** argv) {
    if (argc != 6 && argc != 7) { std::fprintf(stderr, "usage: lb_swap POINTS ANCHORS a b OUT [PLATEAU]\n"); return 2; }
    auto t0 = std::chrono::steady_clock::now();
    A = std::atoll(argv[3]); B = std::atoll(argv[4]);
    if (A <= 0 || 2 * A >= B) { std::fprintf(stderr, "need 0 < eps < 1/2\n"); return 2; }
    FILE* f = std::fopen(argv[1], "r"); if (!f) { std::perror(argv[1]); return 2; }
    for (ll x, y; std::fscanf(f, "%lld %lld", &x, &y) == 2;) { X.push_back(x); Y.push_back(y); }
    std::fclose(f);
    n = (int)X.size();
    std::vector<std::pair<int, int>> S;
    f = std::fopen(argv[2], "r"); if (!f) { std::perror(argv[2]); return 2; }
    for (int p, q; std::fscanf(f, "%d %d", &p, &q) == 2;) S.push_back({p, q});
    std::fclose(f);
    size_t anchors = S.size();

    nearL.resize(n);
    {
        std::vector<std::pair<ll, int>> key(n);
        for (int c = 0; c < n; ++c) {
            for (int z = 0; z < n; ++z) key[z] = {d2(c, z), z};
            std::sort(key.begin(), key.end());
            nearL[c].resize(n);
            for (int z = 0; z < n; ++z) nearL[c][z] = key[z].second;
        }
    }
    t.assign((size_t)n * n, 0); inI.assign((size_t)n * n, false);
    for (auto [p, q] : S) {
        if (p == q || p < 0 || q < 0 || p >= n || q >= n || inI[idx(p, q)]) { std::puts("FAIL bad anchor"); return 1; }
        add(p, q);
    }
    for (auto [p, q] : S) if (t[idx(p, q)]) { std::puts("FAIL compatible anchors"); return 1; }

    // (1) free pairs
    size_t freeAdded = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
        if (!inI[idx(i, j)] && !t[idx(i, j)]) { add(i, j); S.push_back({i, j}); ++freeAdded; }

    // (2) 1-for-2 swaps
    const size_t CAP = 200000;   // compatibility tests per pair of the set and pass
    size_t swaps = 0, swapFree = 0, capHits = 0; int passes = 0; bool converged = false;
    std::vector<std::pair<int, int>> T; std::vector<int> xs, ys;
    std::vector<char> alive(S.size(), 1);
    int plateau = argc == 7 ? std::atoi(argv[6]) : 0, plateauLeft = 0; size_t finNoPlateau = 0; std::mt19937 rng(1);
    int passes0 = 0; size_t swaps0 = 0, swapFree0 = 0; double sec0 = 0;
    auto writeSet = [&](const std::string& path) {
        FILE* g = std::fopen(path.c_str(), "w"); if (!g) { std::perror(path.c_str()); std::exit(2); }
        size_t c = 0;
        for (size_t k = 0; k < S.size(); ++k) if (alive[k]) { std::fprintf(g, "%d %d\n", S[k].first, S[k].second); ++c; }
        if (std::fclose(g)) { std::perror(path.c_str()); std::exit(2); }
        return c;
    };
    auto count = [&] { size_t c = 0; for (char a : alive) c += a; return c; };
    while ((passes < 50 && !converged) || plateauLeft > 0) {
        bool pl = converged || passes >= 50;   // a pass with plateau moves
        if (pl) --plateauLeft;
        ++passes; size_t before = swaps; size_t end = S.size();
        for (size_t k = 0; k < end; ++k) {
            if (!alive[k]) continue;
            auto [x, y] = S[k];
            T.clear();
            forCompatible(x, y, [&](int u, int v) { if (t[idx(u, v)] == 1) T.push_back({u, v}); });
            if (T.empty()) continue;
            auto move = [&] {   // plateau move: e is replaced by a random pair of T, and the pairs of T that become free are added
                if (!pl || (rng() & 1)) return;
                auto w = T[rng() % T.size()];
                removePair(x, y); alive[k] = 0; add(w.first, w.second); S.push_back(w); alive.push_back(1);
                for (auto [u, v] : T) if (!inI[idx(u, v)] && !t[idx(u, v)]) { add(u, v); S.push_back({u, v}); alive.push_back(1); ++swapFree; }
            };
            if (T.size() < 2) { move(); continue; }
            // if the dumbbell of all first and all second endpoints is eps-separated, the pairs of T are pairwise compatible
            xs.clear(); ys.clear();
            for (auto [u, v] : T) { xs.push_back(u); ys.push_back(v); }
            std::sort(xs.begin(), xs.end()); xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
            std::sort(ys.begin(), ys.end()); ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
            ll diam2 = 0, dist2 = INT64_MAX;
            for (size_t i = 0; i < xs.size(); ++i) for (size_t j = i + 1; j < xs.size(); ++j) diam2 = std::max(diam2, d2(xs[i], xs[j]));
            for (size_t i = 0; i < ys.size(); ++i) for (size_t j = i + 1; j < ys.size(); ++j) diam2 = std::max(diam2, d2(ys[i], ys[j]));
            for (int u : xs) for (int v : ys) dist2 = std::min(dist2, d2(u, v));
            if (withinEps(diam2, dist2)) { move(); continue; }
            int fi = -1, gi = -1; size_t tests = 0;
            for (size_t i = 0; i < T.size() && fi < 0 && tests < CAP; ++i)
                for (size_t j = i + 1; j < T.size() && tests < CAP; ++j) {
                    ++tests;
                    if (!compatible(T[i].first, T[i].second, T[j].first, T[j].second)) { fi = (int)i; gi = (int)j; break; }
                }
            if (fi < 0) { if (tests >= CAP) ++capHits; move(); continue; }
            removePair(x, y); alive[k] = 0;
            for (int i : {fi, gi}) { add(T[i].first, T[i].second); S.push_back(T[i]); alive.push_back(1); }
            ++swaps;
            for (auto [u, v] : T) if (!inI[idx(u, v)] && !t[idx(u, v)]) { add(u, v); S.push_back({u, v}); alive.push_back(1); ++swapFree; }
        }
        {   // drop the removed pairs (the order of the others, and so the result, does not change)
            size_t w = 0;
            for (size_t k = 0; k < S.size(); ++k) if (alive[k]) S[w++] = S[k];
            S.resize(w); alive.assign(w, 1);
        }
        if (!pl) { converged = swaps == before; if (converged || passes >= 50) {
            finNoPlateau = count(); plateauLeft = plateau; passes0 = passes; swaps0 = swaps; swapFree0 = swapFree;
            sec0 = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (plateau > 0 && writeSet(std::string(argv[5]) + ".swaps") != finNoPlateau) { std::puts("FAIL count"); return 1; }
        } }
    }
    size_t fin = writeSet(argv[5]);
    std::printf("%d,%zu,%zu,%zu,%zu,%zu,%d,%d,%zu,%.2f", n, anchors, freeAdded, swaps, swapFree, finNoPlateau, passes, (int)converged, capHits,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    if (plateau > 0) std::printf(",%zu,%d,%zu,%zu,%.2f", fin, passes0, swaps0, swapFree0, sec0);
    std::puts("");
    return 0;
}

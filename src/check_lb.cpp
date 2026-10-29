// check_lb: checks the lower bound of the sweep, independently of the code that computes it. This file shares nothing
// with wspd.hpp and needs neither pgl nor CPLEX:   g++ -std=c++17 -O2 -o check_lb check_lb.cpp
//
// usage: check_lb POINTS ANCHORS a b [brute]        (eps = a/b)
//   POINTS    one "x y" line per point (integers); the point index is the line number, from 0
//   ANCHORS   one "p q" line per anchor (point indices): the pairs that start a step of the sweep
// Both files are written by exp with WSPD_ANCHORS=PREFIX (PREFIX.pts and PREFIX.anchors).
//
// Two pairs {p,q} and {p',q'} are compatible if some eps-separated dumbbell {X,Y} covers both, that is, after possibly
// exchanging p' and q', p and p' are in X and q and q' are in Y. Subsets of the sides of an eps-separated dumbbell are
// eps-separated, so this holds if and only if {{p,p'},{q,q'}} or {{p,q'},{q,p'}} is eps-separated (compatible() below).
// If the anchors are pairwise incompatible, every dumbbell of an eps-WSPD covers at most one anchor, and every anchor
// is covered, so every eps-WSPD (cover or partition) of the points has at least as many dumbbells as there are anchors.
//
// The program checks that the anchors are distinct pairs of distinct points and that no two of them are compatible.
// All tests are exact:  max(diam X, diam Y) <= eps dist(X,Y)  <=>  b^2 diam^2 <= a^2 dist^2,  in 128-bit integers.
// Exit status 0 and a line "OK ..." if the lower bound holds, 1 and a line "FAIL ..." otherwise.
//
// Search. With "brute", every pair of anchors is tested (quadratic in the number of anchors; for small inputs and to
// cross-check the default). By default, for every anchor {p,q}, only the anchors {p',q'} with |pp'| <= eps|pq| and
// |qq'| <= eps|pq| are tested. This misses no compatible pair: if {{p,p'},{q,q'}} is eps-separated, then |pp'| and
// |qq'| are at most eps dist({p,p'},{q,q'}) <= eps|pq|. The points within eps|pq| of p and of q are found with a
// uniform grid (the cells of a square around the point that contains the disk, each point then tested exactly), and
// whether {p',q'} is an anchor is read from a table of n^2 bits (1.25 GB for n = 100000).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <utility>
#include <vector>

using ll = long long;
using i128 = __int128;

static std::vector<ll> X, Y;
static ll A, B;   // eps = A/B

static ll d2(int i, int j) { ll dx = X[i] - X[j], dy = Y[i] - Y[j]; return dx * dx + dy * dy; }

// |uv| <= eps * sqrt(D2), for the squared distances duv2 = |uv|^2 and D2
static bool withinEps(ll duv2, ll D2) { return (i128)B * B * duv2 <= (i128)A * A * D2; }

// is the dumbbell {{p,p2},{q,q2}} eps-separated?  (if the two sides share a point, the distance is 0 and it is not)
static bool separated(int p, int p2, int q, int q2) {
    ll diam2 = std::max(d2(p, p2), d2(q, q2));
    ll dist2 = std::min(std::min(d2(p, q), d2(p, q2)), std::min(d2(p2, q), d2(p2, q2)));
    return withinEps(diam2, dist2);
}

// are the two distinct pairs {p,q} and {p2,q2} covered by a common eps-separated dumbbell?
static bool compatible(int p, int q, int p2, int q2) { return separated(p, p2, q, q2) || separated(p, q2, q, p2); }

static int fail(const char* what, int p, int q, int p2 = -1, int q2 = -1) {
    if (p2 < 0) std::printf("FAIL %s: anchor {%d,%d}\n", what, p, q);
    else std::printf("FAIL %s: anchors {%d,%d} and {%d,%d}\n", what, p, q, p2, q2);
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "usage: check_lb POINTS ANCHORS a b [brute]\n"); return 2; }
    A = std::atoll(argv[3]); B = std::atoll(argv[4]);
    bool brute = argc > 5 && std::string_view(argv[5]) == "brute";
    if (A <= 0 || B <= 0 || A > 1000000000 || B > 1000000000) { std::fprintf(stderr, "eps = a/b needs 0 < a, b <= 10^9\n"); return 2; }
    FILE* f = std::fopen(argv[1], "r");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
    for (ll x, y; std::fscanf(f, "%lld %lld", &x, &y) == 2;) { X.push_back(x); Y.push_back(y); }
    std::fclose(f);
    int n = (int)X.size();
    if (n < 2) { std::fprintf(stderr, "fewer than 2 points in %s\n", argv[1]); return 2; }
    ll minx = *std::min_element(X.begin(), X.end()), maxx = *std::max_element(X.begin(), X.end());
    ll miny = *std::min_element(Y.begin(), Y.end()), maxy = *std::max_element(Y.begin(), Y.end());
    // squared distances must fit in 63 bits, and their products with a^2 and b^2 in 127 bits
    if ((i128)maxx - minx > 2000000000 || (i128)maxy - miny > 2000000000) { std::fprintf(stderr, "coordinate range above 2*10^9\n"); return 2; }
    f = std::fopen(argv[2], "r");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", argv[2]); return 2; }
    std::vector<std::pair<int,int>> anchors;
    for (ll p, q; std::fscanf(f, "%lld %lld", &p, &q) == 2;) {
        if (p < 0 || q < 0 || p >= n || q >= n) { std::printf("FAIL point index out of range: anchor {%lld,%lld}\n", p, q); return 1; }
        anchors.push_back({(int)p, (int)q});
    }
    std::fclose(f);

    // the anchors are distinct pairs of distinct points
    std::vector<uint64_t> isAnchor(((size_t)n * n + 63) / 64, 0);   // bit u*n+v, set for both orders of an anchor
    auto bit = [&](int u, int v) { size_t k = (size_t)u * n + v; return (isAnchor[k >> 6] >> (k & 63)) & 1; };
    auto setBit = [&](int u, int v) { size_t k = (size_t)u * n + v; isAnchor[k >> 6] |= (uint64_t)1 << (k & 63); };
    for (auto [p, q] : anchors) {
        if (d2(p, q) == 0) return fail("the two points of an anchor coincide", p, q);
        if (bit(p, q)) return fail("anchor listed twice", p, q);
        setBit(p, q); setBit(q, p);
    }

    // no two anchors are compatible
    if (brute) {
        for (size_t i = 0; i < anchors.size(); ++i) for (size_t j = i + 1; j < anchors.size(); ++j) {
            auto [p, q] = anchors[i]; auto [p2, q2] = anchors[j];
            if (compatible(p, q, p2, q2)) return fail("compatible anchors", p, q, p2, q2);
        }
    } else {
        // uniform grid of G x G cells of side g that contains every point
        ll span = std::max(maxx - minx, maxy - miny) + 1;   // a radius that reaches every point from every point
        ll G = std::max<ll>(1, (ll)std::sqrt((double)n)), g = span / G + 1;
        auto cellOf = [&](ll c, ll lo) { return std::min(G - 1, std::max<ll>(0, (c - lo) / g)); };
        std::vector<int> cellStart(G * G + 1, 0), cellPts(n);
        for (int i = 0; i < n; ++i) ++cellStart[cellOf(Y[i], miny) * G + cellOf(X[i], minx) + 1];
        for (ll c = 0; c < G * G; ++c) cellStart[c + 1] += cellStart[c];
        { std::vector<int> pos(cellStart.begin(), cellStart.end() - 1);
          for (int i = 0; i < n; ++i) cellPts[pos[cellOf(Y[i], miny) * G + cellOf(X[i], minx)]++] = i; }
        // the points z with |cz| <= eps sqrt(D2)
        auto disk = [&](int c, ll D2, std::vector<int>& out) {
            out.clear();
            ll r = std::min(span, (ll)(std::sqrt((long double)D2) * A / B));   // an estimate of eps sqrt(D2), made exact:
            while (r < span && (i128)B * B * r * r < (i128)A * A * D2) ++r;   // now r >= min(span, eps sqrt(D2))
            for (ll cy = cellOf(Y[c] - r, miny); cy <= cellOf(Y[c] + r, miny); ++cy)
                for (ll cx = cellOf(X[c] - r, minx); cx <= cellOf(X[c] + r, minx); ++cx)
                    for (int k = cellStart[cy * G + cx]; k < cellStart[cy * G + cx + 1]; ++k)
                        if (withinEps(d2(c, cellPts[k]), D2)) out.push_back(cellPts[k]);
        };
        std::vector<int> nearP, nearQ;
        for (auto [p, q] : anchors) {
            ll D2 = d2(p, q);
            disk(p, D2, nearP); disk(q, D2, nearQ);
            for (int p2 : nearP) for (int q2 : nearQ) {
                if (p2 == q2 || !bit(p2, q2)) continue;                       // {p2,q2} is not an anchor
                if ((p2 == p && q2 == q) || (p2 == q && q2 == p)) continue;   // it is {p,q} itself
                if (separated(p, p2, q, q2)) return fail("compatible anchors", p, q, p2, q2);
            }
        }
    }
    std::printf("OK lower bound %zu: the anchors are pairwise incompatible, n = %d, eps = %lld/%lld%s\n", anchors.size(), n, A, B,
                brute ? " (brute force)" : "");
    return 0;
}

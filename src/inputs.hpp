// Inputs of the experiments: the synthetic generators (at scale M = 10^6, integer coordinates, duplicates discarded and
// replaced) and the CG:SHOP 2019 loader. GEN is uniform, gauss, cluster, grid (perturbed grid), spiral, annulus, or
// file:PATH.
#pragma once
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include "wspd.hpp"

namespace wspd {

// CG:SHOP 2019 format: lines "index<TAB>x<TAB>y", comments start with '#'. Duplicates are removed.
inline Instance loadFile(const std::string& path, Eps eps) {
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << std::endl; exit(1); }
    std::set<std::pair<ll,ll>> seen; Instance I; I.eps = eps; std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line); ll id, x, y;
        if (!(ss >> id >> x >> y)) continue;
        if (seen.insert({x, y}).second) { I.x.push_back(x); I.y.push_back(y); }
    }
    return I;
}

inline Instance generate(const std::string& g, int n, unsigned seed, Eps eps) {
    if (g.rfind("file:", 0) == 0) return loadFile(g.substr(5), eps);
    std::mt19937_64 rng(seed);
    std::set<std::pair<ll,ll>> S;
    const double M = 1e6;
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> N(0, 1);
    std::vector<std::pair<double,double>> centers;
    if (g == "cluster") for (int k = 0; k < std::max(2, n / 20); ++k) centers.push_back({U(rng) * M, U(rng) * M});
    int side = (int)std::ceil(std::sqrt((double)n));
    int it = 0;
    while ((int)S.size() < n && it < 100 * n) {
        double x = 0, y = 0;
        if (g == "uniform") { x = U(rng) * M; y = U(rng) * M; }
        else if (g == "gauss") { x = M / 2 + N(rng) * M / 8; y = M / 2 + N(rng) * M / 8; }
        else if (g == "cluster") { auto c = centers[rng() % centers.size()]; x = c.first + N(rng) * M / 100; y = c.second + N(rng) * M / 100; }
        else if (g == "grid") { int k = (int)S.size() + it % 7; x = (k % side) * (M / side) + N(rng) * M / side * 0.02; y = (k / side) * (M / side) + N(rng) * M / side * 0.02; }
        else if (g == "spiral") { double t = U(rng) * 6; x = M / 2 + std::exp(t) * 2000 * std::cos(3 * t); y = M / 2 + std::exp(t) * 2000 * std::sin(3 * t); }
        else if (g == "annulus") { double t = U(rng) * 2 * M_PI, r = M / 4 * (1 + 0.1 * U(rng)); x = M / 2 + r * std::cos(t); y = M / 2 + r * std::sin(t); }
        else { std::cerr << "unknown generator " << g << std::endl; exit(1); }
        S.insert({std::llround(x), std::llround(y)}); ++it;
    }
    Instance I; I.eps = eps;
    for (auto& p : S) { I.x.push_back(p.first); I.y.push_back(p.second); }
    return I;
}

} // namespace wspd

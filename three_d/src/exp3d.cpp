// Experiment driver in R^3 (separate from the planar exp.cpp).
// usage: exp3d GEN n seed a b ALGO [ALGO...]      (eps = a/b)
//   GEN  in {uniform, gauss, cluster, grid, spiral, shell} or file:PATH (one "x y z" integer line per point, written by
//        prep_scans.py; a random sample of n distinct points is taken with the seed, all points if n = 0)
//   ALGO in {opt, alpha, sweep, sweepraw, ck, ckexact, ls, qt, qtexact, qtrand, qtexactrand, load}
//     sweep: the greedy sweep followed by redundancy removal (the paper's "sweep"); sweepraw: without the removal
//     ck: CK box; ckexact: CK exact; opt: CPLEX set cover over all maximal dumbbells (WSPD_OPT_TL seconds, default 1800)
//     alpha: the largest set of pairwise incompatible pairs (CPLEX; not a WSPD, printed without the verifier, "unchecked";
//            WSPD_ALPHA=PREFIX writes the pairs to PREFIX.alpha)
//     ls: local search from WSPD_LS_START = ck (CK exact, default) or sweep, for WSPD_LS_TIME seconds (default 10), with
//         WSPD_LS_MOVES = reshape (default) or lns (parameters LSConfig::lnsDefaults; WSPD_LS_K overrides k);
//         WSPD_LS_TRACE: trace file ("time,size"); WSPD_LS_ROUNDS: stop after this many rounds
//     qt*: octree partition WSPD (partition3.hpp; exact: separation tested on the points; rand: randomly shifted and scaled)
//     load: the WSPD in the file WSPD_LOAD (format of WSPD_DUMP) for the same instance; printed as load-<file extension>
// WSPD_PARTITION=SPEC[,SPEC...] (SPEC = order/split/merge, see PartConfig): after each algorithm, its output is converted
// into a partition WSPD and one more line per SPEC is printed, with algo "<algo>+part", the conversion time, and the
// peak memory of the run before the conversion (peakBeforeKB).
// Prints one CSV line per algorithm:   gen,n,seed,eps,algo,size,valid,time,info
// WSPD_ANCHORS=PREFIX: the sweep writes the points to PREFIX.pts ("x y z") and its anchors to PREFIX.anchors ("p q"),
// the input of check_lb3 and lb_swap3. WSPD_DUMP=PREFIX: every output is written to PREFIX.<algo> ("x1 x2 ...|y1 ...").
#include <fstream>
#include <sys/resource.h>
#include <set>
#include <sstream>
#include <tuple>
#include "wspd3.hpp"
#include "ls3.hpp"
#include "partition3.hpp"
using namespace wspd3;

static Instance loadFile(const std::string& path, int n, unsigned seed, Eps eps) {
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << std::endl; exit(1); }
    std::set<std::tuple<ll,ll,ll>> seen; std::vector<std::tuple<ll,ll,ll>> P; std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line); ll x, y, z;
        if (!(ss >> x >> y >> z)) continue;
        if (seen.insert({x, y, z}).second) P.push_back({x, y, z});
    }
    if (n > 0 && n < (int)P.size()) {   // the first n of a random permutation, in file order
        std::mt19937_64 rng(seed);
        std::vector<int> id(P.size()); std::iota(id.begin(), id.end(), 0);
        std::shuffle(id.begin(), id.end(), rng);
        id.resize(n); std::sort(id.begin(), id.end());
        std::vector<std::tuple<ll,ll,ll>> Q; for (int i : id) Q.push_back(P[i]);
        P.swap(Q);
    }
    Instance I; I.eps = eps;
    for (auto [x, y, z] : P) { I.x.push_back(x); I.y.push_back(y); I.z.push_back(z); }
    return I;
}

static Instance generate(const std::string& g, int n, unsigned seed, Eps eps) {
    if (g.rfind("file:", 0) == 0) return loadFile(g.substr(5), n, seed, eps);
    std::mt19937_64 rng(seed);
    std::set<std::tuple<ll,ll,ll>> S;
    const double M = 1e6;
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> N(0, 1);
    std::vector<std::array<double,3>> centers;
    if (g == "cluster") for (int k = 0; k < std::max(2, n / 20); ++k) centers.push_back({U(rng) * M, U(rng) * M, U(rng) * M});
    int side = (int)std::ceil(std::cbrt((double)n) - 1e-9);
    int it = 0;
    while ((int)S.size() < n && it < 100 * n) {
        double x = 0, y = 0, z = 0;
        if (g == "uniform") { x = U(rng) * M; y = U(rng) * M; z = U(rng) * M; }
        else if (g == "gauss") { x = M / 2 + N(rng) * M / 8; y = M / 2 + N(rng) * M / 8; z = M / 2 + N(rng) * M / 8; }
        else if (g == "cluster") { auto c = centers[rng() % centers.size()]; x = c[0] + N(rng) * M / 100; y = c[1] + N(rng) * M / 100; z = c[2] + N(rng) * M / 100; }
        else if (g == "grid") {   // perturbed grid
            int k = (int)S.size() + it % 7; double h = M / side;
            x = (k % side) * h + N(rng) * h * 0.02; y = (k / side % side) * h + N(rng) * h * 0.02; z = (k / side / side) * h + N(rng) * h * 0.02;
        }
        else if (g == "spiral") {  // logarithmic spiral on a cone (distances between neighbours over 2.5 orders of magnitude)
            double t = U(rng) * 6, r = std::exp(t) * 2000;
            x = M / 2 + r * std::cos(3 * t); y = M / 2 + r * std::sin(3 * t); z = M / 2 + r;
        }
        else if (g == "shell") {   // spherical shell of radii M/4 and 1.1 M/4
            double a = N(rng), b = N(rng), c = N(rng), l = std::sqrt(a * a + b * b + c * c), r = M / 4 * (1 + 0.1 * U(rng));
            if (l == 0) { ++it; continue; }
            x = M / 2 + r * a / l; y = M / 2 + r * b / l; z = M / 2 + r * c / l;
        }
        else { std::cerr << "unknown generator " << g << std::endl; exit(1); }
        S.insert({std::llround(x), std::llround(y), std::llround(z)}); ++it;
    }
    Instance I; I.eps = eps;
    for (auto [x, y, z] : S) { I.x.push_back(x); I.y.push_back(y); I.z.push_back(z); }
    return I;
}

int main(int argc, char** argv) {
    if (argc < 7) { std::cerr << "usage: exp3d GEN n seed a b ALGO..." << std::endl; return 1; }
    std::string gen = argv[1]; int n = std::atoi(argv[2]); unsigned seed = std::atoi(argv[3]);
    Eps eps{std::atoll(argv[4]), std::atoll(argv[5])};
    Instance I = generate(gen, n, seed, eps);
    for (const char* var : {"WSPD_DUMP", "WSPD_ANCHORS"}) if (const char* dump = std::getenv(var)) {
        std::ofstream f(std::string(dump) + ".pts");
        for (int i = 0; i < I.n(); ++i) f << I.x[i] << " " << I.y[i] << " " << I.z[i] << "\n";
    }
    const char* anchorDump = std::getenv("WSPD_ANCHORS");
    std::vector<std::pair<int,int>> anchors; auto anchorsOut = anchorDump ? &anchors : nullptr;
    auto env = [](const char* k, const char* d) { const char* v = std::getenv(k); return std::string(v ? v : d); };
    for (int k = 6; k < argc; ++k) {
        std::string algo = argv[k];
        std::vector<Pair> W; std::string info;
        auto t0 = std::chrono::steady_clock::now();
        if (algo == "alpha") {   // not a WSPD: printed here, without the verifier
            Spatial S(I); bool capHit = false;
            auto cols = allMaximalPairs(I, S, &capHit);
            auto pr = maxIncompatibleCplex(I.n(), cols, std::atof(env("WSPD_OPT_TL", "1800").c_str()));
            std::string gname = gen.rfind("file:", 0) == 0 ? gen.substr(gen.find_last_of('/') + 1, gen.find_last_of('.') - gen.find_last_of('/') - 1) : gen;
            std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << ",alpha," << pr.chosen.size() << ",unchecked,"
                      << seconds_since(t0) << ",cols=" << cols.size() << ";optimal=" << (pr.optimal && !capHit) << ";ub=" << pr.ub << std::endl;
            if (const char* dump = std::getenv("WSPD_ALPHA")) {
                std::ofstream f(std::string(dump) + ".alpha");
                for (auto [p, q] : pr.chosen) f << p << " " << q << "\n";
            }
            continue;
        } else if (algo == "opt") {
            Spatial S(I); bool capHit = false;
            auto cols = allMaximalPairs(I, S, &capHit);
            double tEnum = seconds_since(t0);
            std::vector<std::pair<int,int>> E;
            for (int i = 0; i < I.n(); ++i) for (int j = i + 1; j < I.n(); ++j) E.push_back({i, j});
            auto cr = setCoverCplex(I.n(), E, cols, std::atof(env("WSPD_OPT_TL", "1800").c_str()));
            for (int c : cr.chosen) W.push_back(cols[c]);
            info = "cols=" + std::to_string(cols.size()) + ";optimal=" + std::to_string(cr.optimal && !capHit) + ";lb=" + std::to_string(cr.lb)
                 + ";capHit=" + std::to_string(capHit) + ";tEnum=" + std::to_string(tEnum);
        } else if (algo == "sweep" || algo == "sweepraw") {
            SweepStats st; W = sweep(I, st, 16, anchorsOut);
            double tSweep = seconds_since(t0);
            size_t raw = W.size();
            if (algo == "sweep") W = pruneRedundant(I, std::move(W));
            info = "steps=" + std::to_string(st.steps) + ";raw=" + std::to_string(raw) + ";tSweep=" + std::to_string(tSweep)
                 + ";avgR=" + std::to_string(double(st.rTotal) / st.steps) + ";fixedChosen=" + std::to_string(st.fixedChosen) + ";capHits=" + std::to_string(st.capHits);
        } else if (algo == "ck" || algo == "ckexact") {
            W = callahanKosaraju(I, algo == "ckexact");
        } else if (algo == "ls") {
            std::string start = env("WSPD_LS_START", "ck");
            LSConfig cf; cf.timeLimit = std::atof(env("WSPD_LS_TIME", "10").c_str()); cf.traceFile = env("WSPD_LS_TRACE", "");
            cf.seed = seed + 777; cf.maxRounds = std::atol(env("WSPD_LS_ROUNDS", "0").c_str());
            cf.moves = env("WSPD_LS_MOVES", "reshape"); cf.lnsDefaults(I.n());
            cf.lnsK = std::atoi(env("WSPD_LS_K", std::to_string(cf.lnsK).c_str()).c_str());
            std::vector<Pair> W0;
            if (start == "sweep") { SweepStats st; W0 = pruneRedundant(I, sweep(I, st)); }
            else if (start == "ck") W0 = callahanKosaraju(I, true);
            else { std::cerr << "unknown start " << start << std::endl; return 1; }
            double tBuild = seconds_since(t0);
            size_t startSize = W0.size();
            LSResult R; W = runLocalSearch(I, std::move(W0), cf, R);
            auto tPr = std::chrono::steady_clock::now();
            W = pruneRedundant(I, std::move(W));
            double tPrune = seconds_since(tPr);
            algo = "ls-" + start + "-" + cf.moves;
            info = "start=" + std::to_string(startSize) + ";tBuild=" + std::to_string(tBuild) + ";lsSetup=" + std::to_string(R.setup) + ";tPrune=" + std::to_string(tPrune)
                 + ";lsTime=" + std::to_string(R.time) + ";rounds=" + std::to_string(R.rounds) + ";reshapeMoves=" + std::to_string(R.reshapeMoves)
                 + ";reshapeKills=" + std::to_string(R.reshapeKills) + ";extraKills=" + std::to_string(R.extraKills);
            if (cf.moves == "lns") info += ";lnsCalls=" + std::to_string(R.lnsCalls) + ";lnsGain=" + std::to_string(R.lnsGain)
                                         + ";lnsBig=" + std::to_string(R.lnsBig) + ";lnsFail=" + std::to_string(R.lnsFail);
        } else if (algo == "qt" || algo == "qtexact" || algo == "qtrand" || algo == "qtexactrand") {
            QTStats qs; bool rnd = algo.find("rand") != std::string::npos;
            W = quadtreePartition(I, algo.find("exact") != std::string::npos, rnd, seed + 4242, qs);
            std::string why; bool part = verifyPartition(I, W, &why);
            info = "partition=" + std::string(part ? "1" : "0:" + why) + ";calls=" + std::to_string(qs.calls) + ";levels=" + std::to_string(qs.levels);
        } else if (algo == "load") {
            const char* lp = std::getenv("WSPD_LOAD"); if (!lp) { std::cerr << "WSPD_LOAD not set" << std::endl; return 1; }
            std::string path = lp; std::ifstream f(path); std::string line;
            if (!f) { std::cerr << "cannot open WSPD_LOAD" << std::endl; return 1; }
            algo = "load-" + path.substr(path.find_last_of('.') + 1);
            while (std::getline(f, line)) {
                auto bar = line.find('|'); if (bar == std::string::npos) continue;
                Pair p; std::istringstream a(line.substr(0, bar)), b(line.substr(bar + 1)); int v;
                while (a >> v) p.X.push_back(v);
                while (b >> v) p.Y.push_back(v);
                W.push_back(std::move(p));
            }
        } else { std::cerr << "unknown algo " << algo << std::endl; return 1; }
        double t = seconds_since(t0);
        std::string why; bool ok = verify(I, W, &why);
        std::string gname = gen.rfind("file:", 0) == 0 ? gen.substr(gen.find_last_of('/') + 1, gen.find_last_of('.') - gen.find_last_of('/') - 1) : gen;
        std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << "," << algo << "," << W.size() << ","
                  << (ok ? "ok" : "FAIL:" + why) << "," << t << "," << info << std::endl;
        if (anchorDump && !anchors.empty()) {
            std::ofstream f(std::string(anchorDump) + ".anchors");
            for (auto [p, q] : anchors) f << p << " " << q << "\n";
            anchors.clear();
        }
        if (const char* dump = std::getenv("WSPD_DUMP")) {
            std::ofstream f(std::string(dump) + "." + algo);
            for (auto p : W) { p = canonical(p);
                for (size_t i = 0; i < p.X.size(); ++i) { f << (i ? " " : "") << p.X[i]; } f << "|";
                for (size_t i = 0; i < p.Y.size(); ++i) { f << (i ? " " : "") << p.Y[i]; } f << "\n"; }
        }
        if (const char* spec = std::getenv("WSPD_PARTITION")) {
            std::vector<std::string> specs; std::string cur;
            for (const char* c = spec; ; ++c) { if (*c == ',' || !*c) { if (!cur.empty()) specs.push_back(cur); cur.clear(); if (!*c) break; } else cur += *c; }
            // peak memory of the run before the conversion (the conversion may need more): printed in the +part row, so that
            // the memory of the algorithm itself is known when it is converted in the same run
            struct rusage ru; getrusage(RUSAGE_SELF, &ru); long peakBefore = ru.ru_maxrss;
            for (size_t sp = 0; sp < specs.size(); ++sp) {
                PartConfig pc = PartConfig::parse(specs[sp]); pc.seed = seed + 99; PartStats ps;
                auto Wp = toPartition(I, sp + 1 == specs.size() ? std::move(W) : W, pc, ps);
                std::string why2; bool pok = verifyPartition(I, Wp, &why2);
                std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << "," << algo << "+part," << Wp.size() << ","
                          << (pok ? "ok" : "FAIL:" + why2) << "," << ps.time << ",conv=" << pc.name() << ";cover=" << ps.coverSize
                          << ";overlap=" << ps.overlap << ";dropped=" << ps.dropped << ";split=" << ps.split << ";maxPieces=" << ps.maxPieces
                          << ";beforeMerge=" << ps.beforeMerge << ";merges=" << ps.merges << ";tMerge=" << ps.tMerge << ";peakBeforeKB=" << peakBefore << std::endl;
            }
        }
    }
    return 0;
}

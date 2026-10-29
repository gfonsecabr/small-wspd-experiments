// Experiment driver.
// usage: exp GEN n seed a b ALGO [ALGO...]      (eps = a/b; GEN: see inputs.hpp; for file:PATH, n and seed are ignored)
//   ALGO: lazyprune (the sweep: greedy sweep + redundancy pass), lazy (the greedy sweep alone), ck, ckexact (Callahan-Kosaraju with the box or the exact
//         separation test), opt (exact OPT), ls (local search), qtexact, qtexactrand (the quadtree partition WSPD,
//         deterministic or randomly shifted and scaled), load (the WSPD in the file WSPD_LOAD, written by WSPD_DUMP
//         for the same instance; printed as load-<extension of the file>), alpha (not a WSPD: the largest set of pairwise
//         incompatible pairs, a lower bound on OPT, by CPLEX; size = its number of pairs, valid = unchecked; check it with
//         check_lb on the file of WSPD_ALPHA)
// Prints one CSV line per algorithm:  gen,n,seed,eps,algo,size,valid,time,info
// Environment:
//   WSPD_OPT_TL      CPLEX time limit of opt in seconds (default 1800)
//   WSPD_LS_START    start of ls: ck (CK exact, default) | sweep
//   WSPD_LS_MOVES    move of ls: reshape (default) | lns
//   WSPD_LS_TIME     time budget of ls in seconds (default 10)
//   WSPD_LS_TRACE    trace file of ls ("time,size" lines)
//   WSPD_LS_ROUNDS   stop ls after this many rounds (default 0: no limit); with the rounds of a timed run, the run ends
//                    with the same WSPD
//   WSPD_LS_K, WSPD_LS_LNS_RADIUS, WSPD_LS_LNS_UMAX   parameters of LNS: neighbours removed with the target, radius of
//                    the neighbourhood, and largest number of uncovered pairs for which new candidates are generated
//                    (defaults: LSConfig::lnsDefaults in ls.hpp, which depend on n)
//   WSPD_LS_LNS_TL, WSPD_LS_LNS_CAP   CPLEX time limit per LNS call (default 2 s) and number of candidate dumbbells
//                    enumerated per uncovered pair (default 300); used only in the tests of the LNS parameters
//   WSPD_DUMP        PREFIX: write the points to PREFIX.pts and every WSPD to PREFIX.<algo>
//   WSPD_ANCHORS     PREFIX: write the points to PREFIX.pts and the anchors of the steps of lazyprune, whose number is
//                    the lower bound, to PREFIX.anchors (one "p q" line per step), the input of check_lb
//   WSPD_ALPHA       PREFIX: write the points to PREFIX.pts and the pairs of alpha to PREFIX.alpha (one "p q" line per pair)
//   WSPD_PARTITION   if set: after each algorithm, convert its output into a partition WSPD (partition.hpp) and print
//                    one more line, with algo "<algo>+part" and the conversion time
#include <fstream>
#include <sstream>
#include "wspd.hpp"
#include "ls.hpp"
#include "partition.hpp"
#include "inputs.hpp"
using namespace wspd;

static void writeWSPD(const std::string& path, const std::vector<Pair>& W) {   // one dumbbell per line: "x1 x2 ...|y1 y2 ..."
    std::ofstream f(path);
    for (auto p : W) {
        p = canonical(p);
        for (size_t k = 0; k < p.X.size(); ++k) f << (k ? " " : "") << p.X[k];
        f << "|";
        for (size_t k = 0; k < p.Y.size(); ++k) f << (k ? " " : "") << p.Y[k];
        f << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 7) { std::cerr << "usage: exp GEN n seed a b ALGO..." << std::endl; return 1; }
    std::string gen = argv[1]; int n = std::atoi(argv[2]); unsigned seed = std::atoi(argv[3]);
    Eps eps{std::atoll(argv[4]), std::atoll(argv[5])};
    Instance I = generate(gen, n, seed, eps);
    auto env = [](const char* k, const char* d) { const char* v = std::getenv(k); return std::string(v ? v : d); };
    for (const char* var : {"WSPD_DUMP", "WSPD_ANCHORS", "WSPD_ALPHA"}) if (const char* dump = std::getenv(var)) {
        std::ofstream f(std::string(dump) + ".pts");     // point coordinates, one "x y" line per point index
        for (int i = 0; i < I.n(); ++i) f << I.x[i] << " " << I.y[i] << "\n";
    }
    const char* anchorDump = std::getenv("WSPD_ANCHORS");
    std::string gname = gen.rfind("file:", 0) == 0 ? gen.substr(gen.find_last_of('/') + 1) : gen;
    for (int k = 6; k < argc; ++k) {
        std::string algo = argv[k];
        std::vector<Pair> W; std::string info;
        std::vector<std::pair<int,int>> anchors;   // of the sweep (lazyprune), if WSPD_ANCHORS is set
        auto t0 = std::chrono::steady_clock::now();
        if (algo == "alpha") {   // not a WSPD: printed here, without the verifier
            Spatial S(I);
            auto cols = allMaximalPairs(I, S);
            auto pr = maxIncompatibleCplex(I.n(), cols, std::atof(env("WSPD_OPT_TL", "1800").c_str()), 1);
            std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << ",alpha," << pr.chosen.size() << ",unchecked,"
                      << seconds_since(t0) << ",cols=" << cols.size() << ";optimal=" << pr.optimal << ";ub=" << pr.ub << std::endl;
            if (const char* dump = std::getenv("WSPD_ALPHA")) {
                std::ofstream f(std::string(dump) + ".alpha");
                for (auto [p, q] : pr.chosen) f << p << " " << q << "\n";
            }
            continue;
        } else if (algo == "opt") {
            Spatial S(I);
            auto cols = allMaximalPairs(I, S);
            std::vector<std::pair<int,int>> E;
            for (int i = 0; i < I.n(); ++i) for (int j = i + 1; j < I.n(); ++j) E.push_back({i, j});
            auto cr = setCoverCplex(I.n(), E, cols, std::atof(env("WSPD_OPT_TL", "1800").c_str()), 1);
            for (int c : cr.chosen) W.push_back(cols[c]);
            info = "cols=" + std::to_string(cols.size()) + ";optimal=" + std::to_string(cr.optimal) + ";lb=" + std::to_string(cr.lb);
        } else if (algo == "lazyprune" || algo == "lazy") {
            SweepStats st; W = sweep(I, st, anchorDump ? &anchors : nullptr);
            if (algo == "lazyprune") W = pruneRedundant(I, std::move(W));
            info = "steps=" + std::to_string(st.steps) + ";tCand=" + std::to_string(st.tCand) + ";tCover=" + std::to_string(st.tCover)
                 + ";avgR=" + std::to_string(double(st.rTotal)/st.steps) + ";avgCols=" + std::to_string(double(st.candTotal)/st.steps)
                 + ";fixedChosen=" + std::to_string(st.fixedChosen) + ";capHits=" + std::to_string(st.capHits);
        } else if (algo == "ls") {
            std::string start = env("WSPD_LS_START", "ck"), moves = env("WSPD_LS_MOVES", "reshape");
            if ((start != "ck" && start != "sweep") || (moves != "reshape" && moves != "lns")) {
                std::cerr << "WSPD_LS_START must be ck or sweep, WSPD_LS_MOVES reshape or lns" << std::endl; return 1;
            }
            LSConfig cf; cf.lns = moves == "lns"; cf.timeLimit = std::atof(env("WSPD_LS_TIME", "10").c_str());
            cf.traceFile = env("WSPD_LS_TRACE", ""); cf.seed = seed + 777; cf.maxRounds = std::atol(env("WSPD_LS_ROUNDS", "0").c_str());
            LSConfig def; def.lnsDefaults(I.n());   // the parameters of LNS: defaults that depend on n, or the environment
            cf.lnsK = std::atoi(env("WSPD_LS_K", std::to_string(def.lnsK).c_str()).c_str());
            cf.radius = std::atof(env("WSPD_LS_LNS_RADIUS", std::to_string(def.radius).c_str()).c_str());
            cf.lnsUmax = (size_t)std::atol(env("WSPD_LS_LNS_UMAX", std::to_string(def.lnsUmax).c_str()).c_str());
            cf.lnsTL = std::atof(env("WSPD_LS_LNS_TL", std::to_string(def.lnsTL).c_str()).c_str());
            cf.lnsCap = (size_t)std::atol(env("WSPD_LS_LNS_CAP", std::to_string(def.lnsCap).c_str()).c_str());
            std::vector<Pair> W0; auto tb = std::chrono::steady_clock::now();
            if (start == "sweep") { SweepStats st; W0 = pruneRedundant(I, sweep(I, st)); }
            else W0 = callahanKosaraju(I, true);
            double tBuild = seconds_since(tb);
            size_t startSize = W0.size();   // W0 is moved into the local search (one copy of the WSPD in memory)
            LSResult R; W = runLocalSearch(I, std::move(W0), cf, R);
            auto tPr = std::chrono::steady_clock::now();
            W = pruneRedundant(I, std::move(W));
            double tPrune = seconds_since(tPr);
            algo = "ls-" + start + "-" + moves;
            info = "start=" + std::to_string(startSize) + ";tBuild=" + std::to_string(tBuild) + ";lsSetup=" + std::to_string(R.setup) + ";lsFinish=" + std::to_string(R.finish) + ";tPrune=" + std::to_string(tPrune) + ";lsTime=" + std::to_string(R.time)
                 + ";rounds=" + std::to_string(R.rounds) + ";reshapeMoves=" + std::to_string(R.reshapeMoves) + ";reshapeKills=" + std::to_string(R.reshapeKills)
                 + ";lnsCalls=" + std::to_string(R.lnsCalls) + ";lnsGain=" + std::to_string(R.lnsGain) + ";extraKills=" + std::to_string(R.extraKills);
        } else if (algo == "qtexact" || algo == "qtexactrand") {
            QTStats qs;
            W = quadtreePartition(I, algo == "qtexactrand", seed + 4242, qs);
            std::string why; bool part = verifyPartition(I, W, &why);
            info = "partition=" + std::string(part ? "1" : "0:" + why) + ";calls=" + std::to_string(qs.calls) + ";levels=" + std::to_string(qs.levels)
                 + ";S_over_L=" + std::to_string(qs.S / qs.L);
        } else if (algo == "load") {
            std::string path = env("WSPD_LOAD", ""); std::ifstream f(path); std::string line;
            if (!f) { std::cerr << "cannot open WSPD_LOAD" << std::endl; return 1; }
            algo = "load-" + path.substr(path.find_last_of('.') + 1);   // e.g. load-opt, load-ls-sweep-reshape
            while (std::getline(f, line)) {
                auto bar = line.find('|'); if (bar == std::string::npos) continue;
                Pair p; std::istringstream a(line.substr(0, bar)), b(line.substr(bar + 1)); int v;
                while (a >> v) p.X.push_back(v);
                while (b >> v) p.Y.push_back(v);
                W.push_back(std::move(p));
            }
        } else if (algo == "ck" || algo == "ckexact") {
            W = callahanKosaraju(I, algo == "ckexact");
        } else { std::cerr << "unknown algo " << algo << std::endl; return 1; }
        double t = seconds_since(t0);
        std::string why; bool ok = verify(I, W, &why);
        std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << "," << algo << "," << W.size() << ","
                  << (ok ? "ok" : "FAIL:" + why) << "," << t << "," << info << std::endl;
        if (const char* dump = std::getenv("WSPD_DUMP")) writeWSPD(std::string(dump) + "." + algo, W);
        if (anchorDump && algo == "lazyprune") {
            std::ofstream f(std::string(anchorDump) + ".anchors");
            for (auto [p, q] : anchors) f << p << " " << q << "\n";
        }
        if (std::getenv("WSPD_PARTITION")) {
            PartStats ps;
            auto Wp = toPartition(I, std::move(W), ps);
            std::string why; bool pok = verifyPartition(I, Wp, &why);
            std::cout << gname << "," << I.n() << "," << seed << "," << eps.val() << "," << algo << "+part," << Wp.size() << ","
                      << (pok ? "ok" : "FAIL:" + why) << "," << ps.time << ",conv=fewest/min/merge;cover=" << ps.coverSize
                      << ";overlap=" << ps.overlap << ";dropped=" << ps.dropped << ";split=" << ps.split << ";maxPieces=" << ps.maxPieces
                      << ";beforeMerge=" << ps.beforeMerge << ";merges=" << ps.merges << ";tMerge=" << ps.tMerge << std::endl;
        }
    }
    return 0;
}

#include "network.hpp"
#include "sim.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace crn;

namespace {

void usage() {
    std::cout <<
R"(crnsim -- stochastic simulation of chemical reaction networks

Usage:
  crnsim [--network NAME|FILE] [options]

Network:
  --network NAME|FILE   built-in name, or a path to a .crn file (default: lotka-volterra)
  --list                list built-in networks and exit
  --t-end T             simulation horizon (default: per-network)
  --seed N              RNG seed (default: 1)

Output:
  --out FILE            write the SSA trajectory as CSV ('-' for stdout)
  --sample-dt DT        record on a uniform grid instead of at every reaction event

Examples:
  crnsim --network lotka-volterra --out data/lv.csv
  crnsim --network approximate-majority --seed 3 --out data/am.csv
  crnsim --network networks/brusselator.crn --sample-dt 0.01 --out data/br.csv
)";
}

[[noreturn]] void die(const std::string& msg) {
    std::cerr << "crnsim: " << msg << "\n";
    std::exit(1);
}

std::ostream* open_out(const std::string& path, std::ofstream& file) {
    if (path == "-") return &std::cout;
    file.open(path);
    if (!file) die("cannot write '" + path + "'");
    return &file;
}

void write_trajectory(const std::string& path, const Network& net, const Trajectory& tr) {
    std::ofstream file;
    std::ostream& os = *open_out(path, file);

    os << "t";
    for (const auto& s : net.species) os << "," << s;
    os << "\n";

    char buf[64];
    for (std::size_t i = 0; i < tr.rows(); ++i) {
        std::snprintf(buf, sizeof buf, "%.6f", tr.t[i]);
        os << buf;
        const Count* row = tr.row(i);
        for (std::size_t s = 0; s < net.n_species(); ++s) os << "," << row[s];
        os << "\n";
    }
}

void describe(const Network& net, double t_end) {
    std::cerr << "network: " << net.name << "  (t_end " << t_end << ")\n";
    std::cerr << "  species:";
    for (std::size_t s = 0; s < net.n_species(); ++s)
        std::cerr << " " << net.species[s] << "=" << net.x0[s];
    std::cerr << "\n";
    for (const auto& r : net.reactions)
        std::cerr << "  " << r.label << "   k = " << r.k << "\n";
}

std::string need(int argc, char** argv, int& i) {
    if (i + 1 >= argc) die(std::string("missing value for ") + argv[i]);
    return argv[++i];
}

} // namespace

int main(int argc, char** argv) try {
    std::string spec = "lotka-volterra";
    std::string out = "-";
    double t_end = -1.0, sample_dt = 0.0;
    std::uint64_t seed = 1;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") { usage(); return 0; }
        else if (a == "--list") {
            for (const auto& b : builtin_names()) std::cout << b << "\n";
            return 0;
        }
        else if (a == "--network") spec = need(argc, argv, i);
        else if (a == "--out") out = need(argc, argv, i);
        else if (a == "--t-end") t_end = std::stod(need(argc, argv, i));
        else if (a == "--sample-dt") sample_dt = std::stod(need(argc, argv, i));
        else if (a == "--seed") seed = std::stoull(need(argc, argv, i));
        else die("unknown option '" + a + "' (try --help)");
    }

    Network net = load_network(spec);
    if (t_end < 0.0) t_end = net.suggested_t_end;

    describe(net, t_end);

    const SsaResult r = run_ssa(net, t_end, seed, sample_dt);
    std::cerr << "ssa: " << r.events << " reaction events, " << r.traj.rows() << " rows";
    if (!std::isnan(r.quiescent_time))
        std::cerr << "; absorbing state reached at t = " << r.quiescent_time;
    std::cerr << "\n";
    write_trajectory(out, net, r.traj);

    return 0;
} catch (const std::exception& e) {
    std::cerr << "crnsim: " << e.what() << "\n";
    return 1;
}

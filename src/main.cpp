#include "network.hpp"
#include "sim.hpp"

#include <algorithm>
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
  --ode-out FILE        write the RK4 deterministic trajectory as CSV
  --ode-step H          RK4 step size (default: min(0.001, ode-sample-dt/10))
  --ode-sample-dt DT    RK4 output interval (default: t_end/1000)

Examples:
  crnsim --network lotka-volterra --out data/lv.csv --ode-out data/lv_ode.csv
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
        const double* row = tr.row(i);
        for (std::size_t s = 0; s < net.n_species(); ++s) {
            if (tr.integral) {
                std::snprintf(buf, sizeof buf, "%lld", static_cast<long long>(row[s]));
            } else {
                // 10 significant digits: enough that the CSV, not the
                // integrator, is never the limit on accuracy.
                std::snprintf(buf, sizeof buf, "%.10g", row[s]);
            }
            os << "," << buf;
        }
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
    std::string out, ode_out;
    double t_end = -1.0, sample_dt = 0.0;
    double ode_step = -1.0, ode_sample_dt = -1.0;
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
        else if (a == "--ode-out") ode_out = need(argc, argv, i);
        else if (a == "--t-end") t_end = std::stod(need(argc, argv, i));
        else if (a == "--sample-dt") sample_dt = std::stod(need(argc, argv, i));
        else if (a == "--ode-step") ode_step = std::stod(need(argc, argv, i));
        else if (a == "--ode-sample-dt") ode_sample_dt = std::stod(need(argc, argv, i));
        else if (a == "--seed") seed = std::stoull(need(argc, argv, i));
        else die("unknown option '" + a + "' (try --help)");
    }

    Network net = load_network(spec);
    if (t_end < 0.0) t_end = net.suggested_t_end;
    // Output intervals scale with the horizon, so a network that finishes in
    // 0.2 time units gets as many points as one that runs for 30.
    if (ode_sample_dt < 0.0) ode_sample_dt = t_end / 1000.0;
    if (ode_step < 0.0) ode_step = std::min(0.001, ode_sample_dt / 10.0);

    // Nothing asked for: a single SSA run to stdout is the useful default.
    if (out.empty() && ode_out.empty()) out = "-";

    describe(net, t_end);

    if (!out.empty()) {
        const SsaResult r = run_ssa(net, t_end, seed, sample_dt);
        std::cerr << "ssa: " << r.events << " reaction events, " << r.traj.rows() << " rows";
        if (!std::isnan(r.quiescent_time))
            std::cerr << "; absorbing state reached at t = " << r.quiescent_time;
        std::cerr << "\n";
        write_trajectory(out, net, r.traj);
    }

    if (!ode_out.empty()) {
        const Trajectory tr = run_rk4(net, t_end, ode_step, ode_sample_dt);
        std::cerr << "ode: " << tr.rows() << " rows (rk4, h = " << ode_step << ")\n";
        write_trajectory(ode_out, net, tr);
    }

    return 0;
} catch (const std::exception& e) {
    std::cerr << "crnsim: " << e.what() << "\n";
    return 1;
}

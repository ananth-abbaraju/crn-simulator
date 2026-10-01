#pragma once

#include "crn.hpp"

#include <cstdint>
#include <vector>

namespace crn {

// A trajectory stored as a time column plus a row-major state matrix. The SSA
// and the ODE integrator share the type; `integral` records which one produced
// it, so the CSV writer can print exact counts for the SSA and real-valued
// concentrations for the ODE.
struct Trajectory {
    std::vector<double> t;
    std::vector<double> x; // t.size() * n_species, row-major
    std::size_t n_species = 0;
    bool integral = false;

    const double* row(std::size_t i) const { return &x[i * n_species]; }
    std::size_t rows() const { return t.size(); }
};

struct SsaResult {
    Trajectory traj;                // empty when record == false
    std::vector<Count> final_state;
    std::vector<double> first_zero; // per species; NaN if never hit zero
    double quiescent_time = 0.0;    // time total propensity hit 0; NaN if it never did
    unsigned long long events = 0;
};

// Exact stochastic simulation, Gillespie's direct method.
//
// sample_dt > 0 records the state on a uniform grid (bounded memory, and what
// ensemble statistics need); sample_dt == 0 records every reaction event.
SsaResult run_ssa(const Network& net, double t_end, std::uint64_t seed,
                  double sample_dt, bool record = true);

// Fourth-order Runge-Kutta on the mass-action ODEs for the same network.
// dx/dt = sum_j nu_j * flux_j(x). Integrated with step `h`, recorded every
// `sample_dt`.
Trajectory run_rk4(const Network& net, double t_end, double h, double sample_dt);

struct EnsembleStats {
    std::size_t n_runs = 0;
    std::size_t n_species = 0;

    // Per-run summaries, one entry per trajectory.
    std::vector<std::uint64_t> seeds;
    std::vector<double> quiescent_time;  // NaN when still reacting at t_end
    std::vector<double> first_zero;      // n_runs * n_species, row-major
    std::vector<Count> final_state;      // n_runs * n_species, row-major
    std::vector<unsigned long long> events;

    // Grid statistics over the ensemble.
    double grid_dt = 0.0;
    std::vector<double> grid_t;
    std::vector<double> mean;            // grid_t.size() * n_species
    std::vector<double> sd;              // grid_t.size() * n_species

    double wall_seconds = 0.0;
    unsigned threads = 1;
};

// Independent trajectories across threads. Each run owns its RNG stream, so
// results do not depend on the thread count or on scheduling.
EnsembleStats run_ensemble(const Network& net, double t_end, std::uint64_t seed,
                           std::size_t n_runs, unsigned n_threads, double grid_dt);

} // namespace crn

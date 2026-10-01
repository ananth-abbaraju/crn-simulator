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
    Trajectory traj;
    std::vector<Count> final_state;
    double quiescent_time = 0.0; // time total propensity hit 0; NaN if it never did
    unsigned long long events = 0;
};

// Exact stochastic simulation, Gillespie's direct method.
//
// sample_dt > 0 records the state on a uniform grid (bounded memory);
// sample_dt == 0 records every reaction event.
SsaResult run_ssa(const Network& net, double t_end, std::uint64_t seed,
                  double sample_dt);

// Fourth-order Runge-Kutta on the mass-action ODEs for the same network.
// dx/dt = sum_j nu_j * flux_j(x). Integrated with step `h`, recorded every
// `sample_dt`.
Trajectory run_rk4(const Network& net, double t_end, double h, double sample_dt);

} // namespace crn

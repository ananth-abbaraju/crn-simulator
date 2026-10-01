#include "sim.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace crn {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

void push_row(Trajectory& tr, double t, const std::vector<Count>& x) {
    tr.t.push_back(t);
    for (Count c : x) tr.x.push_back(static_cast<double>(c));
}

} // namespace

SsaResult run_ssa(const Network& net, double t_end, std::uint64_t seed,
                  double sample_dt) {
    const std::size_t S = net.n_species();
    const std::size_t R = net.reactions.size();

    SsaResult res;
    res.quiescent_time = kNaN;
    res.traj.n_species = S;
    res.traj.integral = true;

    std::vector<Count> x = net.x0;
    std::vector<double> a(R, 0.0);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    double t = 0.0;
    // Index of the next grid point to emit; grid point i is at i * sample_dt.
    std::size_t next_grid = 0;
    push_row(res.traj, 0.0, x);
    if (sample_dt > 0.0) next_grid = 1;

    while (t < t_end) {
        double a0 = 0.0;
        for (std::size_t j = 0; j < R; ++j) {
            a[j] = propensity(net.reactions[j], x);
            a0 += a[j];
        }

        // Absorbing state: nothing can fire. Without this guard the next line
        // divides by zero and the loop spins forever.
        if (a0 <= 0.0) {
            res.quiescent_time = t;
            break;
        }

        const double tau = -std::log(1.0 - uni(rng)) / a0;
        const double t_next = t + tau;

        // The state is constant on [t, t_next), so all grid points in that
        // half-open interval take the current counts.
        if (sample_dt > 0.0) {
            const double limit = std::min(t_next, t_end);
            while (next_grid * sample_dt < limit) {
                push_row(res.traj, next_grid * sample_dt, x);
                ++next_grid;
            }
        }

        if (t_next > t_end) {
            t = t_end;
            break;
        }
        t = t_next;

        // Pick the reaction: linear scan over the propensity array. With a
        // handful of reactions this beats any indexed structure.
        double target = uni(rng) * a0;
        std::size_t j = 0;
        for (; j + 1 < R; ++j) {
            target -= a[j];
            if (target < 0.0) break;
        }

        const Reaction& r = net.reactions[j];
        for (std::size_t s = 0; s < S; ++s)
            if (r.net[s] != 0) x[s] += r.net[s];
        ++res.events;

        if (sample_dt == 0.0) push_row(res.traj, t, x);
    }

    // Pad the grid out to t_end. Both ways out of the loop above leave the
    // counts frozen -- the horizon, or an absorbing state that can never fire
    // again -- so the padding is the true trajectory, not an extrapolation.
    if (sample_dt > 0.0) {
        while (next_grid * sample_dt <= t_end) {
            push_row(res.traj, next_grid * sample_dt, x);
            ++next_grid;
        }
    } else if (!res.traj.t.empty() && res.traj.t.back() < t_end) {
        push_row(res.traj, t_end, x);
    }

    res.final_state = x;
    return res;
}

Trajectory run_rk4(const Network& net, double t_end, double h, double sample_dt) {
    const std::size_t S = net.n_species();
    const std::size_t R = net.reactions.size();

    std::vector<double> x(S), k1(S), k2(S), k3(S), k4(S), tmp(S);
    for (std::size_t s = 0; s < S; ++s) x[s] = static_cast<double>(net.x0[s]);

    auto deriv = [&](const std::vector<double>& state, std::vector<double>& out) {
        std::fill(out.begin(), out.end(), 0.0);
        for (std::size_t j = 0; j < R; ++j) {
            const double f = flux(net.reactions[j], state);
            if (f == 0.0) continue;
            const Reaction& r = net.reactions[j];
            for (std::size_t s = 0; s < S; ++s)
                if (r.net[s] != 0) out[s] += r.net[s] * f;
        }
    };

    // Integration is in doubles throughout; the recorded trajectory keeps them
    // as doubles too, and only the CSV writer decides on formatting.
    Trajectory tr;
    tr.n_species = S;
    std::vector<double> out_x;

    double t = 0.0;
    const double dt_out = sample_dt > 0.0 ? sample_dt : h;
    double next_out = dt_out;

    tr.t.push_back(0.0);
    out_x.insert(out_x.end(), x.begin(), x.end());

    while (t < t_end) {
        const double dt = std::min(h, t_end - t);
        if (dt <= 0.0) break;

        deriv(x, k1);
        for (std::size_t s = 0; s < S; ++s) tmp[s] = x[s] + 0.5 * dt * k1[s];
        deriv(tmp, k2);
        for (std::size_t s = 0; s < S; ++s) tmp[s] = x[s] + 0.5 * dt * k2[s];
        deriv(tmp, k3);
        for (std::size_t s = 0; s < S; ++s) tmp[s] = x[s] + dt * k3[s];
        deriv(tmp, k4);
        for (std::size_t s = 0; s < S; ++s)
            x[s] += dt / 6.0 * (k1[s] + 2.0 * k2[s] + 2.0 * k3[s] + k4[s]);
        t += dt;

        while (next_out <= t + 1e-12 && next_out <= t_end + 1e-12) {
            tr.t.push_back(next_out);
            out_x.insert(out_x.end(), x.begin(), x.end());
            next_out += dt_out;
        }
    }

    tr.x = std::move(out_x);
    return tr;
}

} // namespace crn

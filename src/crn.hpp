// Core chemical reaction network types.
//
// A network is stored as dense multiplicity vectors indexed by species, so one
// propensity expression covers every reaction with no special-casing.
#pragma once

#include <string>
#include <vector>

namespace crn {

using Count = long long;

struct Reaction {
    std::string label;          // human-readable, e.g. "X + Y -> 2Y"
    std::vector<int> reactants; // multiplicity per species
    std::vector<int> products;  // multiplicity per species
    std::vector<int> net;       // products - reactants, precomputed
    double k = 0.0;

    int order() const {
        int m = 0;
        for (int r : reactants) m += r;
        return m;
    }
};

struct Network {
    std::string name;
    std::vector<std::string> species;
    std::vector<Reaction> reactions;
    std::vector<Count> x0;
    double suggested_t_end = 30.0; // sensible default horizon, overridable with --t-end

    std::size_t n_species() const { return species.size(); }
    void finalize(); // fills in net[], validates sizes
};

// Stochastic propensity, a_j(x).
//
// Unimolecular:            k * n
// Bimolecular, distinct:   k * n_X * n_Y
// Bimolecular, same:       k * n_X * (n_X - 1) / 2
//
// The general form below is k * prod_s C(n_s, m_s) -- the number of distinct
// reactant combinations -- which reduces to all three.
inline double propensity(const Reaction& r, const std::vector<Count>& x) {
    double a = r.k;
    for (std::size_t s = 0; s < r.reactants.size(); ++s) {
        const int m = r.reactants[s];
        if (m == 0) continue;
        const Count n = x[s];
        if (n < m) return 0.0;
        for (int j = 0; j < m; ++j) a *= static_cast<double>(n - j);
        for (int j = 2; j <= m; ++j) a /= static_cast<double>(j);
    }
    return a;
}

} // namespace crn

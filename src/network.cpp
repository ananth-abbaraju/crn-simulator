#include "network.hpp"

#include <cctype>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>

namespace crn {
namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

// "2X" -> {2, "X"};  "X" -> {1, "X"}
std::pair<int, std::string> parse_term(const std::string& term) {
    std::size_t i = 0;
    while (i < term.size() && std::isdigit(static_cast<unsigned char>(term[i]))) ++i;
    const int mult = i > 0 ? std::stoi(term.substr(0, i)) : 1;
    const std::string name = trim(term.substr(i));
    if (name.empty()) throw std::runtime_error("bad term '" + term + "'");
    return {mult, name};
}

// "X + Y" -> multiplicity vector. "0" and "" denote the empty complex, which is
// how buffered sinks are written: `Y -> 0` means Y decays and the product is
// never tracked.
std::vector<int> parse_complex(const std::string& side,
                               std::map<std::string, std::size_t>& index,
                               std::vector<std::string>& species) {
    std::vector<int> mult;
    const std::string s = trim(side);
    if (s.empty() || s == "0" || s == "*") {
        mult.assign(species.size(), 0);
        return mult;
    }
    std::vector<std::pair<int, std::string>> terms;
    std::size_t pos = 0;
    while (pos <= s.size()) {
        const auto plus = s.find('+', pos);
        const std::string term = trim(s.substr(pos, plus - pos));
        if (!term.empty()) terms.push_back(parse_term(term));
        if (plus == std::string::npos) break;
        pos = plus + 1;
    }
    for (const auto& [m, name] : terms) {
        if (!index.count(name)) {
            index[name] = species.size();
            species.push_back(name);
        }
    }
    mult.assign(species.size(), 0);
    for (const auto& [m, name] : terms) mult[index[name]] += m;
    return mult;
}

Reaction make_reaction(const std::string& lhs, const std::string& rhs, double k,
                       std::map<std::string, std::size_t>& index,
                       std::vector<std::string>& species) {
    Reaction r;
    r.k = k;
    r.reactants = parse_complex(lhs, index, species);
    r.products = parse_complex(rhs, index, species);
    r.label = trim(lhs) + " -> " + trim(rhs);
    return r;
}

Network build(const std::string& name,
              const std::vector<std::tuple<std::string, std::string, double>>& rxns,
              const std::vector<std::pair<std::string, Count>>& init,
              double t_end) {
    Network n;
    n.name = name;
    std::map<std::string, std::size_t> index;
    for (const auto& [lhs, rhs, k] : rxns)
        n.reactions.push_back(make_reaction(lhs, rhs, k, index, n.species));
    n.x0.assign(n.species.size(), 0);
    for (const auto& [sp, c] : init) {
        auto it = index.find(sp);
        if (it == index.end())
            throw std::runtime_error("initial count for unknown species '" + sp + "'");
        n.x0[it->second] = c;
    }
    n.suggested_t_end = t_end;
    n.finalize();
    return n;
}

// Lotka-Volterra, A + X -> 2X, X + Y -> 2Y, Y -> B.
//
// A is a chemostat: held at constant concentration by an external reservoir, so
// it is folded into k1 and the reaction becomes pseudo-first-order in X. B is a
// pure sink and is never tracked. Tracking either as a real count kills the
// oscillation immediately -- A depletes and prey growth stops.
//
// Deterministic fixed point is at X = k3/k2 = 120, Y = k1/k2 = 200. The initial
// condition below sits off that point so the orbit is visible from t = 0.
Network lotka_volterra() {
    return build("lotka-volterra",
                 {{"X", "2X", 1.0},        // A + X -> 2X, with [A] folded into k
                  {"X + Y", "2Y", 0.005},
                  {"Y", "0", 0.6}},        // Y -> B, B untracked
                 {{"X", 120}, {"Y", 140}}, 30.0);
}

// Approximate majority: X + Y -> 2B, B + X -> 2X, B + Y -> 2Y.
//
// Total population is conserved. The two consensus states (all X, all Y) are
// the only absorbing states, and the network reaches the one held by the
// initial majority with high probability.
Network approximate_majority() {
    return build("approximate-majority",
                 {{"X + Y", "2B", 1.0},
                  {"B + X", "2X", 1.0},
                  {"B + Y", "2Y", 1.0}},
                 {{"X", 100}, {"Y", 90}, {"B", 0}}, 0.2);
}

} // namespace

void Network::finalize() {
    const std::size_t n = species.size();
    if (x0.size() != n) x0.resize(n, 0);
    for (auto& r : reactions) {
        r.reactants.resize(n, 0);
        r.products.resize(n, 0);
        r.net.assign(n, 0);
        for (std::size_t s = 0; s < n; ++s)
            r.net[s] = r.products[s] - r.reactants[s];
    }
}

std::vector<std::string> builtin_names() {
    return {"lotka-volterra", "approximate-majority"};
}

Network load_network(const std::string& spec) {
    if (spec == "lotka-volterra" || spec == "lv") return lotka_volterra();
    if (spec == "approximate-majority" || spec == "am") return approximate_majority();

    std::string msg = "unknown network '" + spec + "'. built-ins:";
    for (const auto& b : builtin_names()) msg += " " + b;
    throw std::runtime_error(msg);
}

} // namespace crn

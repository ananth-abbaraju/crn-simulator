#pragma once

#include "crn.hpp"

#include <string>
#include <vector>

namespace crn {

// Built-in network names, for --list and error messages.
std::vector<std::string> builtin_names();

// Resolves `spec` as a built-in name. Throws std::runtime_error on failure.
Network load_network(const std::string& spec);

} // namespace crn

#pragma once

#include "crn.hpp"

#include <string>
#include <vector>

namespace crn {

// Built-in network names, for --list and error messages.
std::vector<std::string> builtin_names();

// Resolves `spec` as a built-in name, or as a path to a .crn file.
// Throws std::runtime_error on failure.
Network load_network(const std::string& spec);

Network parse_crn(const std::string& text, const std::string& origin);

} // namespace crn

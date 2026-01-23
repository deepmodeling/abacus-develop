#include "input_help.h"
#include "input_help_data.h" // Generated header
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <vector>

namespace ModuleIO {

// Static member definitions
std::map<std::string, ParameterMetadata> ParameterHelp::registry_;
std::map<std::string, std::string> ParameterHelp::lowercase_to_actual_;
std::once_flag ParameterHelp::init_flag_;

void ParameterHelp::initialize() {
    std::call_once(init_flag_, build_registry);
}

void ParameterHelp::build_registry() {
    // Copy data from generated constexpr array to map for fast lookup
    for (size_t i = 0; i < PARAMETER_COUNT; ++i) {
        const auto& info = PARAMETER_DATA[i];
        ParameterMetadata meta;
        meta.name = info.name;
        meta.category = info.category;
        meta.type = info.type;
        meta.description = info.description;
        meta.default_value = info.default_value;
        meta.unit = info.unit ? info.unit : "";
        meta.availability = info.availability ? info.availability : "";

        // Pre-compute lowercase name for fast fuzzy matching
        meta.name_lowercase = to_lowercase(info.name);

        registry_[meta.name] = meta;

        // Pre-compute lowercase to actual name mapping for O(log n) case-insensitive lookup
        lowercase_to_actual_[meta.name_lowercase] = meta.name;
    }
}

bool ParameterHelp::show_parameter_help(const std::string& key) {
    initialize();

    // Use optimized case-insensitive lookup
    auto it = find_case_insensitive(key);

    if (it == registry_.end()) {
        return false;
    }

    const auto& meta = it->second;

    // Display formatted help information
    std::cout << "\n";
    std::cout << "Parameter: " << meta.name << "\n";
    std::cout << "Type:      " << meta.type << "\n";

    if (!meta.default_value.empty()) {
        std::cout << "Default:   " << meta.default_value << "\n";
    }

    if (!meta.category.empty()) {
        std::cout << "Category:  " << meta.category << "\n";
    }

    if (!meta.unit.empty()) {
        std::cout << "Unit:      " << meta.unit << "\n";
    }

    if (!meta.availability.empty()) {
        std::cout << "Availability: " << meta.availability << "\n";
    }

    std::cout << "\nDescription:\n";

    // Word-wrap description at 70 characters
    std::string desc = meta.description;
    size_t pos = 0;
    size_t line_start = 0;
    const size_t max_width = 70;

    while (pos < desc.length()) {
        // Find next space after max_width characters
        if (pos - line_start >= max_width) {
            size_t space_pos = desc.rfind(' ', pos);
            if (space_pos != std::string::npos && space_pos > line_start) {
                desc[space_pos] = '\n';
                line_start = space_pos + 1;
                pos = space_pos + 1;
            } else {
                pos++;
            }
        } else {
            pos++;
        }
    }

    // Indent each line by 2 spaces
    std::cout << "  ";
    for (char c : desc) {
        std::cout << c;
        if (c == '\n') {
            std::cout << "  ";
        }
    }
    std::cout << "\n\n";

    return true;
}

std::vector<std::string> ParameterHelp::search_parameters(const std::string& query) {
    initialize();

    std::vector<std::string> results;
    std::string query_lower = to_lowercase(query);

    // Search for parameters with case-insensitive substring match
    for (const auto& pair : registry_) {
        std::string name_lower = to_lowercase(pair.first);
        if (name_lower.find(query_lower) != std::string::npos) {
            results.push_back(pair.first);
        }
    }

    // Sort results alphabetically
    std::sort(results.begin(), results.end());

    return results;
}

void ParameterHelp::show_general_help() {
    std::cout << "\n";
    std::cout << "ABACUS - Atomic-orbital Based Ab-initio Computation at UStc\n";
    std::cout << "\n";
    std::cout << "Usage: abacus [options]\n";
    std::cout << "  -v, -V, --version      Display version information\n";
    std::cout << "  -i, -I, --info         Display detailed build information\n";
    std::cout << "  -h, --help [param]     Display help for parameter (or this message)\n";
    std::cout << "  -s, --search <query>   Search for parameters matching query\n";
    std::cout << "  --check-input          Check input file syntax and exit\n";
    std::cout << "\n";
    std::cout << "Common INPUT parameters:\n";
    std::cout << "  calculation    - Calculation type (scf, relax, md, nscf, etc.)\n";
    std::cout << "  basis_type     - Basis set type (pw, lcao)\n";
    std::cout << "  ecutwfc        - Energy cutoff for wavefunctions (Ry)\n";
    std::cout << "  ks_solver      - Kohn-Sham solver (cg, dav, genelpa, etc.)\n";
    std::cout << "  scf_thr        - SCF convergence threshold\n";
    std::cout << "  pseudo_dir     - Directory containing pseudopotential files\n";
    std::cout << "\n";
    std::cout << "For a complete list of parameters, see documentation at:\n";
    std::cout << "https://abacus.deepmodeling.com/\n";
    std::cout << "\n";
    std::cout << "To search for parameters: abacus -s <keyword>\n";
    std::cout << "To get help on a parameter: abacus -h <parameter_name>\n";
    std::cout << "\n";
}

ParameterMetadata ParameterHelp::get_metadata(const std::string& key) {
    initialize();

    // Use optimized case-insensitive lookup
    auto it = find_case_insensitive(key);

    if (it != registry_.end()) {
        return it->second;  // Return copy
    }

    // Return empty metadata to indicate not found
    return ParameterMetadata();
}

std::string ParameterHelp::to_lowercase(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

std::map<std::string, ParameterMetadata>::const_iterator
ParameterHelp::find_case_insensitive(const std::string& key) {
    // Try exact match first
    auto it = registry_.find(key);
    if (it != registry_.end()) {
        return it;
    }

    // Try case-insensitive match using pre-computed mapping (O(log n))
    std::string key_lower = to_lowercase(key);
    auto lower_it = lowercase_to_actual_.find(key_lower);
    if (lower_it != lowercase_to_actual_.end()) {
        return registry_.find(lower_it->second);
    }

    return registry_.end();
}

int ParameterHelp::levenshtein_distance(const std::string& s1, const std::string& s2) {
    const size_t len1 = s1.size();
    const size_t len2 = s2.size();

    // Space-optimized algorithm: only need two rows instead of full matrix
    // This reduces memory usage from O(m*n) to O(n)
    std::vector<int> prev(len2 + 1);
    std::vector<int> curr(len2 + 1);

    // Initialize first row
    for (size_t j = 0; j <= len2; ++j) {
        prev[j] = j;
    }

    // Calculate distances row by row
    for (size_t i = 1; i <= len1; ++i) {
        curr[0] = i;
        for (size_t j = 1; j <= len2; ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;

            curr[j] = std::min({
                prev[j] + 1,        // deletion
                curr[j-1] + 1,      // insertion
                prev[j-1] + cost    // substitution
            });
        }
        std::swap(prev, curr);
    }

    return prev[len2];
}

std::vector<std::string> ParameterHelp::find_similar_parameters(const std::string& query,
                                                                  int max_suggestions,
                                                                  int max_distance) {
    initialize();

    // Store pairs of (distance, parameter_name)
    std::vector<std::pair<int, std::string>> candidates;

    std::string query_lower = to_lowercase(query);

    // Calculate distance for each parameter using pre-computed lowercase names
    for (const auto& pair : registry_) {
        const auto& meta = pair.second;
        int distance = levenshtein_distance(query_lower, meta.name_lowercase);

        // Only consider parameters within max_distance
        if (distance <= max_distance && distance > 0) {
            candidates.push_back({distance, pair.first});
        }
    }

    // Sort by distance (closest first)
    std::sort(candidates.begin(), candidates.end(),
              [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
                  if (a.first != b.first) {
                      return a.first < b.first;
                  }
                  return a.second < b.second;
              });

    // Extract parameter names, limit to max_suggestions
    std::vector<std::string> results;
    for (size_t i = 0; i < candidates.size() && i < static_cast<size_t>(max_suggestions); ++i) {
        results.push_back(candidates[i].second);
    }

    return results;
}

} // namespace ModuleIO

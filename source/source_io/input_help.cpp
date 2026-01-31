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

bool ParameterHelp::show_parameter_help(const std::string& key, std::ostream& os) {
    initialize();

    // Use optimized case-insensitive lookup
    auto it = find_case_insensitive(key);

    if (it == registry_.end()) {
        return false;
    }

    const auto& meta = it->second;

    // Display formatted help information
    os << "\n";
    os << "Parameter: " << meta.name << "\n";
    os << "Type:      " << meta.type << "\n";

    if (!meta.default_value.empty()) {
        os << "Default:   " << meta.default_value << "\n";
    }

    if (!meta.category.empty()) {
        os << "Category:  " << meta.category << "\n";
    }

    if (!meta.unit.empty()) {
        os << "Unit:      " << meta.unit << "\n";
    }

    if (!meta.availability.empty()) {
        os << "Availability: " << meta.availability << "\n";
    }

    os << "\nDescription:\n";

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
    os << "  ";
    for (char c : desc) {
        os << c;
        if (c == '\n') {
            os << "  ";
        }
    }
    os << "\n\n";

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

void ParameterHelp::show_general_help(std::ostream& os) {
    os << "\n";
    os << "ABACUS - Atomic-orbital Based Ab-initio Computation at UStc\n";
    os << "\n";
    os << "Usage: abacus [options]\n";
    os << "  -v, -V, --version      Display version information\n";
    os << "  -i, -I, --info         Display detailed build information\n";
    os << "  -h, --help [param]     Display help for parameter (or this message)\n";
    os << "  -s, --search <query>   Search for parameters matching query\n";
    os << "  --check-input          Check input file syntax and exit\n";
    os << "\n";
    os << "Common INPUT parameters:\n";
    os << "  calculation    - Calculation type (scf, relax, md, nscf, etc.)\n";
    os << "  basis_type     - Basis set type (pw, lcao)\n";
    os << "  ecutwfc        - Energy cutoff for wavefunctions (Ry)\n";
    os << "  ks_solver      - Kohn-Sham solver (cg, dav, genelpa, etc.)\n";
    os << "  scf_thr        - SCF convergence threshold\n";
    os << "  pseudo_dir     - Directory containing pseudopotential files\n";
    os << "\n";
    os << "For a complete list of parameters, see documentation at:\n";
    os << "https://abacus.deepmodeling.com/\n";
    os << "\n";
    os << "To search for parameters: abacus -s <keyword>\n";
    os << "To get help on a parameter: abacus -h <parameter_name>\n";
    os << "\n";
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

    // If max_distance is 0, return nothing (exact matches are excluded by design)
    if (max_distance == 0) {
        return std::vector<std::string>();
    }

    // Store tuples of (effective_distance, parameter_name)
    // Effective distance prioritizes prefix/substring matches over pure edit distance
    std::vector<std::pair<int, std::string>> candidates;

    std::string query_lower = to_lowercase(query);

    // Calculate distance for each parameter using pre-computed lowercase names
    for (const auto& pair : registry_) {
        const auto& meta = pair.second;
        const std::string& name_lower = meta.name_lowercase;

        int effective_distance;

        // Priority 1: Exact prefix match (e.g., "relax" matches "relax_new")
        // Give these the lowest effective distance (0)
        if (name_lower.size() > query_lower.size() &&
            name_lower.compare(0, query_lower.size(), query_lower) == 0 &&
            name_lower[query_lower.size()] == '_') {
            effective_distance = 0;
        }
        // Priority 2: Substring match (e.g., "cut" matches "ecutwfc")
        // Give these a low effective distance (1)
        else if (name_lower.find(query_lower) != std::string::npos) {
            effective_distance = 1;
        }
        // Priority 3: Use Levenshtein distance for fuzzy matching
        else {
            int distance = levenshtein_distance(query_lower, name_lower);
            // Only consider parameters within max_distance
            if (distance > max_distance || distance == 0) {
                continue;  // Skip exact matches (distance 0) and too-distant matches
            }
            effective_distance = distance + 10;  // Add offset to prioritize after prefix/substring
        }

        candidates.push_back({effective_distance, pair.first});
    }

    // Sort by effective distance (closest first), then alphabetically
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

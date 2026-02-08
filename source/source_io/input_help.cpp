#include "input_help.h"
#include "module_parameter/read_input.h" // For accessing Input_Item documentation
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace ModuleIO {

namespace {
// Constants for display formatting
constexpr size_t MAX_WIDTH = 70;
constexpr size_t INDENT_SIZE = 2;
constexpr size_t LIST_INDENT = 4;
constexpr size_t NESTED_LIST_INDENT = 6;

/**
 * Word-wrap text at specified width with given indentation.
 * @param text Text to wrap
 * @param width Maximum line width
 * @param indent Number of spaces to indent each line
 * @param first_indent Number of spaces for the first line (if different)
 * @return Word-wrapped text
 */
std::string word_wrap(const std::string& text, size_t width, size_t indent, size_t first_indent) {
    if (text.empty()) {
        return "";
    }

    std::string result;
    std::string indent_str(indent, ' ');
    std::string first_indent_str(first_indent, ' ');
    size_t col = first_indent;

    std::istringstream iss(text);
    std::string word;
    bool first = true;

    while (iss >> word) {
        if (first) {
            result = first_indent_str + word;
            col = first_indent + word.length();
            first = false;
        } else if (col + 1 + word.length() > width) {
            result += "\n" + indent_str + word;
            col = indent + word.length();
        } else {
            result += " " + word;
            col += 1 + word.length();
        }
    }
    return result;
}

/**
 * Word-wrap with same indent for all lines.
 */
std::string word_wrap(const std::string& text, size_t width, size_t indent) {
    return word_wrap(text, width, indent, indent);
}

/**
 * Format a structured description that may contain markers for lists and paragraphs.
 *
 * Markers:
 *   \n\n     - paragraph break (blank line)
 *   \n*      - top-level list item
 *   \n  *    - nested list item
 *   [NOTE]   - note/blockquote
 *
 * @param desc The description string with embedded markers
 * @return Formatted string for terminal display
 */
std::string format_structured_description(const std::string& desc) {
    // Check if description contains any structure markers
    bool has_structure = (desc.find("\n\n") != std::string::npos ||
                         desc.find("\n*") != std::string::npos ||
                         desc.find("[NOTE]") != std::string::npos);

    if (!has_structure) {
        // Simple case: just word-wrap with basic indent
        return word_wrap(desc, MAX_WIDTH, INDENT_SIZE);
    }

    std::string result;
    size_t pos = 0;
    std::string current_text;
    bool at_line_start = true;  // Track if we're at the start of a logical line

    while (pos < desc.length()) {
        // Check for paragraph break (\n\n)
        if (pos + 1 < desc.length() && desc[pos] == '\n' && desc[pos + 1] == '\n') {
            // Flush current text
            if (!current_text.empty()) {
                if (!result.empty() && result.back() != '\n') {
                    result += "\n";
                }
                result += word_wrap(current_text, MAX_WIDTH, INDENT_SIZE);
                current_text.clear();
            }
            result += "\n\n";  // Two newlines: one to end current line, one for blank line
            pos += 2;
            at_line_start = true;
            continue;
        }

        // Check for nested list item (\n  * or at line start with leading spaces and *)
        if ((pos + 3 < desc.length() && desc[pos] == '\n' &&
            desc[pos + 1] == ' ' && desc[pos + 2] == ' ' && desc[pos + 3] == '*') ||
            (at_line_start && pos + 2 < desc.length() &&
             desc[pos] == ' ' && desc[pos + 1] == ' ' && desc[pos + 2] == '*')) {
            // Flush current text
            if (!current_text.empty()) {
                if (!result.empty() && result.back() != '\n') {
                    result += "\n";
                }
                result += word_wrap(current_text, MAX_WIDTH, INDENT_SIZE);
                current_text.clear();
            }
            // Skip the marker
            if (desc[pos] == '\n') {
                pos += 4;  // Skip "\n  *"
            } else {
                pos += 3;  // Skip "  *"
            }
            // Skip any whitespace after *
            while (pos < desc.length() && desc[pos] == ' ') {
                pos++;
            }
            // Collect list item text until next marker or end
            std::string item_text;
            while (pos < desc.length()) {
                if (desc[pos] == '\n') {
                    break;  // Stop at any newline marker
                }
                item_text += desc[pos++];
            }
            // Format nested list item with deeper indentation (indent=6 for continuation, first_indent=4)
            if (!result.empty() && result.back() != '\n') {
                result += "\n";
            }
            std::string prefix = "- ";
            result += word_wrap(prefix + item_text, MAX_WIDTH, NESTED_LIST_INDENT, LIST_INDENT);
            at_line_start = false;
            continue;
        }

        // Check for top-level list item (\n* or * at line start)
        if ((pos + 1 < desc.length() && desc[pos] == '\n' && desc[pos + 1] == '*') ||
            (at_line_start && desc[pos] == '*')) {
            // Flush current text
            if (!current_text.empty()) {
                if (!result.empty() && result.back() != '\n') {
                    result += "\n";
                }
                result += word_wrap(current_text, MAX_WIDTH, INDENT_SIZE);
                current_text.clear();
            }
            // Skip the marker
            if (desc[pos] == '\n') {
                pos += 2;  // Skip "\n*"
            } else {
                pos += 1;  // Skip "*"
            }
            // Skip any whitespace after *
            while (pos < desc.length() && desc[pos] == ' ') {
                pos++;
            }
            // Collect list item text until next marker or end
            std::string item_text;
            while (pos < desc.length()) {
                if (desc[pos] == '\n') {
                    break;  // Stop at any newline marker
                }
                item_text += desc[pos++];
            }
            // Format list item with "  - " prefix (indent=4 for continuation, first_indent=2)
            if (!result.empty() && result.back() != '\n') {
                result += "\n";
            }
            std::string prefix = "- ";
            result += word_wrap(prefix + item_text, MAX_WIDTH, LIST_INDENT, INDENT_SIZE);
            at_line_start = false;
            continue;
        }

        // Check for [NOTE] marker
        if (pos + 6 <= desc.length() && desc.substr(pos, 6) == "[NOTE]") {
            // Flush current text
            if (!current_text.empty()) {
                if (!result.empty() && result.back() != '\n') {
                    result += "\n";
                }
                result += word_wrap(current_text, MAX_WIDTH, INDENT_SIZE);
                current_text.clear();
            }
            pos += 6;  // Skip "[NOTE]"
            // Skip any whitespace after [NOTE]
            while (pos < desc.length() && desc[pos] == ' ') {
                pos++;
            }
            // Collect note text until next marker or end
            std::string note_text;
            while (pos < desc.length()) {
                if (desc[pos] == '\n') {
                    break;  // Stop at any newline marker
                }
                note_text += desc[pos++];
            }
            // Format note with "Note: " prefix
            if (!result.empty() && result.back() != '\n') {
                result += "\n";
            }
            result += word_wrap("Note: " + note_text, MAX_WIDTH, INDENT_SIZE + 6, INDENT_SIZE);
            at_line_start = false;
            continue;
        }

        // Regular character - accumulate
        at_line_start = false;
        current_text += desc[pos++];
    }

    // Flush any remaining text
    if (!current_text.empty()) {
        if (!result.empty() && result.back() != '\n') {
            result += "\n";
        }
        result += word_wrap(current_text, MAX_WIDTH, INDENT_SIZE);
    }

    return result;
}
}  // anonymous namespace

// ---- YAML serialization helpers (anonymous namespace) ----
namespace {

/**
 * Emit a string as a YAML literal block scalar (|).
 * Each line is prefixed with `indent` spaces. Empty lines are preserved.
 */
void emit_block_scalar(std::ostream& os, const std::string& text, int indent) {
    os << "|\n";
    std::string prefix(indent, ' ');
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.empty()) {
            os << "\n";
        } else {
            os << prefix << line << "\n";
        }
    }
    // If the text doesn't end with a newline, there's nothing extra to emit
}

/**
 * Return the value double-quoted with escaping if it contains
 * YAML-special characters, or unquoted if safe.
 */
std::string yaml_quote_if_needed(const std::string& value) {
    if (value.empty()) {
        return "\"\"";
    }

    // Check for YAML boolean keywords (case-insensitive)
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower == "true" || lower == "false" || lower == "yes" || lower == "no"
        || lower == "on" || lower == "off" || lower == "null" || lower == "~"
        || lower == ".inf" || lower == "-.inf" || lower == ".nan") {
        return "\"" + value + "\"";
    }

    // Quote numeric-looking values so YAML parsers keep them as strings.
    // Matches integers (0, -1, 0x1a, 0o17), floats (1.0, -3.14, 1.0e-6), etc.
    {
        bool all_numeric_chars = true;
        for (char c : value) {
            if (!std::isdigit(static_cast<unsigned char>(c))
                && c != '.' && c != '-' && c != '+' && c != 'e' && c != 'E'
                && c != 'x' && c != 'X' && c != 'o' && c != 'O'
                && c != 'a' && c != 'b' && c != 'c' && c != 'd' && c != 'f'
                && c != 'A' && c != 'B' && c != 'C' && c != 'D' && c != 'F') {
                all_numeric_chars = false;
                break;
            }
        }
        if (all_numeric_chars) {
            return "\"" + value + "\"";
        }
    }

    // Check for characters that need quoting
    bool needs_quoting = false;
    for (char c : value) {
        if (c == ':' || c == '#' || c == '[' || c == ']' ||
            c == '{' || c == '}' || c == '\\' || c == '*' ||
            c == '&' || c == '!' || c == '|' || c == '>' ||
            c == '\'' || c == '"' || c == '%' || c == '@' ||
            c == '`' || c == ',' || c == '\n' || c == '\r') {
            needs_quoting = true;
            break;
        }
    }

    // Also quote if starts/ends with whitespace or starts with special chars
    if (!needs_quoting) {
        if (value.front() == ' ' || value.back() == ' ' ||
            value.front() == '-' || value.front() == '?' ||
            value.front() == '{' || value.front() == '[') {
            needs_quoting = true;
        }
    }

    if (!needs_quoting) {
        return value;
    }

    // Double-quote with escaping
    std::string result = "\"";
    for (char c : value) {
        if (c == '"') {
            result += "\\\"";
        } else if (c == '\\') {
            result += "\\\\";
        } else if (c == '\n') {
            result += "\\n";
        } else if (c == '\r') {
            result += "\\r";
        } else if (c == '\t') {
            result += "\\t";
        } else {
            result += c;
        }
    }
    result += "\"";
    return result;
}

}  // anonymous namespace (YAML helpers)

// Static member definitions
std::map<std::string, ParameterMetadata> ParameterHelp::registry_;
std::map<std::string, std::string> ParameterHelp::lowercase_to_actual_;
std::once_flag ParameterHelp::init_flag_;

void ParameterHelp::initialize() {
    std::call_once(init_flag_, build_registry);
}

void ParameterHelp::build_registry() {
    // Create a ReadInput instance to access Input_Item documentation
    // Use rank -1 to indicate help-system mode (no MPI operations)
    ReadInput reader(-1);

    // Build registry from Input_Item objects
    const auto& input_lists = reader.get_input_lists();

    for (const auto& pair : input_lists) {
        const auto& item = pair.second;
        ParameterMetadata meta;
        meta.name = item.label;
        meta.category = item.category;
        meta.type = item.type;
        meta.description = item.description;
        meta.default_value = item.default_value;
        meta.unit = item.unit;
        meta.availability = item.availability;

        // Pre-compute lowercase name for fast fuzzy matching
        meta.name_lowercase = to_lowercase(item.label);

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

    // Use structured formatting for description
    os << format_structured_description(meta.description) << "\n\n";

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
    os << "  --generate-parameters-yaml\n";
    os << "                         Dump all parameter metadata as YAML\n";
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

void ParameterHelp::generate_yaml(std::ostream& os) {
    // Create a ReadInput instance to access Input_Item documentation
    // Use rank -1 to indicate help-system mode (no MPI operations)
    ReadInput reader(-1);

    const auto& input_lists = reader.get_input_lists();

    os << "# Auto-generated by: abacus --generate-parameters-yaml\n";
    os << "# Do not edit manually.\n";
    os << "parameters:\n";

    for (const auto& pair : input_lists) {
        const auto& item = pair.second;

        // Skip items without a category (undocumented internal items)
        if (item.category.empty()) {
            continue;
        }

        os << "  - name: " << yaml_quote_if_needed(item.label) << "\n";
        os << "    category: " << yaml_quote_if_needed(item.category) << "\n";
        os << "    type: " << yaml_quote_if_needed(item.type) << "\n";

        // Description uses block scalar
        os << "    description: ";
        if (item.description.empty()) {
            os << "\"\"\n";
        } else {
            emit_block_scalar(os, item.description, 6);
        }

        os << "    default_value: " << yaml_quote_if_needed(item.default_value) << "\n";
        os << "    unit: " << yaml_quote_if_needed(item.unit) << "\n";
        os << "    availability: " << yaml_quote_if_needed(item.availability) << "\n";
    }
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

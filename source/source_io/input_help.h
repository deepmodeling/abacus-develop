#ifndef INPUT_HELP_H
#define INPUT_HELP_H

#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace ModuleIO {

/**
 * @brief Metadata for a single INPUT parameter
 */
struct ParameterMetadata {
    std::string name;
    std::string type;
    std::string description;
    std::string default_value;
    std::string category;
    std::string unit;         // Empty string if no unit
    std::string availability; // Empty string if always available
    std::string name_lowercase; // Pre-computed lowercase for fast fuzzy matching
};

/**
 * @brief Help system for ABACUS INPUT parameters
 *
 * This class provides functionality to search for and display help information
 * about INPUT parameters. The parameter data is loaded from auto-generated
 * code that parses the documentation at build time.
 */
class ParameterHelp {
public:
    /**
     * @brief Initialize the help registry from generated data
     *
     * This function is called automatically on first use. It builds the
     * parameter registry from the generated PARAMETER_DATA array.
     */
    static void initialize();

    /**
     * @brief Display detailed help for a specific parameter
     *
     * @param key The parameter name to look up (case-sensitive)
     * @return true if parameter was found and help was displayed, false otherwise
     */
    static bool show_parameter_help(const std::string& key);

    /**
     * @brief Search for parameters matching a query string
     *
     * Performs case-insensitive substring matching on parameter names.
     *
     * @param query The search query string
     * @return Vector of matching parameter names (sorted alphabetically)
     */
    static std::vector<std::string> search_parameters(const std::string& query);

    /**
     * @brief Display general help message
     *
     * Shows usage information and lists commonly used parameters.
     */
    static void show_general_help();

    /**
     * @brief Get metadata for a specific parameter
     *
     * Returns a copy of the parameter metadata. Check if the returned
     * metadata has a non-empty name to verify the parameter was found.
     *
     * @param key The parameter name to look up (case-sensitive)
     * @return ParameterMetadata with empty name if not found, otherwise the parameter metadata
     *
     * Example:
     *   auto meta = ParameterHelp::get_metadata("ecutwfc");
     *   if (!meta.name.empty()) {
     *       // Parameter found, use meta.description, etc.
     *   }
     */
    static ParameterMetadata get_metadata(const std::string& key);

    /**
     * @brief Find similar parameter names for fuzzy matching
     *
     * Uses Levenshtein distance to find parameters with similar names.
     * Useful for suggesting corrections when a parameter is not found.
     *
     * @param query The parameter name to find similar matches for
     * @param max_suggestions Maximum number of suggestions to return (default: 5)
     * @param max_distance Maximum edit distance to consider (default: 3)
     * @return Vector of similar parameter names sorted by similarity
     */
    static std::vector<std::string> find_similar_parameters(const std::string& query,
                                                             int max_suggestions = 5,
                                                             int max_distance = 3);

private:
    static std::map<std::string, ParameterMetadata> registry_;
    static std::map<std::string, std::string> lowercase_to_actual_;
    static std::once_flag init_flag_;

    /**
     * @brief Build the registry from generated PARAMETER_DATA
     *
     * This is called once during initialization to populate the registry
     * from the static constexpr data array. Thread-safe via std::call_once.
     */
    static void build_registry();

    /**
     * @brief Find parameter with case-insensitive matching
     *
     * Uses pre-computed lowercase mappings for O(log n) performance.
     *
     * @param key The parameter name to look up (any case)
     * @return Iterator to the parameter in registry_, or registry_.end() if not found
     */
    static std::map<std::string, ParameterMetadata>::const_iterator
    find_case_insensitive(const std::string& key);

    /**
     * @brief Convert string to lowercase for case-insensitive comparison
     */
    static std::string to_lowercase(const std::string& str);

    /**
     * @brief Calculate Levenshtein distance between two strings
     *
     * Returns the minimum number of single-character edits (insertions,
     * deletions, or substitutions) required to change one string into another.
     *
     * @param s1 First string
     * @param s2 Second string
     * @return Edit distance between the strings
     */
    static int levenshtein_distance(const std::string& s1, const std::string& s2);
};

} // namespace ModuleIO

#endif // INPUT_HELP_H

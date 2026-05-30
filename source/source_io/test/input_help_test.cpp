#include "gtest/gtest.h"
#include "source_io/input_help.h"
#include <sstream>
#include <iostream>

// RAII helper class to safely redirect std::cout and restore it even if exceptions occur
class ScopedCoutRedirect {
public:
    explicit ScopedCoutRedirect(std::streambuf* new_buf)
        : old_buf_(std::cout.rdbuf(new_buf)) {}

    ~ScopedCoutRedirect() {
        std::cout.rdbuf(old_buf_);
    }

    // Prevent copying
    ScopedCoutRedirect(const ScopedCoutRedirect&) = delete;
    ScopedCoutRedirect& operator=(const ScopedCoutRedirect&) = delete;

private:
    std::streambuf* old_buf_;
};

// Test fixture for ParameterHelp tests
class ParameterHelpTest : public testing::Test {
protected:
    void SetUp() override {
        // Initialize help system before each test
        ModuleIO::ParameterHelp::initialize();
    }
};

// Test: Initialization works correctly
TEST_F(ParameterHelpTest, Initialization) {
    // If we get here, initialize() didn't crash
    SUCCEED();
}

// Test: Get metadata for known parameter
TEST_F(ParameterHelpTest, GetMetadataKnownParameter) {
    auto meta = ModuleIO::ParameterHelp::get_metadata("ecutwfc");
    ASSERT_FALSE(meta.name.empty());
    EXPECT_EQ(meta.name, "ecutwfc");
    EXPECT_EQ(meta.type, "Real");
    EXPECT_FALSE(meta.description.empty());
    EXPECT_FALSE(meta.unit.empty());
    EXPECT_EQ(meta.unit, "Ry");
}

// Test: Get metadata for unknown parameter
TEST_F(ParameterHelpTest, GetMetadataUnknownParameter) {
    auto meta = ModuleIO::ParameterHelp::get_metadata("nonexistent_parameter_xyz");
    EXPECT_TRUE(meta.name.empty());
}

// Test: Show parameter help for known parameter
TEST_F(ParameterHelpTest, ShowParameterHelpKnown) {
    std::ostringstream captured;
    bool result;

    {
        // RAII ensures cout is restored even if an exception is thrown
        ScopedCoutRedirect redirect(captured.rdbuf());
        result = ModuleIO::ParameterHelp::show_parameter_help("calculation");
    }

    EXPECT_TRUE(result);
    std::string output = captured.str();
    EXPECT_NE(output.find("Parameter: calculation"), std::string::npos);
    EXPECT_NE(output.find("Type:"), std::string::npos);
    EXPECT_NE(output.find("Description:"), std::string::npos);
}

// Test: Show parameter help for unknown parameter
TEST_F(ParameterHelpTest, ShowParameterHelpUnknown) {
    bool result = ModuleIO::ParameterHelp::show_parameter_help("this_param_does_not_exist");
    EXPECT_FALSE(result);
}

// Test: Search for parameters (exact match)
TEST_F(ParameterHelpTest, SearchParametersExact) {
    auto results = ModuleIO::ParameterHelp::search_parameters("ecutwfc");
    ASSERT_EQ(results.size(), 1);
    EXPECT_EQ(results[0], "ecutwfc");
}

// Test: Search for parameters (partial match)
TEST_F(ParameterHelpTest, SearchParametersPartial) {
    auto results = ModuleIO::ParameterHelp::search_parameters("ecut");
    EXPECT_GT(results.size(), 1); // Should find multiple matches like ecutwfc, ecutrho, etc.

    // Check that results are sorted
    for (size_t i = 1; i < results.size(); ++i) {
        EXPECT_LT(results[i-1], results[i]);
    }
}

// Test: Search for parameters (case insensitive)
TEST_F(ParameterHelpTest, SearchParametersCaseInsensitive) {
    auto results_lower = ModuleIO::ParameterHelp::search_parameters("ecutwfc");
    auto results_upper = ModuleIO::ParameterHelp::search_parameters("ECUTWFC");
    auto results_mixed = ModuleIO::ParameterHelp::search_parameters("EcUtWfC");

    EXPECT_EQ(results_lower.size(), results_upper.size());
    EXPECT_EQ(results_lower.size(), results_mixed.size());
    EXPECT_EQ(results_lower, results_upper);
    EXPECT_EQ(results_lower, results_mixed);
}

// Test: Search for parameters (no matches)
TEST_F(ParameterHelpTest, SearchParametersNoMatch) {
    auto results = ModuleIO::ParameterHelp::search_parameters("xyzabc123notfound");
    EXPECT_EQ(results.size(), 0);
}

// Test: Show general help
TEST_F(ParameterHelpTest, ShowGeneralHelp) {
    std::ostringstream captured;

    {
        // RAII ensures cout is restored even if an exception is thrown
        ScopedCoutRedirect redirect(captured.rdbuf());
        ModuleIO::ParameterHelp::show_general_help();
    }

    std::string output = captured.str();
    EXPECT_NE(output.find("ABACUS"), std::string::npos);
    EXPECT_NE(output.find("Usage:"), std::string::npos);
    EXPECT_NE(output.find("-h"), std::string::npos);
    EXPECT_NE(output.find("-s"), std::string::npos);
    EXPECT_NE(output.find("Common INPUT parameters:"), std::string::npos);
}

// Test: Verify multiple common parameters exist
TEST_F(ParameterHelpTest, CommonParametersExist) {
    std::vector<std::string> common_params = {
        "calculation", "basis_type", "ecutwfc", "ks_solver",
        "scf_thr", "pseudo_dir", "nspin", "nbands"
    };

    for (const auto& param : common_params) {
        auto meta = ModuleIO::ParameterHelp::get_metadata(param);
        EXPECT_FALSE(meta.name.empty()) << "Parameter " << param << " should exist";
        if (!meta.name.empty()) {
            EXPECT_EQ(meta.name, param);
        }
    }
}

// Test: Verify metadata fields are populated
TEST_F(ParameterHelpTest, MetadataFieldsPopulated) {
    auto meta = ModuleIO::ParameterHelp::get_metadata("symmetry_prec");
    ASSERT_FALSE(meta.name.empty());

    EXPECT_FALSE(meta.name.empty());
    EXPECT_FALSE(meta.type.empty());
    EXPECT_FALSE(meta.description.empty());
    EXPECT_FALSE(meta.default_value.empty());
    EXPECT_FALSE(meta.category.empty());

    // Unit should be "Bohr" for symmetry_prec
    EXPECT_EQ(meta.unit, "Bohr");
}

// Test: Case-insensitive parameter lookup
TEST_F(ParameterHelpTest, CaseInsensitiveLookup) {
    // Test with different capitalizations
    auto meta_lower = ModuleIO::ParameterHelp::get_metadata("ecutwfc");
    auto meta_upper = ModuleIO::ParameterHelp::get_metadata("ECUTWFC");
    auto meta_mixed = ModuleIO::ParameterHelp::get_metadata("EcUtWfC");

    ASSERT_FALSE(meta_lower.name.empty());
    ASSERT_FALSE(meta_upper.name.empty());
    ASSERT_FALSE(meta_mixed.name.empty());

    EXPECT_EQ(meta_lower.name, "ecutwfc");
    EXPECT_EQ(meta_upper.name, "ecutwfc");
    EXPECT_EQ(meta_mixed.name, "ecutwfc");
}

// Test: Fuzzy matching - single character typo
TEST_F(ParameterHelpTest, FuzzyMatchingSingleCharTypo) {
    // Missing 'c' at the end
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("ecutwf", 5, 3);
    EXPECT_GT(results.size(), 0);
    EXPECT_EQ(results[0], "ecutwfc"); // Should be the closest match
}

// Test: Fuzzy matching - extra character
TEST_F(ParameterHelpTest, FuzzyMatchingExtraChar) {
    // Extra 's' at the end
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("exx_hybrid_steps", 5, 3);
    ASSERT_GT(results.size(), 0);
    EXPECT_EQ(results[0], "exx_hybrid_step");
}

// Test: Fuzzy matching - swapped characters
TEST_F(ParameterHelpTest, FuzzyMatchingSwappedChars) {
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("scf_htr", 5, 3);
    EXPECT_GT(results.size(), 0);
    // scf_thr should be in the results
    bool found = false;
    for (const auto& r : results) {
        if (r == "scf_thr") {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

// Test: Fuzzy matching - case insensitive
TEST_F(ParameterHelpTest, FuzzyMatchingCaseInsensitive) {
    auto results_lower = ModuleIO::ParameterHelp::find_similar_parameters("ecutwf", 5, 3);
    auto results_upper = ModuleIO::ParameterHelp::find_similar_parameters("ECUTWF", 5, 3);
    auto results_mixed = ModuleIO::ParameterHelp::find_similar_parameters("EcUtWf", 5, 3);

    EXPECT_EQ(results_lower.size(), results_upper.size());
    EXPECT_EQ(results_lower.size(), results_mixed.size());
    EXPECT_EQ(results_lower, results_upper);
    EXPECT_EQ(results_lower, results_mixed);
}

// Test: Fuzzy matching - multiple suggestions
TEST_F(ParameterHelpTest, FuzzyMatchingMultipleSuggestions) {
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("relax_met", 5, 3);
    EXPECT_GT(results.size(), 1); // Should find multiple matches
    // Results should be sorted by distance (closest first)
    // "relax_met" to "relax_new": distance 2 (m->n, t->w)
    // "relax_met" to "relax_method": distance 3 (insert h, o, d)
    // Note: Actual results depend on which parameters exist in the parameter database
}

// Test: Fuzzy matching - no close matches
TEST_F(ParameterHelpTest, FuzzyMatchingNoCloseMatches) {
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("completely_wrong_parameter_xyz", 5, 3);
    EXPECT_EQ(results.size(), 0); // Should find nothing within distance 3
}

// Test: Fuzzy matching - max suggestions limit
TEST_F(ParameterHelpTest, FuzzyMatchingMaxSuggestions) {
    // Use a parameter that has many similar matches
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("md", 3, 2);
    EXPECT_LE(results.size(), 3); // Should not exceed max_suggestions
}

// Test: Fuzzy matching - max distance limit
TEST_F(ParameterHelpTest, FuzzyMatchingMaxDistance) {
    // With max_distance=1, should only find very close matches
    auto results = ModuleIO::ParameterHelp::find_similar_parameters("ecutwf", 5, 1);
    EXPECT_GT(results.size(), 0);

    // With max_distance=0, should find nothing (exact match excluded)
    auto results_zero = ModuleIO::ParameterHelp::find_similar_parameters("ecutwfc", 5, 0);
    EXPECT_EQ(results_zero.size(), 0);
}

// Test: Parameter with list items displays with bullets
TEST_F(ParameterHelpTest, ListFormattingDisplayed) {
    std::ostringstream captured;
    bool result;

    {
        ScopedCoutRedirect redirect(captured.rdbuf());
        result = ModuleIO::ParameterHelp::show_parameter_help("calculation");
    }

    EXPECT_TRUE(result);
    std::string output = captured.str();

    // Check that the output contains list item markers
    // The calculation parameter has items like "scf:", "nscf:", etc.
    EXPECT_NE(output.find("- scf:"), std::string::npos)
        << "Expected list item '- scf:' in output:\n" << output;
    EXPECT_NE(output.find("- relax:"), std::string::npos)
        << "Expected list item '- relax:' in output:\n" << output;
}

// Test: Lines do not exceed maximum width
TEST_F(ParameterHelpTest, WordWrapMaxWidth) {
    std::ostringstream captured;

    {
        ScopedCoutRedirect redirect(captured.rdbuf());
        ModuleIO::ParameterHelp::show_parameter_help("calculation");
    }

    std::string output = captured.str();
    std::istringstream iss(output);
    std::string line;

    // Check each line is within reasonable width (allowing some margin for list markers)
    const size_t max_line_width = 80;  // Allow some margin beyond 70
    while (std::getline(iss, line)) {
        EXPECT_LE(line.length(), max_line_width)
            << "Line exceeds max width: \"" << line << "\"";
    }
}

// Test: esolver_type has list items properly formatted
TEST_F(ParameterHelpTest, EsolverTypeListFormatting) {
    std::ostringstream captured;
    bool result;

    {
        ScopedCoutRedirect redirect(captured.rdbuf());
        result = ModuleIO::ParameterHelp::show_parameter_help("esolver_type");
    }

    EXPECT_TRUE(result);
    std::string output = captured.str();

    // esolver_type has items like "ksdft:", "ofdft:", etc.
    EXPECT_NE(output.find("ksdft:"), std::string::npos)
        << "Expected 'ksdft:' in output:\n" << output;
    EXPECT_NE(output.find("tddft:"), std::string::npos)
        << "Expected 'tddft:' in output:\n" << output;
}

// Test: symmetry parameter displays with list items
TEST_F(ParameterHelpTest, SymmetryListFormatting) {
    std::ostringstream captured;
    bool result;

    {
        ScopedCoutRedirect redirect(captured.rdbuf());
        result = ModuleIO::ParameterHelp::show_parameter_help("symmetry");
    }

    EXPECT_TRUE(result);
    std::string output = captured.str();

    // symmetry has items like "-1:", "0:", "1:"
    // These should be preserved in some form
    EXPECT_NE(output.find("-1:"), std::string::npos)
        << "Expected '-1:' in output:\n" << output;
}

// Test: Description content is preserved (not truncated)
TEST_F(ParameterHelpTest, DescriptionContentPreserved) {
    auto meta = ModuleIO::ParameterHelp::get_metadata("calculation");
    ASSERT_FALSE(meta.name.empty());

    // The calculation parameter should mention scf, relax, md, etc.
    EXPECT_NE(meta.description.find("scf"), std::string::npos)
        << "Expected 'scf' in description: " << meta.description;
    EXPECT_NE(meta.description.find("relax"), std::string::npos)
        << "Expected 'relax' in description: " << meta.description;
    EXPECT_NE(meta.description.find("md"), std::string::npos)
        << "Expected 'md' in description: " << meta.description;
}

// Test: generate_yaml() produces valid YAML header
TEST_F(ParameterHelpTest, GenerateYamlHeader) {
    std::ostringstream captured;
    ModuleIO::ParameterHelp::generate_yaml(captured);
    std::string output = captured.str();

    // Must contain the top-level "parameters:" key
    EXPECT_NE(output.find("parameters:"), std::string::npos)
        << "YAML output must contain 'parameters:' header";
}

// Test: generate_yaml() contains known parameters
TEST_F(ParameterHelpTest, GenerateYamlContainsKnownParams) {
    std::ostringstream captured;
    ModuleIO::ParameterHelp::generate_yaml(captured);
    std::string output = captured.str();

    EXPECT_NE(output.find("name: ecutwfc"), std::string::npos)
        << "YAML output must contain parameter 'ecutwfc'";
    EXPECT_NE(output.find("name: calculation"), std::string::npos)
        << "YAML output must contain parameter 'calculation'";
    EXPECT_NE(output.find("name: basis_type"), std::string::npos)
        << "YAML output must contain parameter 'basis_type'";
}

// Test: generate_yaml() has sufficient parameter count
TEST_F(ParameterHelpTest, GenerateYamlParameterCount) {
    std::ostringstream captured;
    ModuleIO::ParameterHelp::generate_yaml(captured);
    std::string output = captured.str();

    // Count occurrences of "  - name: " which marks each parameter entry
    size_t count = 0;
    size_t pos = 0;
    std::string marker = "  - name: ";
    while ((pos = output.find(marker, pos)) != std::string::npos) {
        ++count;
        pos += marker.size();
    }

    EXPECT_GT(count, 400u)
        << "YAML output should contain > 400 parameters, found " << count;
}

// Test: generate_yaml() uses block scalar syntax for descriptions
TEST_F(ParameterHelpTest, GenerateYamlBlockScalar) {
    std::ostringstream captured;
    ModuleIO::ParameterHelp::generate_yaml(captured);
    std::string output = captured.str();

    // Block scalar markers "description: |" should be present
    EXPECT_NE(output.find("description: |"), std::string::npos)
        << "YAML output must use block scalar syntax (|) for descriptions";
}

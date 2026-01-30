#include "gtest/gtest.h"
#include "../input_help.h"
#include <sstream>
#include <iostream>

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
    // Redirect stdout to capture output
    std::streambuf* old_cout = std::cout.rdbuf();
    std::ostringstream captured;
    std::cout.rdbuf(captured.rdbuf());

    bool result = ModuleIO::ParameterHelp::show_parameter_help("calculation");

    // Restore stdout
    std::cout.rdbuf(old_cout);

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
    // Redirect stdout to capture output
    std::streambuf* old_cout = std::cout.rdbuf();
    std::ostringstream captured;
    std::cout.rdbuf(captured.rdbuf());

    ModuleIO::ParameterHelp::show_general_help();

    // Restore stdout
    std::cout.rdbuf(old_cout);

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
    // Results should be sorted by distance
    // relax_method has distance 3 (replace 'hod' with nothing = 3 deletions)
    // relax_new has distance 3
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

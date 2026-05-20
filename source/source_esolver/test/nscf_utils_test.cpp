#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <sstream>

/************************************************
 *  unit test of NSCF (non-self-consistent field)
 *  related utility functions and logic
 ***********************************************/

/**
 * Tested functionality:
 * - NSCF mode detection (calculation == "nscf")
 * - skip_charge logic
 * - init_chg file vs atomic handling
 * - out_band parameter parsing
 * - NSCF-specific validation rules
 */

// ===================================================================
// 1. NSCF Mode Detection Tests
// ===================================================================

class NscfModeDetectionTest : public ::testing::Test
{
  protected:
    bool is_nscf(const std::string& calculation)
    {
        return calculation == "nscf";
    }

    bool should_skip_charge(const std::string& calculation)
    {
        return calculation == "nscf";
    }
};

TEST_F(NscfModeDetectionTest, ScfMode_IsNotNscf)
{
    EXPECT_FALSE(is_nscf("scf"));
    EXPECT_FALSE(should_skip_charge("scf"));
}

TEST_F(NscfModeDetectionTest, NscfMode_IsNscf)
{
    EXPECT_TRUE(is_nscf("nscf"));
    EXPECT_TRUE(should_skip_charge("nscf"));
}

TEST_F(NscfModeDetectionTest, OtherModes_AreNotNscf)
{
    EXPECT_FALSE(is_nscf("relax"));
    EXPECT_FALSE(is_nscf("cell-relax"));
    EXPECT_FALSE(is_nscf("md"));
    EXPECT_FALSE(is_nscf(""));
}

// ===================================================================
// 2. Charge Initialization Tests
// ===================================================================

class ChargeInitTest : public ::testing::Test
{
  protected:
    bool should_read_charge_from_file(const std::string& init_chg,
                                       const std::string& calculation)
    {
        // NSCF always reads charge from file
        if (calculation == "nscf") return true;
        // SCF with init_chg="file" reads from file
        if (init_chg == "file") return true;
        return false;
    }

    bool requires_charge_file(const std::string& calculation)
    {
        return calculation == "nscf";
    }
};

TEST_F(ChargeInitTest, Nscf_AlwaysReadsFromFile)
{
    EXPECT_TRUE(should_read_charge_from_file("atomic", "nscf"));
    EXPECT_TRUE(should_read_charge_from_file("file", "nscf"));
    EXPECT_TRUE(requires_charge_file("nscf"));
}

TEST_F(ChargeInitTest, Scf_WithAtomicInit)
{
    EXPECT_FALSE(should_read_charge_from_file("atomic", "scf"));
    EXPECT_FALSE(requires_charge_file("scf"));
}

TEST_F(ChargeInitTest, Scf_WithFileInit)
{
    EXPECT_TRUE(should_read_charge_from_file("file", "scf"));
    EXPECT_FALSE(requires_charge_file("scf"));
}

// ===================================================================
// 3. Charge Density Filename Tests
// ===================================================================

class ChargeFilenameTest : public ::testing::Test
{
  protected:
    std::string make_charge_filename(const std::string& suffix)
    {
        return suffix + "-CHARGE-DENSITY.restart";
    }

    std::string extract_suffix_from_input(const std::string& input_content)
    {
        // Simplified: look for "suffix" keyword
        size_t pos = input_content.find("suffix");
        if (pos == std::string::npos) return "OUT";
        size_t start = input_content.find_first_not_of(" \t", pos + 6);
        size_t end = input_content.find_first_of(" \t\n", start);
        return input_content.substr(start, end - start);
    }
};

TEST_F(ChargeFilenameTest, DefaultSuffix)
{
    EXPECT_EQ(make_charge_filename("OUT"), "OUT-CHARGE-DENSITY.restart");
}

TEST_F(ChargeFilenameTest, CustomSuffix)
{
    EXPECT_EQ(make_charge_filename("autotest"), "autotest-CHARGE-DENSITY.restart");
    EXPECT_EQ(make_charge_filename("my_calc"), "my_calc-CHARGE-DENSITY.restart");
}

TEST_F(ChargeFilenameTest, ExtractSuffixFromInput)
{
    std::string input = "INPUT_PARAMETERS\nsuffix    autotest\necutwfc    20\n";
    EXPECT_EQ(extract_suffix_from_input(input), "autotest");
}

TEST_F(ChargeFilenameTest, ExtractSuffix_DefaultWhenMissing)
{
    std::string input = "INPUT_PARAMETERS\necutwfc    20\n";
    EXPECT_EQ(extract_suffix_from_input(input), "OUT");
}

// ===================================================================
// 4. Band Output Parameter Tests
// ===================================================================

class BandOutputParamTest : public ::testing::Test
{
  protected:
    struct BandParams
    {
        bool enabled = false;
        int precision = 8;
    };

    BandParams parse_out_band(const std::vector<std::string>& values)
    {
        BandParams params;
        if (values.empty()) return params;

        // Parse boolean
        std::string val = values[0];
        params.enabled = (val == "1" || val == "true" || val == "on");

        if (values.size() > 1)
        {
            params.precision = std::stoi(values[1]);
        }
        return params;
    }
};

TEST_F(BandOutputParamTest, DisabledByDefault)
{
    auto params = parse_out_band({});
    EXPECT_FALSE(params.enabled);
    EXPECT_EQ(params.precision, 8);
}

TEST_F(BandOutputParamTest, EnabledWithDefaultPrecision)
{
    auto params = parse_out_band({"1"});
    EXPECT_TRUE(params.enabled);
    EXPECT_EQ(params.precision, 8);
}

TEST_F(BandOutputParamTest, EnabledWithCustomPrecision)
{
    auto params = parse_out_band({"1", "12"});
    EXPECT_TRUE(params.enabled);
    EXPECT_EQ(params.precision, 12);
}

TEST_F(BandOutputParamTest, VariousEnableValues)
{
    EXPECT_TRUE(parse_out_band({"true"}).enabled);
    EXPECT_TRUE(parse_out_band({"on"}).enabled);
    EXPECT_FALSE(parse_out_band({"0"}).enabled);
    EXPECT_FALSE(parse_out_band({"false"}).enabled);
}

// ===================================================================
// 5. NSCF Input Validation Tests
// ===================================================================

class NscfValidationTest : public ::testing::Test
{
  protected:
    struct ValidationIssue
    {
        std::string field;
        std::string message;
    };

    std::vector<ValidationIssue> validate_nscf_input(
        const std::string& calculation,
        const std::string& init_chg,
        const std::string& read_file_dir,
        bool has_charge_file)
    {
        std::vector<ValidationIssue> issues;

        if (calculation == "nscf")
        {
            if (init_chg != "file")
            {
                issues.push_back({"init_chg",
                    "NSCF calculation requires init_chg=file"});
            }
            if (!has_charge_file && read_file_dir.empty())
            {
                issues.push_back({"read_file_dir",
                    "NSCF calculation requires charge density file"});
            }
        }
        return issues;
    }
};

TEST_F(NscfValidationTest, ValidNscfInput)
{
    auto issues = validate_nscf_input("nscf", "file", "./", true);
    EXPECT_TRUE(issues.empty());
}

TEST_F(NscfValidationTest, Nscf_WrongInitChg)
{
    auto issues = validate_nscf_input("nscf", "atomic", "./", true);
    EXPECT_EQ(issues.size(), 1);
    EXPECT_EQ(issues[0].field, "init_chg");
}

TEST_F(NscfValidationTest, ScfInput_NoIssues)
{
    auto issues = validate_nscf_input("scf", "atomic", "", false);
    EXPECT_TRUE(issues.empty());
}

// ===================================================================
// 6. K-Point Grid Tests for NSCF
// ===================================================================

class NscfKPointTest : public ::testing::Test
{
  protected:
    struct KPointGrid
    {
        int nx, ny, nz;
        int offset_x, offset_y, offset_z;
    };

    // Parse Monkhorst-Pack grid from KPT line
    KPointGrid parse_mp_grid(const std::string& grid_line)
    {
        KPointGrid grid = {0, 0, 0, 0, 0, 0};
        std::istringstream iss(grid_line);
        iss >> grid.nx >> grid.ny >> grid.nz
            >> grid.offset_x >> grid.offset_y >> grid.offset_z;
        return grid;
    }

    int total_kpoints(const KPointGrid& grid)
    {
        return grid.nx * grid.ny * grid.nz;
    }
};

TEST_F(NscfKPointTest, ParseMonkhorstPackGrid)
{
    auto grid = parse_mp_grid("2 2 2 0 0 0");
    EXPECT_EQ(grid.nx, 2);
    EXPECT_EQ(grid.ny, 2);
    EXPECT_EQ(grid.nz, 2);
    EXPECT_EQ(total_kpoints(grid), 8);
}

TEST_F(NscfKPointTest, ParseGammaCenteredGrid)
{
    auto grid = parse_mp_grid("4 4 4 0 0 0");
    EXPECT_EQ(total_kpoints(grid), 64);
}

TEST_F(NscfKPointTest, ParseOffsetGrid)
{
    auto grid = parse_mp_grid("3 3 3 1 1 1");
    EXPECT_EQ(grid.offset_x, 1);
    EXPECT_EQ(grid.offset_y, 1);
    EXPECT_EQ(grid.offset_z, 1);
}

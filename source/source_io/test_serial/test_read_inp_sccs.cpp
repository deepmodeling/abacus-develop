#include "source_io/module_parameter/read_inp_sccs.h"
#include "source_io/module_parameter/input_parameter.h"
#include <gtest/gtest.h>

TEST(ReadInpSccs, LegacyBooleansAndModelSelection)
{
    std::string error;
    const std::string yes[] = {"true", "TRUE", "yes", "On", ".true.", "1", "t", "y"};
    const std::string no[] = {"false", "FALSE", "no", "Off", ".false.", "0", "f", "n"};
    int model = 0;
    for (const std::string& value : yes)
    {
        ASSERT_TRUE(ModuleIO::parse_solvation_model(value, model, error));
        EXPECT_EQ(model, 1);
    }
    for (const std::string& value : no)
    {
        ASSERT_TRUE(ModuleIO::parse_solvation_model(value, model, error));
        EXPECT_EQ(model, 0);
    }
    ASSERT_TRUE(ModuleIO::parse_solvation_model("2", model, error));
    EXPECT_EQ(model, 2);
    EXPECT_FALSE(ModuleIO::parse_solvation_model("3", model, error));
    EXPECT_EQ(model, 2);
    EXPECT_FALSE(ModuleIO::parse_solvation_model("2.0", model, error));
}

TEST(ReadInpSccs, SupportedScopeAndNumericalValidation)
{
    Input_para input;
    input.imp_sol = 2;
    input.device = "cpu";
    input.basis_type = "pw";
    std::string error;
    ASSERT_TRUE(ModuleIO::validate_sccs_input(input, error)) << error;
    input.cal_force = true;
    EXPECT_FALSE(ModuleIO::validate_sccs_input(input, error));
    input.cal_force = false;
    input.assume_isolated = "pcc_0d";
    EXPECT_FALSE(ModuleIO::validate_sccs_input(input, error));
    input.assume_isolated = "none";
    input.sccs_rho_min = input.sccs_rho_max;
    EXPECT_FALSE(ModuleIO::validate_sccs_input(input, error));
    input.sccs_rho_min = 1e-4;
    input.sccs_preset = "unknown";
    EXPECT_FALSE(ModuleIO::validate_sccs_input(input, error));
    input.sccs_preset = "vacuum";
    input.basis_type = "lcao";
    EXPECT_TRUE(ModuleIO::validate_sccs_input(input, error));
    input.imp_sol = 1;
    input.sccs_preset = "unused";
    EXPECT_TRUE(ModuleIO::validate_sccs_input(input, error));
}

#include "source_io/module_parameter/read_inp_sccs.h"
#include "source_io/module_parameter/input_parameter.h"
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

TEST(ReadInpSccs, LegacyBooleansAndModelSelection)
{
    std::string error;
    const std::string yes[] = {"true", "TRUE", "yes", "On", ".true.", "1", "t", "y"};
    const std::string no[] = {"false", "FALSE", "no", "Off", ".false.", "0", "f", "n"};
    int model = 0;
    for (const std::string& value : yes)
    {
        ASSERT_TRUE(ModuleIO::parse_imp_sol(value, model, error));
        EXPECT_EQ(model, 1);
    }
    for (const std::string& value : no)
    {
        ASSERT_TRUE(ModuleIO::parse_imp_sol(value, model, error));
        EXPECT_EQ(model, 0);
    }
    ASSERT_TRUE(ModuleIO::parse_imp_sol("2", model, error));
    EXPECT_EQ(model, 2);
    EXPECT_FALSE(ModuleIO::parse_imp_sol("3", model, error));
    EXPECT_EQ(model, 2);
    EXPECT_FALSE(ModuleIO::parse_imp_sol("2.0", model, error));
}

TEST(ReadInpSccs, SupportedScopeAndNumericalValidation)
{
    Input_para input;
    input.imp_sol = 2;
    input.device = "cpu";
    input.basis_type = "pw";
    std::string error;
    ASSERT_TRUE(ModuleIO::check_sccs_inp(input, error)) << error;
    input.cal_force = true;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.cal_stress = true;
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.cal_stress = false;
    input.cal_force = false;
    input.assume_isolated = "pcc_0d";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.cal_force = true;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.cal_force = false;
    input.assume_isolated = "pcc_2d";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.cal_force = true;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.cal_force = false;
    input.calculation = "relax";
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.calculation = "scf";
    input.assume_isolated = "none";
    input.sccs_rho_min = input.sccs_rho_max;
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_rho_min = 1e-4;
    input.sccs_preset = "unknown";
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_preset = "vacuum";
    input.basis_type = "lcao";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.imp_sol = 1;
    input.sccs_preset = "unused";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
}

TEST(ReadInpSccs, LowpassRequiresPairedValuesAndPcc)
{
    Input_para input;
    input.imp_sol = 2;
    std::string error;
    EXPECT_DOUBLE_EQ(input.sccs_lowpass_p1, -1.0);
    EXPECT_DOUBLE_EQ(input.sccs_lowpass_p2, -1.0);
    input.sccs_lowpass_p1 = 10.0;
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_lowpass_p2 = 5.0;
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.assume_isolated = "pcc_0d";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.assume_isolated = "pcc_2d";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_lowpass_p1 = 0.0;
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_lowpass_p2 = 0.0;
    input.assume_isolated = "none";
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
}

TEST(ReadInpSccs, FullCavityWidthsAndModeValidation)
{
    Input_para input;
    input.imp_sol = 2;
    input.sccs_solvent_mode = "full";
    std::string error;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_corespread = {0.7, 0.0, -1.0};
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_corespread = {0.0};
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_corespread.clear();
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
    input.sccs_corespread = {0.5};
    input.sccs_solvent_mode = "invalid";
    EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
}

TEST(ReadInpSccs, CoreSpreadWordsMustBeFiniteNumbers)
{
    std::vector<double> spreads = {0.5};
    std::string error;
    ASSERT_TRUE(ModuleIO::parse_spreads({"0.7", "-1", "1e-1"}, spreads, error)) << error;
    ASSERT_EQ(spreads.size(), 3u);
    EXPECT_DOUBLE_EQ(spreads[0], 0.7);
    EXPECT_DOUBLE_EQ(spreads[1], -1.0);
    EXPECT_DOUBLE_EQ(spreads[2], 0.1);
    const std::vector<std::string> invalid[] = {{"abc"}, {"0.5abc"}, {"0.5", "x"}, {"nan"}, {"inf"}, {""}};
    for (const std::vector<std::string>& words : invalid)
    {
        EXPECT_FALSE(ModuleIO::parse_spreads(words, spreads, error));
        EXPECT_NE(error.find("sccs_corespread"), std::string::npos);
        EXPECT_EQ(spreads.size(), 3u);
    }
}

// nan and inf pass the ordered range checks, so check_sccs_inp tests every real
// SCCS key for finiteness. Signed gamma and pressure stay valid, and imp_sol 0
// and 1 skip the SCCS checks.
TEST(ReadInpSccs, RejectsNonfiniteReals)
{
    Input_para input;
    input.imp_sol = 2;
    input.device = "cpu";
    input.basis_type = "pw";
    input.sccs_gamma = -1.0;
    input.sccs_pressure = -0.36;
    std::string error;
    ASSERT_TRUE(ModuleIO::check_sccs_inp(input, error)) << error;

    double* const reals[] = {&input.sccs_epsilon,     &input.sccs_rho_min,    &input.sccs_rho_max,
                             &input.sccs_gamma,       &input.sccs_pressure,   &input.sccs_tol_rms,
                             &input.sccs_tol_max,     &input.sccs_surface_eta, &input.sccs_lowpass_p1,
                             &input.sccs_lowpass_p2};
    const double infinity = std::numeric_limits<double>::infinity();
    const double nonfinite[] = {std::nan(""), infinity, -infinity};
    for (double* value : reals)
    {
        const double saved = *value;
        for (double bad : nonfinite)
        {
            *value = bad;
            EXPECT_FALSE(ModuleIO::check_sccs_inp(input, error));
            const std::size_t position = error.find("must be finite");
            EXPECT_NE(position, std::string::npos) << error;
        }
        *value = saved;
    }

    input.sccs_surface_eta = std::nan("");
    input.imp_sol = 0;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error)) << error;
    input.imp_sol = 1;
    EXPECT_TRUE(ModuleIO::check_sccs_inp(input, error)) << error;
}

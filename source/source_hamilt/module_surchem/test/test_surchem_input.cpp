#include "../surchem_input.h"
#include "../surchem.h"
#include "source_io/module_parameter/input_parameter.h"

#include <gtest/gtest.h>

TEST(SurchemInput, SelectsVacuumPccIndependentlyOfSolvent)
{
    Input_para input;
    UnitCell cell;
    input.imp_sol = 0;
    input.assume_isolated = "pcc_0d";
    const SurchemParameters parameters = ModuleSurchem::make_parameters(input, cell, 0.0, false, 4);
    EXPECT_FALSE(parameters.use_sccs);
    EXPECT_FALSE(parameters.use_legacy_solvent);
    EXPECT_EQ(parameters.pcc_boundary, ModuleSccs::Boundary::Pcc0d);
    EXPECT_EQ(parameters.pool_process_count, 4);
    EXPECT_TRUE(input.uses_surchem_correction());
}

TEST(SurchemInput, TransfersPresetAndSolverControls)
{
    Input_para input;
    UnitCell cell;
    input.imp_sol = 2;
    input.assume_isolated = "pcc_2d";
    input.sccs_preset = "vacuum";
    input.sccs_debug = 2;
    input.sccs_mixing_type = "pulay";
    input.sccs_mixing = 0.2;
    input.sccs_start_drho = 0.01;
    input.sccs_start_nmax = 12;
    const SurchemParameters parameters = ModuleSurchem::make_parameters(input, cell, 0.0, false, 2);
    EXPECT_TRUE(parameters.use_sccs);
    EXPECT_EQ(parameters.sccs_config.boundary, parameters.pcc_boundary);
    EXPECT_DOUBLE_EQ(parameters.sccs_config.cavity.epsilon_bulk, 1.0);
    EXPECT_EQ(parameters.sccs_config.mixing_method, "pulay");
    EXPECT_DOUBLE_EQ(parameters.sccs_config.mixing, 0.2);
    EXPECT_DOUBLE_EQ(parameters.start_drho, 0.01);
    EXPECT_EQ(parameters.start_nmax, 12);
    EXPECT_EQ(parameters.debug, 2);
}

TEST(SurchemInput, PreservesLegacyParametersAndOrdinaryVacuum)
{
    Input_para input;
    UnitCell cell;
    EXPECT_FALSE(input.uses_surchem_correction());
    input.imp_sol = 1;
    input.eb_k = 80.0;
    input.tau = 0.00002;
    const SurchemParameters parameters = ModuleSurchem::make_parameters(input, cell, 0.0, false, 1);
    EXPECT_TRUE(parameters.use_legacy_solvent);
    EXPECT_FALSE(parameters.use_sccs);
    EXPECT_EQ(parameters.pcc_boundary, ModuleSccs::Boundary::Periodic);
    EXPECT_DOUBLE_EQ(parameters.eb_k, input.eb_k);
    EXPECT_DOUBLE_EQ(parameters.tau, input.tau);
}

TEST(SurchemInput, UltrasoftPseudopotentialsWarnWithoutRejectingSccsOrPcc)
{
    UnitCell cell;
    for (const int model : {0, 2})
    {
        for (const std::string boundary : {"pcc_0d", "pcc_2d"})
        {
            Input_para input;
            input.imp_sol = model;
            input.assume_isolated = boundary;
            const SurchemParameters parameters = ModuleSurchem::make_parameters(input, cell, 0.0, true, 1);
            EXPECT_EQ(parameters.use_sccs, model == 2);
            EXPECT_NE(parameters.pcc_boundary, ModuleSccs::Boundary::Periodic);
        }
    }
    Input_para periodic;
    periodic.imp_sol = 2;
    const SurchemParameters parameters = ModuleSurchem::make_parameters(periodic, cell, 0.0, true, 1);
    EXPECT_TRUE(parameters.use_sccs);
    EXPECT_EQ(parameters.pcc_boundary, ModuleSccs::Boundary::Periodic);
}

TEST(SurchemInput, PresetsOverrideOnlyPhysicalParameters)
{
    struct PresetValues
    {
        const char* name;
        double epsilon;
        double rho_min;
        double rho_max;
        double gamma;
        double pressure;
    };
    const PresetValues values[] = {
        {"custom", 12.0, 0.001, 0.02, 2.0, 0.2},
        {"vacuum", 1.0, 1.0e-4, 5.0e-3, 0.0, 0.0},
        {"water-neutral", 78.3, 1.0e-4, 5.0e-3, 47.9, -0.36},
        {"water-cation", 78.3, 2.0e-4, 3.5e-3, 5.0, 0.125},
        {"water-anion", 78.3, 2.4e-3, 1.55e-2, 0.0, 0.45}};
    UnitCell cell;
    for (const PresetValues& preset : values)
    {
        Input_para input;
        input.imp_sol = 2;
        input.sccs_preset = preset.name;
        input.sccs_epsilon = 12.0;
        input.sccs_rho_min = 0.001;
        input.sccs_rho_max = 0.02;
        input.sccs_gamma = 2.0;
        input.sccs_pressure = 0.2;
        input.sccs_mixing_type = "anderson";
        input.sccs_mixing = 0.25;
        input.sccs_surface_eta = 2.0e-8;
        const SurchemParameters parameters = ModuleSurchem::make_parameters(input, cell, 0.0, false, 1);
        const double gamma = ModuleSccs::dyn_per_cm_to_hartree_per_bohr2(preset.gamma);
        const double pressure = ModuleSccs::gpa_to_hartree_per_bohr3(preset.pressure);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.cavity.epsilon_bulk, preset.epsilon);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.cavity.density_min, preset.rho_min);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.cavity.density_max, preset.rho_max);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.surface_tension, gamma);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.pressure, pressure);
        EXPECT_EQ(parameters.sccs_config.mixing_method, "anderson");
        EXPECT_DOUBLE_EQ(parameters.sccs_config.mixing, 0.25);
        EXPECT_DOUBLE_EQ(parameters.sccs_config.surface_regularization, 2.0e-8);
    }
}

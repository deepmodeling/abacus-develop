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

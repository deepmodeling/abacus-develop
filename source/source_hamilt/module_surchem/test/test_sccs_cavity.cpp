#include "../sccs/sccs_cavity.h"

#include "gtest/gtest.h"

#include <cmath>
#include <stdexcept>

namespace
{

ModuleSccs::CavityParameters water_parameters()
{
    ModuleSccs::CavityParameters parameters;
    parameters.density_min = 1.0e-4;
    parameters.density_max = 5.0e-3;
    parameters.epsilon_bulk = 78.3;
    return parameters;
}

TEST(SccsCavity, ReachesBulkAndSoluteLimits)
{
    const ModuleSccs::CavityParameters parameters = water_parameters();
    const ModuleSccs::CavityPoint solvent = ModuleSccs::evaluate_cavity(0.0, parameters);
    const ModuleSccs::CavityPoint solute = ModuleSccs::evaluate_cavity(0.1, parameters);

    EXPECT_DOUBLE_EQ(solvent.solute, 0.0);
    EXPECT_DOUBLE_EQ(solvent.epsilon, parameters.epsilon_bulk);
    EXPECT_DOUBLE_EQ(solvent.dsolute_drho, 0.0);
    EXPECT_DOUBLE_EQ(solvent.depsilon_drho, 0.0);
    EXPECT_DOUBLE_EQ(solute.solute, 1.0);
    EXPECT_DOUBLE_EQ(solute.epsilon, 1.0);
    EXPECT_DOUBLE_EQ(solute.dsolute_drho, 0.0);
    EXPECT_DOUBLE_EQ(solute.depsilon_drho, 0.0);
}

TEST(SccsCavity, IsContinuousAtDensityThresholds)
{
    const ModuleSccs::CavityParameters parameters = water_parameters();
    const double relative_step = 1.0e-8;
    const ModuleSccs::CavityPoint below
        = ModuleSccs::evaluate_cavity(parameters.density_min * (1.0 + relative_step), parameters);
    const ModuleSccs::CavityPoint above
        = ModuleSccs::evaluate_cavity(parameters.density_max * (1.0 - relative_step), parameters);

    EXPECT_NEAR(below.solute, 0.0, 1.0e-14);
    EXPECT_NEAR(below.epsilon, parameters.epsilon_bulk, 1.0e-12);
    EXPECT_NEAR(below.dsolute_drho, 0.0, 1.0e-10);
    EXPECT_NEAR(above.solute, 1.0, 1.0e-14);
    EXPECT_NEAR(above.epsilon, 1.0, 1.0e-14);
    EXPECT_NEAR(above.dsolute_drho, 0.0, 1.0e-12);
}

TEST(SccsCavity, AnalyticDerivativesMatchCentralDifferences)
{
    const ModuleSccs::CavityParameters parameters = water_parameters();
    const double density = std::sqrt(parameters.density_min * parameters.density_max);
    const double step = density * 1.0e-5;
    const ModuleSccs::CavityPoint center = ModuleSccs::evaluate_cavity(density, parameters);
    const ModuleSccs::CavityPoint plus = ModuleSccs::evaluate_cavity(density + step, parameters);
    const ModuleSccs::CavityPoint minus = ModuleSccs::evaluate_cavity(density - step, parameters);

    const double finite_solute = (plus.solute - minus.solute) / (2.0 * step);
    const double finite_epsilon = (plus.epsilon - minus.epsilon) / (2.0 * step);
    EXPECT_NEAR(center.dsolute_drho, finite_solute, std::abs(finite_solute) * 1.0e-9);
    EXPECT_NEAR(center.depsilon_drho, finite_epsilon, std::abs(finite_epsilon) * 1.0e-9);
    EXPECT_GT(center.dsolute_drho, 0.0);
    EXPECT_LT(center.depsilon_drho, 0.0);
}

TEST(SccsCavity, VacuumPermittivityKeepsTheCavityDefinition)
{
    ModuleSccs::CavityParameters parameters = water_parameters();
    parameters.epsilon_bulk = 1.0;
    const double density = std::sqrt(parameters.density_min * parameters.density_max);
    const ModuleSccs::CavityPoint result = ModuleSccs::evaluate_cavity(density, parameters);

    EXPECT_GT(result.solute, 0.0);
    EXPECT_LT(result.solute, 1.0);
    EXPECT_DOUBLE_EQ(result.epsilon, 1.0);
    EXPECT_DOUBLE_EQ(result.depsilon_drho, 0.0);
}

TEST(SccsCavity, RejectsInvalidInputs)
{
    ModuleSccs::CavityParameters parameters = water_parameters();
    EXPECT_THROW(ModuleSccs::evaluate_cavity(std::nan(""), parameters), std::domain_error);
    parameters.density_min = parameters.density_max;
    EXPECT_THROW(ModuleSccs::evaluate_cavity(1.0e-3, parameters), std::invalid_argument);
}

TEST(SccsCavity, TreatsNegativeGridRingingAsBulkSolvent)
{
    const ModuleSccs::CavityParameters parameters = water_parameters();
    const ModuleSccs::CavityPoint result
        = ModuleSccs::evaluate_cavity(-1.0e-6, parameters);
    EXPECT_DOUBLE_EQ(result.solute, 0.0);
    EXPECT_DOUBLE_EQ(result.epsilon, parameters.epsilon_bulk);
    EXPECT_DOUBLE_EQ(result.dsolute_drho, 0.0);
    EXPECT_DOUBLE_EQ(result.depsilon_drho, 0.0);
}

} // namespace

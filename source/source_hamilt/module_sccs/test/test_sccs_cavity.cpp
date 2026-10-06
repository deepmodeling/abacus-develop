#include "../sccs_cavity.h"
#include <gtest/gtest.h>

#include <cmath>

namespace
{
class SccsCavityTest : public ::testing::Test
{
protected:
    ModuleSccs::CavityParameters parameters;

    void SetUp() override
    {
        parameters.density_min = 1e-4;
        parameters.density_max = 5e-3;
        parameters.epsilon_bulk = 78.3;
    }
};

TEST_F(SccsCavityTest, BulkAndSoluteLimits)
{
    const double bulk_densities[] = {-1e-6, 0.0, 1e-4};
    for (double density : bulk_densities)
    {
        ModuleSccs::CavityPoint point;
        point = ModuleSccs::evaluate_cavity(density, parameters);
        EXPECT_DOUBLE_EQ(point.solute, 0.0);
        EXPECT_DOUBLE_EQ(point.epsilon, 78.3);
        EXPECT_DOUBLE_EQ(point.dsolute_drho, 0.0);
        EXPECT_DOUBLE_EQ(point.depsilon_drho, 0.0);
    }
    const double solute_densities[] = {5e-3, 0.1};
    for (double density : solute_densities)
    {
        ModuleSccs::CavityPoint point;
        point = ModuleSccs::evaluate_cavity(density, parameters);
        EXPECT_DOUBLE_EQ(point.solute, 1.0);
        EXPECT_DOUBLE_EQ(point.epsilon, 1.0);
        EXPECT_DOUBLE_EQ(point.dsolute_drho, 0.0);
        EXPECT_DOUBLE_EQ(point.depsilon_drho, 0.0);
    }
}

TEST_F(SccsCavityTest, LogarithmicMidpoint)
{
    const double product = parameters.density_min * parameters.density_max;
    const double density = std::sqrt(product);
    const double ratio = parameters.density_max / parameters.density_min;
    const double width = std::log(ratio);
    const double epsilon = std::sqrt(parameters.epsilon_bulk);
    const double log_epsilon = std::log(parameters.epsilon_bulk);
    const double dsolute = 2.0 / (width * density);
    const double depsilon = -epsilon * log_epsilon * dsolute;
    ModuleSccs::CavityPoint point;
    point = ModuleSccs::evaluate_cavity(density, parameters);
    EXPECT_NEAR(point.solute, 0.5, 1e-14);
    EXPECT_NEAR(point.epsilon, epsilon, 1e-13);
    EXPECT_NEAR(point.dsolute_drho, dsolute, 1e-10);
    EXPECT_NEAR(point.depsilon_drho, depsilon, 1e-9);
}

TEST_F(SccsCavityTest, DerivativesAcrossTransition)
{
    const double fractions[] = {0.2, 0.5, 0.8};
    const double ratio = parameters.density_max / parameters.density_min;
    for (double fraction : fractions)
    {
        const double scale = std::pow(ratio, fraction);
        const double density = parameters.density_min * scale;
        const double step = density * 1e-5;
        const double lower_density = density - step;
        const double upper_density = density + step;
        ModuleSccs::CavityPoint point;
        ModuleSccs::CavityPoint lower;
        ModuleSccs::CavityPoint upper;
        point = ModuleSccs::evaluate_cavity(density, parameters);
        lower = ModuleSccs::evaluate_cavity(lower_density, parameters);
        upper = ModuleSccs::evaluate_cavity(upper_density, parameters);
        const double solute_fd = (upper.solute - lower.solute) / (2.0 * step);
        const double epsilon_fd = (upper.epsilon - lower.epsilon) / (2.0 * step);
        const double solute_ratio = solute_fd / point.dsolute_drho;
        const double epsilon_ratio = epsilon_fd / point.depsilon_drho;
        EXPECT_NEAR(solute_ratio, 1.0, 1e-8);
        EXPECT_NEAR(epsilon_ratio, 1.0, 1e-8);
        EXPECT_GT(point.dsolute_drho, 0.0);
        EXPECT_LT(point.depsilon_drho, 0.0);
    }
}

TEST_F(SccsCavityTest, SmoothThresholds)
{
    const double near_min = parameters.density_min * (1.0 + 1e-8);
    const double near_max = parameters.density_max * (1.0 - 1e-8);
    ModuleSccs::CavityPoint lower;
    ModuleSccs::CavityPoint upper;
    lower = ModuleSccs::evaluate_cavity(near_min, parameters);
    upper = ModuleSccs::evaluate_cavity(near_max, parameters);
    EXPECT_NEAR(lower.solute, 0.0, 1e-14);
    EXPECT_NEAR(upper.solute, 1.0, 1e-14);
    EXPECT_NEAR(lower.epsilon, 78.3, 1e-12);
    EXPECT_NEAR(upper.epsilon, 1.0, 1e-12);
    EXPECT_NEAR(lower.dsolute_drho, 0.0, 1e-10);
    EXPECT_NEAR(upper.dsolute_drho, 0.0, 1e-10);
    EXPECT_NEAR(lower.depsilon_drho, 0.0, 1e-8);
    EXPECT_NEAR(upper.depsilon_drho, 0.0, 1e-8);
}

TEST_F(SccsCavityTest, VacuumHasNoDielectricResponse)
{
    parameters.epsilon_bulk = 1.0;
    ModuleSccs::CavityPoint point;
    point = ModuleSccs::evaluate_cavity(1e-3, parameters);
    EXPECT_GT(point.solute, 0.0);
    EXPECT_LT(point.solute, 1.0);
    EXPECT_DOUBLE_EQ(point.epsilon, 1.0);
    EXPECT_DOUBLE_EQ(point.depsilon_drho, 0.0);
}
} // namespace

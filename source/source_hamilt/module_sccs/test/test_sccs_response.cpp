#include "sccs_test.h"
#include "../sccs_response.h"
#include "../sccs_parameters.h"

#include "source_base/parallel_reduce.h"

#include <algorithm>

using SccsResponseTest = SccsTest::PwTest;

TEST_F(SccsResponseTest, VacuumAndUniformDielectricAnalyticLimits)
{
    const double dielectrics[] = {1.0, 5.0};
    const std::vector<double> charge = cosine_mode(0);
    const std::vector<double> density(basis.nrxx, 0.0);
    const std::vector<double> cold_start;
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    for (double epsilon : dielectrics)
    {
        cavity.epsilon_bulk = epsilon;
        ModuleSccs::SccsResponse response;
        ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                                   basis, tpiba, response, error)) << error;
        const double kernel = ModuleBase::FOUR_PI / (epsilon * tpiba * tpiba);
        EXPECT_EQ(response.polarization.iterations, 1);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const double expected = kernel * charge[ir];
            EXPECT_NEAR(response.polarization.potential[ir], expected, 1e-11);
            EXPECT_DOUBLE_EQ(response.cavity_potential[ir], 0.0);
        }
    }
}

TEST_F(SccsResponseTest, ManufacturedNonuniformDielectricSolution)
{
    const std::vector<double> mode_x = cosine_mode(0);
    const std::vector<double> mode_y = cosine_mode(1);
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 78.3;
    std::vector<double> density(basis.nrxx);
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 1e-3 + 2e-5 * mode_x[ir];
        ModuleSccs::CavityPoint point;
        ASSERT_TRUE(ModuleSccs::evaluate_cavity(density[ir], cavity, point, error));
        // v=cos(ky), eps=eps(x): -div(eps grad v)/(4 pi)=eps k^2 cos(ky)/(4 pi).
        charge[ir] = point.epsilon * tpiba * tpiba * mode_y[ir] / ModuleBase::FOUR_PI;
    }
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    const std::vector<double> cold_start;
    ModuleSccs::SccsResponse response;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, response, error)) << error;
    double maximum_error = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double difference = response.polarization.potential[ir] - mode_y[ir];
        const double magnitude = std::abs(difference);
        maximum_error = std::max(maximum_error, magnitude);
        EXPECT_NEAR(response.polarization.gradient[ir].x, 0.0, 1e-9);
        EXPECT_NEAR(response.polarization.gradient[ir].z, 0.0, 1e-9);
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, maximum_error);
    EXPECT_LT(maximum_error, 1e-9);
    EXPECT_GT(response.polarization.iterations, 1);
}

TEST_F(SccsResponseTest, ConvergenceFailureAndExplicitWarmStart)
{
    const std::vector<double> mode_x = cosine_mode(0);
    std::vector<double> charge = cosine_mode(1);
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 0.009 + 0.008 * mode_x[ir];
        charge[ir] *= 1e-3;
    }
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 0.0024;
    cavity.density_max = 0.0155;
    cavity.epsilon_bulk = 78.3;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-11;
    solver.tolerance_max = 1e-10;
    const std::vector<double> cold_start;
    ModuleSccs::SccsResponse response;
    solver.max_iterations = 1;
    response.polarization.potential.assign(1, 12.0);
    EXPECT_FALSE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, response, error));
    EXPECT_NE(error.find("iteration limit"), std::string::npos);
    EXPECT_DOUBLE_EQ(response.polarization.potential[0], 12.0);
    solver.max_iterations = 200;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, response, error)) << error;
    EXPECT_LE(response.polarization.residual_rms, solver.tolerance_rms);
    EXPECT_LE(response.polarization.residual_max, solver.tolerance_max);
    for (double& value : density)
    {
        value *= 1.001;
    }
    ModuleSccs::SccsResponse cold;
    ModuleSccs::SccsResponse warm;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, cold, error)) << error;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, response.restart_potential,
                                               basis, tpiba, warm, error)) << error;
    EXPECT_TRUE(warm.polarization.warm_started);
    EXPECT_LT(warm.polarization.iterations, cold.polarization.iterations);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(warm.polarization.potential[ir], cold.polarization.potential[ir], 1e-9);
    }
}

TEST_F(SccsResponseTest, RankLocalInvalidParametersAndRestartReturnCollectively)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    const std::vector<double> density(basis.nrxx, 0.0);
    const std::vector<double> charge = cosine_mode(0);
    const std::vector<double> cold_start;
    ModuleSccs::SccsResponse response;
    response.polarization.potential.assign(1, 12.0);
    if (basis.poolrank == 0)
    {
        solver.max_iterations = 0;
    }
    EXPECT_FALSE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, response, error));
    EXPECT_DOUBLE_EQ(response.polarization.potential[0], 12.0);
    solver.max_iterations = 200;
    std::vector<double> initial(basis.nrxx, 0.0);
    if (basis.poolrank == 0)
    {
        initial.pop_back();
    }
    EXPECT_FALSE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, initial,
                                               basis, tpiba, response, error));
    EXPECT_DOUBLE_EQ(response.polarization.potential[0], 12.0);
}

TEST_F(SccsResponseTest, RequiresBothResidualTolerances)
{
    const std::vector<double> mode_x = cosine_mode(0);
    std::vector<double> charge = cosine_mode(1);
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 0.009 + 0.008 * mode_x[ir];
        charge[ir] *= 1e-3;
    }
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 0.0024;
    cavity.density_max = 0.0155;
    cavity.epsilon_bulk = 78.3;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-5;
    solver.tolerance_max = 1e-5;
    const std::vector<double> cold_start;
    ModuleSccs::SccsResponse loose;
    ModuleSccs::SccsResponse tight_maximum;
    ModuleSccs::SccsResponse tight_rms;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, loose, error)) << error;
    solver.tolerance_max = 1e-11;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, tight_maximum, error)) << error;
    solver.tolerance_rms = 1e-11;
    solver.tolerance_max = 1.0;
    ASSERT_TRUE(ModuleSccs::solve_sccs_response(density, charge, cavity, solver, cold_start,
                                               basis, tpiba, tight_rms, error)) << error;
    EXPECT_LE(tight_maximum.polarization.residual_max, 1e-11);
    EXPECT_LE(tight_rms.polarization.residual_rms, 1e-11);
    EXPECT_GT(tight_maximum.polarization.iterations, loose.polarization.iterations);
    EXPECT_GT(tight_rms.polarization.iterations, loose.polarization.iterations);
}

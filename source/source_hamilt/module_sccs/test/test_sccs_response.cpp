#include "sccs_test.h"
#include "../sccs_functional.h"
#include "../sccs_response.h"
#include "../sccs_parameters.h"

#include "source_base/parallel_reduce.h"

#include <algorithm>
#include <cmath>
#include <vector>

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
        ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                                   basis, tpiba, response);
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
        const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(density[ir], cavity);
        // v=cos(ky), eps=eps(x): -div(eps grad v)/(4 pi)=eps k^2 cos(ky)/(4 pi).
        charge[ir] = point.epsilon * tpiba * tpiba * mode_y[ir] / ModuleBase::FOUR_PI;
    }
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-12;
    solver.tolerance_max = 1e-12;
    const std::vector<double> cold_start;
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                               basis, tpiba, response);
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

TEST_F(SccsResponseTest, ConvergenceFailureBothTolerancesAndExplicitWarmStart)
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
    // A missed tolerance stops the run; a forked death test cannot share the
    // pool collectives of a multi-rank run.
    if (SccsTest::pool_size == 1)
    {
        solver.max_iterations = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleSccs::solve_response(density, charge, cavity, solver, cold_start, basis, tpiba,
                                               response),
                    ::testing::ExitedWithCode(1), "");
        testing::internal::GetCapturedStdout();
    }
    solver.max_iterations = 200;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start, basis, tpiba, response);
    EXPECT_LE(response.polarization.residual_rms, solver.tolerance_rms);
    EXPECT_LE(response.polarization.residual_max, solver.tolerance_max);
    for (double& value : density)
    {
        value *= 1.001;
    }
    ModuleSccs::SccsResponse cold;
    ModuleSccs::SccsResponse warm;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                               basis, tpiba, cold);
    ModuleSccs::solve_response(density, charge, cavity, solver, response.restart_potential,
                               basis, tpiba, warm);
    EXPECT_TRUE(warm.polarization.warm_started);
    EXPECT_LT(warm.polarization.iterations, cold.polarization.iterations);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(warm.polarization.potential[ir], cold.polarization.potential[ir], 1e-9);
    }

    // Each of the two tolerances alone drives the iteration count.
    solver.tolerance_rms = 1e-5;
    solver.tolerance_max = 1e-5;
    ModuleSccs::SccsResponse loose;
    ModuleSccs::SccsResponse tight_maximum;
    ModuleSccs::SccsResponse tight_rms;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                               basis, tpiba, loose);
    solver.tolerance_max = 1e-11;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                               basis, tpiba, tight_maximum);
    solver.tolerance_rms = 1e-11;
    solver.tolerance_max = 1.0;
    ModuleSccs::solve_response(density, charge, cavity, solver, cold_start,
                               basis, tpiba, tight_rms);
    EXPECT_LE(tight_maximum.polarization.residual_max, 1e-11);
    EXPECT_LE(tight_rms.polarization.residual_rms, 1e-11);
    EXPECT_GT(tight_maximum.polarization.iterations, loose.polarization.iterations);
    EXPECT_GT(tight_rms.polarization.iterations, loose.polarization.iterations);
}

namespace
{
// cos(2 pi harmonic x / L) on a basis' local real-space grid.
std::vector<double> x_harmonic(const ModulePW::PW_Basis& basis, int harmonic)
{
    std::vector<double> values(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const double angle = ModuleBase::TWO_PI * harmonic * ix / basis.nx;
        values[ir] = std::cos(angle);
    }
    return values;
}
} // namespace

// The cavity potential -eps'|grad v|^2/(8 pi) against a central finite
// difference of the reaction energy along 1 + cos(3 g x). The base density is
// uniform inside the switching interval, so v = 0.1 cos(6 g x) solves the base
// response exactly; the 80 Ry grid resolves |grad v|^2.
TEST(SccsResponseRefinement, CavityDerivativeMatchesFiniteDifferenceOnAResolvedGrid)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(SccsTest::pool_size, SccsTest::pool_rank, POOL_WORLD);
#endif
    const double length = 10.0;
    const double gridecut = 80.0;
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    basis.initgrids(length, lattice, gridecut);
    basis.initparameters(false, gridecut, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    const double tpiba = ModuleBase::TWO_PI / length;

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1e-4;
    config.cavity.density_max = 5e-3;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1e-8;
    // The CG stalls at a uniform residual (~2e-9) outside the range of the
    // zero-mean preconditioner; it does not change the energy.
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-8;
    solver.tolerance_max = 1e-8;
    const double base_density = 7e-4;
    const double base_epsilon = ModuleSccs::evaluate_cavity(base_density, config.cavity).epsilon;
    const double wavevector = 6.0 * tpiba;
    const std::vector<double> mode6 = x_harmonic(basis, 6);
    const std::vector<double> mode3 = x_harmonic(basis, 3);
    std::vector<double> charge(basis.nrxx);
    std::vector<double> direction(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        charge[ir] = 0.1 * base_epsilon * wavevector * wavevector * mode6[ir] / ModuleBase::FOUR_PI;
        direction[ir] = 1.0 + mode3[ir];
    }
    const std::vector<double> cold_start;
    auto reaction_energy = [&](double step, ModuleSccs::SccsResponse& response) -> double {
        std::vector<double> density(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            density[ir] = base_density + step * direction[ir];
        }
        ModuleSccs::solve_response(density, charge, config.cavity, solver, cold_start, basis, tpiba, response);
        ModuleSccs::FunctionalResult functional;
        ModuleSccs::eval_functional(charge, response, config, basis, tpiba, functional);
        return functional.reaction_energy;
    };
    ModuleSccs::SccsResponse center;
    ModuleSccs::SccsResponse plus;
    ModuleSccs::SccsResponse minus;
    reaction_energy(0.0, center);
    const double step = 1e-7;
    const double backward_step = -step;
    const double energy_plus = reaction_energy(step, plus);
    const double energy_minus = reaction_energy(backward_step, minus);
    const double finite_difference = (energy_plus - energy_minus) / (2.0 * step);
    double analytic = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        analytic += center.cavity_potential[ir] * direction[ir];
    }
    Parallel_Reduce::reduce_pool(analytic);
    analytic *= basis.omega / basis.nxyz;
    const double tolerance = 1e-6 * std::abs(finite_difference);
    EXPECT_NEAR(analytic, finite_difference, tolerance);
}

#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_periodic.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

TEST(SccsPeriodic, ContinuumSourceScreensUniformDielectricAndIncludesInterfaceField)
{
    const std::vector<double> charge = {1.0, -2.0};
    ModuleSccs::PeriodicSccsResult response;
    response.epsilon.assign(2, 5.0);
    response.grad_log_epsilon.resize(2);
    response.polarization.field.gradient.resize(2);
    const std::vector<double> uniform = ModuleSccs::continuum_polarization_charge(charge, response);
    EXPECT_NEAR(uniform[0], -0.8, 1.0e-15);
    EXPECT_NEAR(uniform[1], 1.6, 1.0e-15);
    response.grad_log_epsilon[0].y = ModuleBase::FOUR_PI;
    response.polarization.field.gradient[0].y = 0.5;
    const std::vector<double> interface = ModuleSccs::continuum_polarization_charge(charge, response);
    EXPECT_NEAR(interface[0], -0.3, 1.0e-15);
    EXPECT_NEAR(interface[1], uniform[1], 1.0e-15);
    response.epsilon[0] = 0.0;
    EXPECT_THROW(ModuleSccs::continuum_polarization_charge(charge, response), std::domain_error);
    response.epsilon.pop_back();
    EXPECT_THROW(ModuleSccs::continuum_polarization_charge(charge, response), std::invalid_argument);
}

TEST(SccsPeriodic, UniformDielectricScreensSingleFourierShell)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    basis.initgrids(10.0, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    std::vector<std::complex<double>> charge_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const double gx = basis.gdirect[ig].x;
        const double gy = basis.gdirect[ig].y;
        const double gz = basis.gdirect[ig].z;
        if (std::abs(std::abs(gx) - 1.0) < 1.0e-12 && std::abs(gy) < 1.0e-12
            && std::abs(gz) < 1.0e-12)
        {
            charge_g[ig] = 0.5;
        }
    }
    std::vector<double> solute_charge(basis.nrxx);
    basis.recip2real(charge_g.data(), solute_charge.data());

    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1.0e-4;
    cavity.density_max = 5.0e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 100;
    solver.mixing = 0.7;
    solver.tolerance_rms = 1.0e-12;
    solver.tolerance_max = 1.0e-12;

    const std::vector<double> cavity_density(basis.nrxx, 0.0);
    const ModuleSccs::PeriodicSccsResult result
        = ModuleSccs::solve_periodic_sccs(cavity_density,
                                          solute_charge,
                                          cavity,
                                          solver,
                                          std::vector<double>(),
                                          basis,
                                          ModuleBase::TWO_PI / 10.0,
                                          1);

    EXPECT_EQ(result.polarization.status, ModuleSccs::PolarizationStatus::Converged);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(result.epsilon[ir], 5.0);
        EXPECT_NEAR(result.polarization.polarization_charge[ir], -0.8 * solute_charge[ir], 2.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].x, 0.0, 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].z, 0.0, 1.0e-12);
    }
}

TEST(SccsPeriodic, ChainGradientMatchesAnalyticDensityModeAcrossCavityEdges)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    const double tpiba = ModuleBase::TWO_PI / length;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 0.0024;
    cavity.density_max = 0.0155;
    cavity.epsilon_bulk = 78.3;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 100;
    solver.mixing = 0.5;
    solver.tolerance_rms = 1.0e-13;
    solver.tolerance_max = 1.0e-11;
    std::vector<double> density(basis.nrxx);
    std::vector<double> expected_gradient(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const double phase = ModuleBase::TWO_PI * ix / basis.nx;
        density[ir] = 0.009 + 0.008 * std::cos(phase);
        const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(density[ir], cavity);
        expected_gradient[ir] = -0.008 * tpiba * std::sin(phase)
                                * point.depsilon_drho / point.epsilon;
    }
    const std::vector<double> charge(basis.nrxx, 0.0);
    const std::vector<double> initial;
    const ModuleSccs::PeriodicSccsResult result
        = ModuleSccs::solve_periodic_sccs(density, charge, cavity, solver, initial, basis, tpiba, 1);
    ASSERT_EQ(result.polarization.status, ModuleSccs::PolarizationStatus::Converged);
    ASSERT_EQ(result.density_gradient.size(), density.size());
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(result.grad_log_epsilon[ir].x, expected_gradient[ir], 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].z, 0.0, 1.0e-12);
        const int ix = ir / (basis.ny * basis.nplane);
        const double phase = ModuleBase::TWO_PI * ix / basis.nx;
        const double expected_density_gradient = -0.008 * tpiba * std::sin(phase);
        EXPECT_NEAR(result.density_gradient[ir].x, expected_density_gradient, 1.0e-12);
        EXPECT_NEAR(result.density_gradient[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.density_gradient[ir].z, 0.0, 1.0e-12);
    }
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}

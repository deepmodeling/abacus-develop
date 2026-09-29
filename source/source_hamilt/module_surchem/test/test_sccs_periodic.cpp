#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_periodic.h"
#include "../sccs/sccs_pw_coulomb.h"
#include "../common/charge_reduction.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

// The SCCS driver's periodic path: plain Coulomb preconditioner, pool reduction.
ModuleSccs::PeriodicSccsResult solve_periodic(const std::vector<double>& cavity_density,
                                              const std::vector<double>& solute_charge,
                                              const ModuleSccs::CavityParameters& cavity,
                                              const ModuleSccs::PolarizationSolverParameters& solver,
                                              const std::vector<double>& initial_potential,
                                              const ModulePW::PW_Basis& basis,
                                              const double tpiba)
{
    const ModuleSccs::PeriodicCoulombOperator coulomb(basis, tpiba);
    const ModuleSurchem::PoolChargeReduction reduction(1);
    return ModuleSccs::solve_chain_sccs_response(cavity_density,
                                                 solute_charge,
                                                 cavity,
                                                 solver,
                                                 initial_potential,
                                                 basis,
                                                 tpiba,
                                                 coulomb,
                                                 reduction);
}

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
    solver.tolerance_rms = 1.0e-12;
    solver.tolerance_max = 1.0e-12;

    const std::vector<double> cavity_density(basis.nrxx, 0.0);
    const std::vector<double> cold_start;
    const double tpiba = ModuleBase::TWO_PI / 10.0;
    const ModuleSccs::PeriodicSccsResult result
        = solve_periodic(cavity_density, solute_charge, cavity, solver, cold_start, basis, tpiba);

    EXPECT_DOUBLE_EQ(result.far_field_polarization_charge, 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(result.epsilon[ir], 5.0);
        EXPECT_NEAR(result.polarization.polarization_charge[ir], -0.8 * solute_charge[ir], 2.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].x, 0.0, 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.grad_log_epsilon[ir].z, 0.0, 1.0e-12);
    }
}

// A y charge mode inside an x-modulated dielectric needs several CG steps.
class SqrtCgFixture : public testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        basis.initmpi(1, 0, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0);
        basis.initgrids(length, lattice, 20.0);
        basis.initparameters(false, 20.0, 1, false);
        basis.setuptransform();
        basis.collect_local_pw();
        cavity.density_min = 0.0024;
        cavity.density_max = 0.0155;
        cavity.epsilon_bulk = 78.3;
        solver.max_iterations = 200;
        density.resize(basis.nrxx);
        charge.resize(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const int ix = ir / (basis.ny * basis.nplane);
            const int iy = ir / basis.nplane - ix * basis.ny;
            density[ir] = 0.009 + 0.008 * std::cos(ModuleBase::TWO_PI * ix / basis.nx);
            charge[ir] = 1.0e-3 * std::cos(ModuleBase::TWO_PI * iy / basis.ny);
        }
    }

    ModuleSccs::PeriodicSccsResult solve_from(const std::vector<double>& initial_potential) const
    {
        const double tpiba = ModuleBase::TWO_PI / length;
        return solve_periodic(density, charge, cavity, solver, initial_potential, basis, tpiba);
    }

    ModuleSccs::PeriodicSccsResult solve() const
    {
        const std::vector<double> cold_start;
        return solve_from(cold_start);
    }

    const double length = 10.0;
    ModulePW::PW_Basis basis{"cpu", "double"};
    ModuleSccs::CavityParameters cavity;
    ModuleSccs::PolarizationSolverParameters solver;
    std::vector<double> density;
    std::vector<double> charge;
};

TEST_F(SqrtCgFixture, StopsOnlyWhenRmsAndMaximumResidualsPass)
{
    // The initial residual has maximum 1e-3, so every case needs CG steps.
    solver.tolerance_rms = 1.0e-5;
    solver.tolerance_max = 1.0e-5;
    const ModuleSccs::PeriodicSccsResult loose = solve();
    solver.tolerance_max = 1.0e-11;
    const ModuleSccs::PeriodicSccsResult tight_maximum = solve();
    solver.tolerance_rms = 1.0e-11;
    solver.tolerance_max = 1.0;
    const ModuleSccs::PeriodicSccsResult tight_rms = solve();

    EXPECT_LE(loose.polarization.residual_rms, 1.0e-5);
    EXPECT_LE(loose.polarization.residual_max, 1.0e-5);
    EXPECT_LE(tight_maximum.polarization.residual_max, 1.0e-11);
    EXPECT_LE(tight_rms.polarization.residual_rms, 1.0e-11);
    EXPECT_GT(loose.polarization.iterations, 0);
    EXPECT_GT(tight_maximum.polarization.iterations, loose.polarization.iterations);
    EXPECT_GT(tight_rms.polarization.iterations, loose.polarization.iterations);
}

TEST_F(SqrtCgFixture, VerifiesPreconditionedFixedPointOnlyOnRequest)
{
    solver.tolerance_rms = 1.0e-11;
    solver.tolerance_max = 1.0e-10;
    const ModuleSccs::PeriodicSccsResult unchecked = solve();
    EXPECT_FALSE(unchecked.polarization.fixed_point_checked);

    solver.check_fixed_point = true;
    const ModuleSccs::PeriodicSccsResult checked = solve();
    ASSERT_TRUE(checked.polarization.fixed_point_checked);
    EXPECT_LT(checked.polarization.fixed_point_defect_rms, 1.0e-8);
    EXPECT_LT(checked.polarization.fixed_point_defect_max, 1.0e-7);
    EXPECT_EQ(checked.polarization.iterations, unchecked.polarization.iterations);
}

TEST_F(SqrtCgFixture, WarmStartFromPreviousPotentialReachesSameSolutionFaster)
{
    solver.tolerance_rms = 1.0e-11;
    solver.tolerance_max = 1.0e-10;
    const ModuleSccs::PeriodicSccsResult cold = solve();
    ASSERT_GT(cold.polarization.iterations, 1);
    EXPECT_FALSE(cold.polarization.warm_started);

    // A slightly changed cavity mimics the next SCF step.
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] *= 1.001;
    }
    const ModuleSccs::PeriodicSccsResult reference = solve();
    ASSERT_EQ(cold.restart_potential.size(), static_cast<std::size_t>(basis.nrxx));
    const ModuleSccs::PeriodicSccsResult warm = solve_from(cold.restart_potential);
    EXPECT_TRUE(warm.polarization.warm_started);
    EXPECT_LT(warm.polarization.iterations, reference.polarization.iterations);
    EXPECT_LE(warm.polarization.residual_rms, 1.0e-11);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(warm.polarization.field.potential[ir],
                    reference.polarization.field.potential[ir],
                    1.0e-9);
    }
}

TEST_F(SqrtCgFixture, RejectsWarmStartWorseThanColdStart)
{
    solver.tolerance_rms = 1.0e-11;
    solver.tolerance_max = 1.0e-10;
    const ModuleSccs::PeriodicSccsResult cold = solve();
    std::vector<double> poor_guess(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        poor_guess[ir] = 1.0e3 * std::cos(0.37 * ir);
    }
    const ModuleSccs::PeriodicSccsResult rejected = solve_from(poor_guess);
    EXPECT_FALSE(rejected.polarization.warm_started);
    EXPECT_EQ(rejected.polarization.iterations, cold.polarization.iterations);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(rejected.polarization.field.potential[ir],
                         cold.polarization.field.potential[ir]);
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
        = solve_periodic(density, charge, cavity, solver, initial, basis, tpiba);
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
    // Error-path tests throw inside timed functions and leave their
    // ModuleBase::timer entries running; production turns these exceptions
    // into WARNING_QUIT, so the timers are not under test here.
    ModuleBase::timer::disable();
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}

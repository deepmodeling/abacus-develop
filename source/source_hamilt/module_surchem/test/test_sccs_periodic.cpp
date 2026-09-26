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

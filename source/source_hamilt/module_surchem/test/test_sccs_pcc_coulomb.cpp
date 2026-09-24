#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs_pcc_coulomb.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <vector>

namespace
{

TEST(SccsPccCoulomb, ChargedUniformDielectricScreensPccPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double cube_length = 10.0;
    const double volume = cube_length * cube_length * cube_length;
    basis.initgrids(cube_length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions(basis.nrxx);
    ModuleSccs::PccGeometry geometry
        = ModuleSccs::pcc_geometry(lattice, cube_length, 1.0e-10);
    geometry.origin = ModuleBase::Vector3<double>();
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::PccCoulombOperator coulomb(basis,
                                                 ModuleBase::TWO_PI / cube_length,
                                                 positions,
                                                 volume_element,
                                                 geometry,
                                                 charge_reduction);

    const double solute_charge_value = 1.0;
    const std::vector<double> solute_charge(basis.nrxx, solute_charge_value / volume);
    const std::vector<double> epsilon(basis.nrxx, 5.0);
    const std::vector<ModuleBase::Vector3<double>> grad_log_epsilon(basis.nrxx);
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 100;
    solver.mixing = 0.7;
    solver.tolerance_rms = 1.0e-14;
    solver.tolerance_max = 1.0e-14;

    const ModuleSccs::PolarizationResult result
        = ModuleSccs::solve_polarization(solute_charge,
                                         epsilon,
                                         grad_log_epsilon,
                                         std::vector<double>(),
                                         solver,
                                         coulomb);

    ASSERT_EQ(result.status, ModuleSccs::PolarizationStatus::Converged);
    const double screened_charge = solute_charge_value / 5.0;
    const double expected_potential
        = geometry.parameters.madelung * screened_charge / cube_length;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(result.polarization_charge[ir], -0.8 / volume, 1.0e-14);
        EXPECT_NEAR(result.field.potential[ir], expected_potential, 2.0e-12);
        EXPECT_NEAR(result.field.gradient[ir].x, 0.0, 1.0e-12);
        EXPECT_NEAR(result.field.gradient[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.field.gradient[ir].z, 0.0, 1.0e-12);
    }
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &POOL_WORLD);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Comm_free(&POOL_WORLD);
    MPI_Finalize();
#endif
    return result;
}

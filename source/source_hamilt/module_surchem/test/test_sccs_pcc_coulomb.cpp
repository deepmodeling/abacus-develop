#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../pcc/sccs_pcc_coulomb.h"
#include "../sccs/sccs_charge.h"
#include "../sccs/sccs_periodic.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace
{

// The production sqrt-CG keeps the PCC monopole gauge: no zero-mean shift, and
// the continuum polarization charge of a uniform dielectric is -(1 - 1/eps) q.
TEST(SccsPccCoulomb, SqrtCgKeepsChargedUniformDielectricPccGauge)
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
    const double tpiba = ModuleBase::TWO_PI / cube_length;
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
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    const ModuleSccs::PccCoulombOperator coulomb(basis,
                                                 tpiba,
                                                 positions,
                                                 volume_element,
                                                 geometry,
                                                 charge_reduction);

    // A zero electron density puts the whole cell in bulk solvent.
    const std::vector<double> cavity_density(basis.nrxx, 0.0);
    const std::vector<double> solute_charge(basis.nrxx, 1.0 / volume);
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1.0e-4;
    cavity.density_max = 5.0e-3;
    cavity.epsilon_bulk = 5.0;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.max_iterations = 10;
    solver.tolerance_rms = 1.0e-14;
    solver.tolerance_max = 1.0e-14;
    const std::vector<double> cold_start;
    const ModuleSccs::PeriodicSccsResult result
        = ModuleSccs::solve_chain_sccs_response(cavity_density,
                                                solute_charge,
                                                cavity,
                                                solver,
                                                cold_start,
                                                basis,
                                                tpiba,
                                                coulomb,
                                                polarization_reduction);

    ASSERT_EQ(result.polarization.status, ModuleSccs::PolarizationStatus::Converged);
    EXPECT_EQ(result.polarization.iterations, 1);
    EXPECT_NEAR(result.far_field_polarization_charge,
                -(1.0 - 1.0 / cavity.epsilon_bulk),
                1.0e-12);
    const double expected_potential
        = geometry.parameters.madelung * (1.0 / cavity.epsilon_bulk) / cube_length;
    ASSERT_GT(std::abs(expected_potential), 1.0e-3);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_NEAR(result.polarization.field.potential[ir], expected_potential, 2.0e-12);
        EXPECT_DOUBLE_EQ(result.restart_potential[ir], result.polarization.field.potential[ir]);
        EXPECT_NEAR(result.polarization.polarization_charge[ir], -0.8 / volume, 1.0e-14);
        EXPECT_NEAR(result.polarization.field.gradient[ir].x, 0.0, 1.0e-12);
        EXPECT_NEAR(result.polarization.field.gradient[ir].y, 0.0, 1.0e-12);
        EXPECT_NEAR(result.polarization.field.gradient[ir].z, 0.0, 1.0e-12);
    }
}

TEST(SccsPccCoulomb, ScalarPotentialMatchesFullFieldAndAvoidsGradientTransforms)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const double length = 10.0;
    const double tpiba = ModuleBase::TWO_PI / length;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    const double dv = length * length * length / basis.nxyz;
    std::vector<ModuleBase::Vector3<double>> positions(basis.nrxx);
    for (int index = 0; index < basis.nrxx; ++index)
    {
        const int ix = index / (basis.ny * basis.nplane);
        const int iy = (index / basis.nplane) % basis.ny;
        const int iz = index % basis.nplane + basis.startz_current;
        const double x = length * ix / basis.nx;
        const double y = length * iy / basis.ny;
        const double z = length * iz / basis.nz;
        positions[index] = ModuleBase::Vector3<double>(x, y, z);
    }
    const ModuleSccs::PccGeometry geometry = ModuleSccs::pcc_geometry(lattice, length, 1.0e-10);
    const ModuleSccs::SerialChargeReduction reduction;
    const ModuleSccs::PccCoulombOperator coulomb(basis, tpiba, positions, dv, geometry, reduction);
    const std::vector<double> charge(basis.nrxx, 0.001);
    ModuleSccs::ElectrostaticField field;
    coulomb.apply(charge, field);
    std::vector<double> potential;
    coulomb.apply_potential(charge, potential);
    EXPECT_EQ(potential, field.potential);
    EXPECT_EQ(coulomb.transform_counts().forward_calls, 2);
    EXPECT_EQ(coulomb.transform_counts().inverse_calls, 5);
    const std::vector<double> zero(basis.nrxx, 0.0);
    coulomb.apply_potential(zero, potential);
    for (double value : potential)
    {
        EXPECT_DOUBLE_EQ(value, 0.0);
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

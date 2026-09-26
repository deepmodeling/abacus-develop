#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_pw_charge.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace
{

TEST(SccsPwCharge, ReconstructsNonzeroModesAndRestoresPhysicalZeroMode)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double cell_length = 10.0;
    const double cell_volume = cell_length * cell_length * cell_length;
    const double tpiba = ModuleBase::TWO_PI / cell_length;
    basis.initgrids(cell_length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double ionic_charge = 3.0;
    std::vector<std::complex<double>> ionic_g(basis.npw);
    std::vector<std::complex<double>> potential_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (ig == basis.ig_gge0 || basis.gg[ig] == 0.0)
        {
            ionic_g[ig] = ionic_charge / cell_volume;
        }
        else if (std::abs(std::abs(basis.gdirect[ig].x) - 1.0) < 1.0e-12
                 && std::abs(basis.gdirect[ig].y) < 1.0e-12
                 && std::abs(basis.gdirect[ig].z) < 1.0e-12)
        {
            ionic_g[ig] = 2.0e-3;
            const double kernel
                = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba * tpiba * basis.gg[ig]);
            potential_g[ig] = -kernel * ionic_g[ig];
        }
    }
    std::vector<double> expected_density(basis.nrxx);
    std::vector<double> local_potential(basis.nrxx);
    basis.recip2real(ionic_g.data(), expected_density.data());
    basis.recip2real(potential_g.data(), local_potential.data());

    const std::vector<double> reconstructed
        = ModuleSccs::ionic_charge_from_local_potential(local_potential,
                                                        ionic_charge,
                                                        cell_volume,
                                                        tpiba,
                                                        basis);
    ASSERT_EQ(reconstructed.size(), expected_density.size());
    double maximum_error = 0.0;
    for (std::size_t index = 0; index < reconstructed.size(); ++index)
    {
        maximum_error
            = std::max(maximum_error, std::abs(reconstructed[index] - expected_density[index]));
    }
    EXPECT_LT(maximum_error, 1.0e-12);
}

TEST(SccsPwCharge, ValidatesCubeAndBuildsCellCenteredGrid)
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

    EXPECT_DOUBLE_EQ(ModuleSccs::validate_cubic_cell(lattice, 10.0, 1.0e-10), 10.0);
    const ModuleBase::Vector3<double> center = ModuleSccs::cell_center(lattice, 10.0);
    EXPECT_DOUBLE_EQ(center.x, 5.0);
    EXPECT_DOUBLE_EQ(center.y, 5.0);
    EXPECT_DOUBLE_EQ(center.z, 5.0);

    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, 10.0);
    ASSERT_EQ(positions.size(), static_cast<std::size_t>(basis.nrxx));
    EXPECT_NEAR(positions[0].x, 5.0 / static_cast<double>(basis.nx), 1.0e-14);
    EXPECT_NEAR(positions[0].y, 5.0 / static_cast<double>(basis.ny), 1.0e-14);
    EXPECT_NEAR(positions[0].z, 5.0 / static_cast<double>(basis.nz), 1.0e-14);

    const ModuleBase::Matrix3 orthorhombic(1.0, 0.0, 0.0,
                                           0.0, 1.1, 0.0,
                                           0.0, 0.0, 1.0);
    EXPECT_THROW(ModuleSccs::validate_cubic_cell(orthorhombic, 10.0, 1.0e-10),
                 std::invalid_argument);
}

TEST(SccsPwCharge, KeepsYCoordinatesIndependentOfDistributedZSlab)
{
    ModulePW::PW_Basis basis("cpu", "double");
    basis.nx = 3;
    basis.ny = 4;
    basis.nz = 7;
    basis.nplane = 2;
    basis.startz_current = 3;
    basis.nrxx = basis.nx * basis.ny * basis.nplane;
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 1.0,
                                      0.0, 5.0, 0.0,
                                      1.0, 0.0, 3.0);
    const double lattice_scale = 2.0;
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, lattice_scale);

    const int ix = 1;
    const int iy = 2;
    const int iz_local = 1;
    const int index = (ix * basis.ny + iy) * basis.nplane + iz_local;
    const double fractional_x = 1.5 / 3.0;
    const double fractional_y = 2.5 / 4.0;
    const double fractional_z = 4.5 / 7.0;
    EXPECT_NEAR(positions[index].x,
                lattice_scale * (2.0 * fractional_x + fractional_z),
                1.0e-14);
    EXPECT_NEAR(positions[index].y,
                lattice_scale * 5.0 * fractional_y,
                1.0e-14);
    EXPECT_NEAR(positions[index].z,
                lattice_scale * (fractional_x + 3.0 * fractional_z),
                1.0e-14);
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

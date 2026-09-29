#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../common/pw_grid.h"

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

TEST(SccsPwCharge, BuildsCellCenterAndIntegerNodeGrid)
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

    const ModuleBase::Vector3<double> center = ModuleSccs::cell_center(lattice, 10.0);
    EXPECT_DOUBLE_EQ(center.x, 5.0);
    EXPECT_DOUBLE_EQ(center.y, 5.0);
    EXPECT_DOUBLE_EQ(center.z, 5.0);

    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, 10.0);
    ASSERT_EQ(positions.size(), static_cast<std::size_t>(basis.nrxx));
    EXPECT_NEAR(positions[0].x, 0.0, 1.0e-14);
    EXPECT_NEAR(positions[0].y, 0.0, 1.0e-14);
    EXPECT_NEAR(positions[0].z, 0.0, 1.0e-14);
}

TEST(SccsPwCharge, GridCoordinatesMatchInverseFourierPhase)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, length);
    // Construct cos(G.r) + 0.4 sin(G.r) directly in reciprocal space.
    // Checking all axes detects a half-grid offset independently of the
    // coordinate implementation used by PCC moments and correction fields.
    std::vector<std::complex<double>> coefficients(basis.npw, 0.0);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const ModuleBase::Vector3<double>& g = basis.gdirect[ig];
        // Direct reciprocal coordinates are integer mode indices.
        if (g.x == 1.0 && g.y == 1.0 && g.z == 1.0)
        {
            coefficients[ig] = std::complex<double>(0.5, -0.2);
        }
        else if (g.x == -1.0 && g.y == -1.0 && g.z == -1.0)
        {
            coefficients[ig] = std::complex<double>(0.5, 0.2);
        }
    }
    std::vector<double> values(basis.nrxx);
    basis.recip2real(coefficients.data(), values.data());
    double maximum_error = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double phase = ModuleBase::TWO_PI
                             * (positions[ir].x + positions[ir].y + positions[ir].z) / length;
        const double expected = std::cos(phase) + 0.4 * std::sin(phase);
        const double difference = values[ir] - expected;
        const double error = std::abs(difference);
        maximum_error = std::max(maximum_error, error);
    }
    EXPECT_LT(maximum_error, 1.0e-12);
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
    const double fractional_x = 1.0 / 3.0;
    const double fractional_y = 2.0 / 4.0;
    const double fractional_z = 4.0 / 7.0;
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

#include "source_cell/unitcell.h"

// Minimal mocks to satisfy the linker without pulling in heavy dependencies.
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}

#include "source_base/parallel_reduce.h"
namespace Parallel_Reduce
{
// Serial stub: single-rank Allreduce is a no-op.
template <>
void reduce_all<double>(double* data, int n)
{
}
} // namespace Parallel_Reduce

#include "../dftu_nao_fs_r.h"

#include "source_base/matrix.h"
#include "gtest/gtest.h"

#include <vector>

/***********************************************************************
 * Unit tests for dftu_nao_fs_reduce.cpp.
 *
 * The file provides two post-processing helpers used at the end of the
 * real-space force/stress path:
 *
 *   reduce_force_impl(force, nspin)
 *     - MPI Allreduce over ranks (stubbed no-op here)
 *     - multiply by 2 for nspin != 4 (spin-degeneracy factor)
 *
 *   reduce_stress_impl(ucell, stress_tmp, stress)
 *     - MPI Allreduce over ranks (stubbed no-op here)
 *     - scale by lat0/omega and inflate the 6-component Voigt form
 *       (xx, xy, xz, yy, yz, zz) into the full 3x3 symmetric tensor
 ***********************************************************************/

namespace DFTU_LCAO
{
namespace
{

class DftuFsReduceTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ucell.lat0 = 2.0;
        ucell.omega = 10.0;
    }

    UnitCell ucell;
};

// nspin = 1: spin-degenerate, force must be doubled.
TEST_F(DftuFsReduceTest, ReduceForceNspin1Doubles)
{
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0; force(0, 1) = -2.0; force(0, 2) = 3.5;
    force(1, 0) = 0.0; force(1, 1) = 4.0;  force(1, 2) = -1.0;

    reduce_force_impl(force, 1);

    EXPECT_DOUBLE_EQ(force(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(force(0, 1), -4.0);
    EXPECT_DOUBLE_EQ(force(0, 2), 7.0);
    EXPECT_DOUBLE_EQ(force(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(force(1, 1), 8.0);
    EXPECT_DOUBLE_EQ(force(1, 2), -2.0);
}

// nspin = 2: also doubled (two explicit spin channels summed).
TEST_F(DftuFsReduceTest, ReduceForceNspin2Doubles)
{
    ModuleBase::matrix force(1, 3);
    force(0, 0) = 1.5;
    force(0, 1) = -0.5;
    force(0, 2) = 2.0;

    reduce_force_impl(force, 2);

    EXPECT_DOUBLE_EQ(force(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(force(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(force(0, 2), 4.0);
}

// nspin = 4: non-collinear, no extra scaling.
TEST_F(DftuFsReduceTest, ReduceForceNspin4Unchanged)
{
    ModuleBase::matrix force(1, 3);
    force(0, 0) = 1.5;
    force(0, 1) = -0.5;
    force(0, 2) = 2.0;

    reduce_force_impl(force, 4);

    EXPECT_DOUBLE_EQ(force(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(force(0, 1), -0.5);
    EXPECT_DOUBLE_EQ(force(0, 2), 2.0);
}

// Stress: Voigt (xx, xy, xz, yy, yz, zz) -> 3x3 tensor, scaled by lat0/omega.
TEST_F(DftuFsReduceTest, ReduceStressVoigtToTensor)
{
    std::vector<double> stress_tmp = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    ModuleBase::matrix stress(3, 3);

    reduce_stress_impl(&ucell, stress_tmp, stress);

    const double w = ucell.lat0 / ucell.omega; // 0.2
    EXPECT_DOUBLE_EQ(stress(0, 0), 1.0 * w);
    EXPECT_DOUBLE_EQ(stress(0, 1), 2.0 * w);
    EXPECT_DOUBLE_EQ(stress(0, 2), 3.0 * w);
    EXPECT_DOUBLE_EQ(stress(1, 0), 2.0 * w);
    EXPECT_DOUBLE_EQ(stress(1, 1), 4.0 * w);
    EXPECT_DOUBLE_EQ(stress(1, 2), 5.0 * w);
    EXPECT_DOUBLE_EQ(stress(2, 0), 3.0 * w);
    EXPECT_DOUBLE_EQ(stress(2, 1), 5.0 * w);
    EXPECT_DOUBLE_EQ(stress(2, 2), 6.0 * w);
}

// Stress: result must be symmetric regardless of input values.
TEST_F(DftuFsReduceTest, ReduceStressSymmetric)
{
    std::vector<double> stress_tmp = {-1.0, 0.5, 2.5, -3.0, 1.5, 0.0};
    ModuleBase::matrix stress(3, 3);

    reduce_stress_impl(&ucell, stress_tmp, stress);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            EXPECT_DOUBLE_EQ(stress(i, j), stress(j, i));
        }
    }
}

// Stress: unit cell scaling weight is applied to every component.
TEST_F(DftuFsReduceTest, ReduceStressScalingWeight)
{
    ucell.lat0 = 1.0;
    ucell.omega = 4.0;
    std::vector<double> stress_tmp = {8.0, 0.0, 0.0, 8.0, 0.0, 8.0};
    ModuleBase::matrix stress(3, 3);

    reduce_stress_impl(&ucell, stress_tmp, stress);

    // weight = 1.0 / 4.0 = 0.25; hydrostatic input stays hydrostatic
    EXPECT_DOUBLE_EQ(stress(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(stress(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(stress(2, 2), 2.0);
    EXPECT_DOUBLE_EQ(stress(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(stress(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(stress(1, 2), 0.0);
}

} // namespace
} // namespace DFTU_LCAO

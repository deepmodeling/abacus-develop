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
template <>
void reduce_all<double>(double* data, int n)
{
}
} // namespace Parallel_Reduce

#include "../dftu_nao_fs_r.h"

#include "source_base/matrix.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>

namespace DFTU_LCAO
{
namespace
{

class DftuFsRTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ucell.lat0 = 2.0;
        ucell.omega = 10.0;
    }

    UnitCell ucell;
};

TEST_F(DftuFsRTest, ReduceForceNspin1)
{
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0; force(0, 1) = 2.0; force(0, 2) = 3.0;
    force(1, 0) = 4.0; force(1, 1) = 5.0; force(1, 2) = 6.0;

    reduce_force_impl(force, 1);

    // nspin != 4: force *= 2
    EXPECT_DOUBLE_EQ(force(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(force(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(force(1, 2), 12.0);
}

TEST_F(DftuFsRTest, ReduceForceNspin4)
{
    ModuleBase::matrix force(2, 3);
    force(0, 0) = 1.0; force(0, 1) = 2.0; force(0, 2) = 3.0;

    reduce_force_impl(force, 4);

    // nspin == 4: no scaling
    EXPECT_DOUBLE_EQ(force(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(force(0, 1), 2.0);
}

TEST_F(DftuFsRTest, ReduceStressBasic)
{
    std::vector<double> stress_tmp = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    ModuleBase::matrix stress(3, 3);

    reduce_stress_impl(&ucell, stress_tmp, stress);

    // weight = lat0 / omega = 2.0 / 10.0 = 0.2
    const double w = 0.2;
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

} // namespace
} // namespace DFTU_LCAO

#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include "source_base/matrix3.h"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(SternheimerABACUSFDAdapter, BuildsOrthogonalGridFromABACUSLattice)
{
    const ModuleBase::Matrix3 latvec(3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0);
    const auto grid_data = ModuleRI::make_sternheimer_fd_grid_from_lattice(6, 8, 10, 6 * 8 * 10, 2.0, latvec);

    EXPECT_EQ(grid_data.grid.nx, 6);
    EXPECT_EQ(grid_data.grid.ny, 8);
    EXPECT_EQ(grid_data.grid.nz, 10);
    EXPECT_TRUE(grid_data.grid.periodic);
    EXPECT_DOUBLE_EQ(grid_data.grid.hx, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.grid.hy, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.grid.hz, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.volume_element, 1.0);
}

TEST(SternheimerABACUSFDAdapter, RejectsDistributedRealSpaceGridForDensePrototype)
{
    const ModuleBase::Matrix3 latvec(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    EXPECT_THROW(ModuleRI::make_sternheimer_fd_grid_from_lattice(4, 4, 4, 32, 8.0, latvec), std::invalid_argument);
}

TEST(SternheimerABACUSFDAdapter, RejectsNonOrthogonalLatticeForCurrentStencil)
{
    const ModuleBase::Matrix3 latvec(1.0, 0.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    EXPECT_THROW(ModuleRI::make_sternheimer_fd_grid_from_lattice(4, 4, 4, 64, 8.0, latvec), std::invalid_argument);
}

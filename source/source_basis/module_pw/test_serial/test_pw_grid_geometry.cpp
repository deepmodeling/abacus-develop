#include "source_basis/module_pw/pw_grid_geometry.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_base/matrix3.h"

#include <gtest/gtest.h>

TEST(PwGridGeometry, UsesGlobalSlabOffsetAndRotatedLattice)
{
    ModulePW::PW_Basis basis("cpu", "double");
    basis.nx = 2;
    basis.ny = 2;
    basis.nz = 4;
    basis.nplane = 2;
    basis.startz_current = 2;
    basis.nrxx = 8;
    const ModuleBase::Matrix3 lattice(0.0, 1.0, 0.0,
                                     -2.0, 0.0, 0.0,
                                     0.0, 0.0, 3.0);
    std::vector<ModuleBase::Vector3<double>> positions;
    std::string error;
    ASSERT_TRUE(ModulePW::grid_positions(basis, lattice, 10.0, positions, error));
    ASSERT_EQ(positions.size(), 8u);
    EXPECT_DOUBLE_EQ(positions[0].x, 0.0);
    EXPECT_DOUBLE_EQ(positions[0].y, 0.0);
    EXPECT_DOUBLE_EQ(positions[0].z, 15.0);
    EXPECT_DOUBLE_EQ(positions[7].x, -10.0);
    EXPECT_DOUBLE_EQ(positions[7].y, 5.0);
    EXPECT_DOUBLE_EQ(positions[7].z, 22.5);
}

TEST(PwGridGeometry, PermitsEmptySlabsAndRejectsInvalidOffsets)
{
    ModulePW::PW_Basis basis("cpu", "double");
    basis.nx = 2;
    basis.ny = 2;
    basis.nz = 4;
    basis.nplane = 0;
    basis.startz_current = 4;
    basis.nrxx = 0;
    const ModuleBase::Matrix3 lattice;
    std::vector<ModuleBase::Vector3<double>> positions;
    std::string error;
    ASSERT_TRUE(ModulePW::grid_positions(basis, lattice, 10.0, positions, error));
    EXPECT_TRUE(positions.empty());
    basis.nplane = 1;
    basis.nrxx = 4;
    EXPECT_FALSE(ModulePW::grid_positions(basis, lattice, 10.0, positions, error));
    EXPECT_TRUE(positions.empty());
    EXPECT_FALSE(error.empty());
}

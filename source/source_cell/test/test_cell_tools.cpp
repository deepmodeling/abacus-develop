#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>
#include <string>
#include <vector>

#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "source_cell/unitcell.h"
#include "source_cell/cell_tools.h"
#include "prepare_unitcell.h"

// The test-only cell_info object library does not contain magnetism.cpp,
// so the Magnetism constructor/destructor must be provided locally.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of cell_tools.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - IfCellCanChange
 *     - if_cell_can_change(): truth table over the three lattice-axis flags
 *   - GetAtomCounts
 *     - get_atomCounts(): number of atoms per type as a vector
 *   - GetLnchiCounts
 *     - get_lnchiCounts(): number of chi functions per L per type
 *   - SelectiveDynamics
 *     - if_atoms_can_move(): true if any atom has a movable coordinate
 */

class CellToolsTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
};

TEST_F(CellToolsTest, IfCellCanChange)
{
    // Mirror the fixed_axes -> lat_axis_free mapping produced by
    // UnitCell::setup_from_input: the cell can change whenever at least one
    // lattice axis is free; only "abc" (all axes fixed) returns false.
    std::vector<std::vector<int>> axes_free = {
        {1, 1, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0},
        {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
    for (int i = 0; i < 7; ++i)
    {
        EXPECT_TRUE(unitcell::if_cell_can_change(axes_free[i]));
    }
    EXPECT_FALSE(unitcell::if_cell_can_change(axes_free[7]));
}

TEST_F(CellToolsTest, GetAtomCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    ucell->set_iat2itia();
    std::vector<int> atomCounts = unitcell::get_atomCounts(ucell->atoms, ucell->ntype);
    EXPECT_EQ(atomCounts[0], 1);
    EXPECT_EQ(atomCounts[1], 2);
}

TEST_F(CellToolsTest, GetLnchiCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    ucell->set_iat2itia();
    std::vector<std::vector<int>> lnchiCounts = unitcell::get_lnchiCounts(ucell->atoms, ucell->ntype);
    EXPECT_EQ(lnchiCounts[0][0], 1);
    EXPECT_EQ(lnchiCounts[0][1], 1);
    EXPECT_EQ(lnchiCounts[0][2], 1);
    EXPECT_EQ(lnchiCounts[1][0], 1);
    EXPECT_EQ(lnchiCounts[1][1], 1);
    EXPECT_EQ(lnchiCounts[1][2], 1);
}

TEST_F(CellToolsTest, SelectiveDynamics)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-SD"];
    ucell = utp.SetUcellInfo();
    EXPECT_TRUE(unitcell::if_atoms_can_move(ucell->atoms, ucell->ntype));
}

TEST(CellTools, ExtractsPositionsMassesAndValenceChargesInAtomOrder)
{
    Atom atoms[2];
    atoms[0].na = 2;
    atoms[0].mass = 12.0;
    atoms[0].ncpp.zv = 4.0;
    atoms[0].tau = {ModuleBase::Vector3<double>(0.1, 0.2, 0.3),
                    ModuleBase::Vector3<double>(0.4, 0.5, 0.6)};
    atoms[1].na = 1;
    atoms[1].mass = 1.0;
    atoms[1].ncpp.zv = 1.0;
    atoms[1].tau = {ModuleBase::Vector3<double>(0.7, 0.8, 0.9)};

    const std::vector<unitcell::AtomData> data = unitcell::get_atom_data(atoms, 2, 10.0);
    ASSERT_EQ(data.size(), 3u);
    EXPECT_DOUBLE_EQ(data[0].position.x, 1.0);
    EXPECT_DOUBLE_EQ(data[0].position.y, 2.0);
    EXPECT_DOUBLE_EQ(data[0].position.z, 3.0);
    EXPECT_DOUBLE_EQ(data[1].position.x, 4.0);
    EXPECT_DOUBLE_EQ(data[1].position.y, 5.0);
    EXPECT_DOUBLE_EQ(data[1].position.z, 6.0);
    EXPECT_DOUBLE_EQ(data[2].position.x, 7.0);
    EXPECT_DOUBLE_EQ(data[2].position.y, 8.0);
    EXPECT_DOUBLE_EQ(data[2].position.z, 9.0);
    EXPECT_DOUBLE_EQ(data[0].mass, 12.0);
    EXPECT_DOUBLE_EQ(data[1].mass, 12.0);
    EXPECT_DOUBLE_EQ(data[2].mass, 1.0);
    EXPECT_DOUBLE_EQ(data[0].valence_charge, 4.0);
    EXPECT_DOUBLE_EQ(data[1].valence_charge, 4.0);
    EXPECT_DOUBLE_EQ(data[2].valence_charge, 1.0);
    EXPECT_DOUBLE_EQ(atoms[0].tau[0].x, 0.1);
}

TEST(CellTools, AcceptsNoAtomTypes)
{
    const std::vector<unitcell::AtomData> data = unitcell::get_atom_data(nullptr, 0, 10.0);
    EXPECT_TRUE(data.empty());
}

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

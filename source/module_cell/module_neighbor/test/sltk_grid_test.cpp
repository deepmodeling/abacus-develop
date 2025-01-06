#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private

#define private public
#include "../sltk_grid.h"
#include "prepare_unitcell.h"
#undef private
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

/************************************************
 *  unit test of sltk_grid
 ***********************************************/

/**
 * - Tested Functions:
 *   - Init: Grid::init()
 *     - setMemberVariables: really set member variables
 *       (like dx, dy, dz and d_minX, d_minY, d_minZ) by
 *       reading from getters of Atom_input, and construct the
 *       member Cell as a 3D array of CellSet
 */

void SetGlobalV()
{
    PARAM.input.test_grid = 0;
}

class SltkGridTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::ofstream ofs;
    std::ifstream ifs;
    bool pbc = true;
    double radius = ((8 + 5.01) * 2.0 + 0.01) / 10.2;
    int test_atom_in = 0;
    std::string output;
    void SetUp()
    {
        SetGlobalV();
        ucell = utp.SetUcellInfo();
    }
    void TearDown()
    {
        delete ucell;
    }
};

using SltkGridDeathTest = SltkGridTest;

TEST_F(SltkGridTest, Init)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, radius, pbc);
    EXPECT_EQ(LatGrid.getGlayerX(), 6);
    EXPECT_EQ(LatGrid.getGlayerY(), 6);
    EXPECT_EQ(LatGrid.getGlayerZ(), 6);
    EXPECT_EQ(LatGrid.getGlayerX_minus(), 5);
    EXPECT_EQ(LatGrid.getGlayerY_minus(), 5);
    EXPECT_EQ(LatGrid.getGlayerZ_minus(), 5);
    ofs.close();
    remove("test.out");
}

TEST_F(SltkGridTest, InitSmall)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    radius = 0.5;
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, radius, pbc);
    LatGrid.setMemberVariables(ofs,  *ucell);
    EXPECT_EQ(LatGrid.pbc, true);
    EXPECT_TRUE(LatGrid.pbc);
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, radius * radius);
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, 0.5 * 0.5);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, radius);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, 0.5);
    /*
    // minimal value of x, y, z
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_x, 1);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_y, 1);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_z, 1);
    // number of cells in x, y, z
    EXPECT_EQ(LatGrid.cell_nx, 3);
    EXPECT_EQ(LatGrid.cell_ny, 3);
    EXPECT_EQ(LatGrid.cell_nz, 3);
    */
    ofs.close();
    remove("test.out");
}

/*
// This test cannot pass because setAtomLinkArray() is unsuccessful
// if expand_flag is false
TEST_F(SltkGridTest, InitNoExpand)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    PARAM.input.test_grid = 1;
    double radius = 1e-1000;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    ofs.close();
}
*/

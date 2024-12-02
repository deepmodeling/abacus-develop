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
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    EXPECT_TRUE(LatGrid.init_cell_flag);
    EXPECT_EQ(LatGrid.getCellX(), 11);
    EXPECT_EQ(LatGrid.getCellY(), 11);
    EXPECT_EQ(LatGrid.getCellZ(), 11);
    EXPECT_EQ(LatGrid.getTrueCellX(), 5);
    EXPECT_EQ(LatGrid.getTrueCellY(), 5);
    EXPECT_EQ(LatGrid.getTrueCellZ(), 5);
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
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(PARAM.input.test_grid);
    LatGrid.setMemberVariables(ofs, Atom_inp);
    EXPECT_EQ(LatGrid.pbc, Atom_inp.getBoundary());
    EXPECT_TRUE(LatGrid.pbc);
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, Atom_inp.getRadius() * Atom_inp.getRadius());
    EXPECT_DOUBLE_EQ(LatGrid.sradius2, 0.5 * 0.5);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, Atom_inp.getRadius());
    EXPECT_DOUBLE_EQ(LatGrid.sradius, 0.5);
    
    // minimal value of x, y, z
    EXPECT_DOUBLE_EQ(LatGrid.d_minX, Atom_inp.minX());
    EXPECT_DOUBLE_EQ(LatGrid.d_minY, Atom_inp.minY());
    EXPECT_DOUBLE_EQ(LatGrid.d_minZ, Atom_inp.minZ());
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_x, 2);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_y, 2);
    EXPECT_DOUBLE_EQ(LatGrid.true_cell_z, 2);
    // number of cells in x, y, z
    EXPECT_EQ(LatGrid.cell_nx, Atom_inp.getCell_nX());
    EXPECT_EQ(LatGrid.cell_ny, Atom_inp.getCell_nY());
    EXPECT_EQ(LatGrid.cell_nz, Atom_inp.getCell_nZ());
    EXPECT_EQ(LatGrid.cell_nx, 4);
    EXPECT_EQ(LatGrid.cell_ny, 4);
    EXPECT_EQ(LatGrid.cell_nz, 4);
    // init cell flag
    EXPECT_TRUE(LatGrid.init_cell_flag);

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

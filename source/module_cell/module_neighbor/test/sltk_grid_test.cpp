#include "gmock/gmock.h"
#include "gtest/gtest.h"

#define private public
#include "../sltk_grid.h"
#include "prepare_unitcell.h"

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
    GlobalV::test_grid = 0;
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
    GlobalV::test_grid = 1;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(GlobalV::test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    EXPECT_TRUE(LatGrid.init_cell_flag);
    EXPECT_EQ(LatGrid.getCellX(), 11);
    EXPECT_EQ(LatGrid.getCellY(), 11);
    EXPECT_EQ(LatGrid.getCellZ(), 11);
    EXPECT_EQ(LatGrid.getD_minX(), -5);
    EXPECT_EQ(LatGrid.getD_minY(), -5);
    EXPECT_EQ(LatGrid.getD_minZ(), -5);
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
    GlobalV::test_grid = 1;
    double radius = 1e-1000;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(GlobalV::test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    ofs.close();
}
*/

#undef private

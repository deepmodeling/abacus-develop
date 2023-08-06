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

TEST_F(SltkGridTest, InitSmall)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    GlobalV::test_grid = 1;
    radius = 0.5;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(GlobalV::test_grid);
    LatGrid.setMemberVariables(ofs, Atom_inp);
    EXPECT_EQ(LatGrid.natom, Atom_inp.getAmount());
    EXPECT_EQ(LatGrid.natom, 128);
    EXPECT_EQ(LatGrid.pbc, Atom_inp.getBoundary());
    EXPECT_TRUE(LatGrid.pbc);
    EXPECT_DOUBLE_EQ(LatGrid.sradius, Atom_inp.getRadius());
    EXPECT_DOUBLE_EQ(LatGrid.sradius, 0.5);
    for (int i = 0; i < 3; i++)
    {
        EXPECT_DOUBLE_EQ(LatGrid.vec1[i], Atom_inp.vec1[i]);
        EXPECT_DOUBLE_EQ(LatGrid.vec2[i], Atom_inp.vec2[i]);
        EXPECT_DOUBLE_EQ(LatGrid.vec3[i], Atom_inp.vec3[i]);
    }
    EXPECT_EQ(LatGrid.vec1[0], -0.5);
    EXPECT_EQ(LatGrid.vec1[1], 0.0);
    EXPECT_EQ(LatGrid.vec1[2], 0.5);
    EXPECT_EQ(LatGrid.vec2[0], 0.0);
    EXPECT_EQ(LatGrid.vec2[1], 0.5);
    EXPECT_EQ(LatGrid.vec2[2], 0.5);
    EXPECT_EQ(LatGrid.vec3[0], -0.5);
    EXPECT_EQ(LatGrid.vec3[1], 0.5);
    EXPECT_EQ(LatGrid.vec3[2], 0.0);
    // lattice grid length
    EXPECT_DOUBLE_EQ(LatGrid.grid_length[0], Atom_inp.Clength0());
    EXPECT_DOUBLE_EQ(LatGrid.grid_length[1], Atom_inp.Clength1());
    EXPECT_DOUBLE_EQ(LatGrid.grid_length[2], Atom_inp.Clength2());
    EXPECT_NEAR(LatGrid.grid_length[0], 2.8284271247461903, 1e-15);
    EXPECT_NEAR(LatGrid.grid_length[1], 2.8284271247461903, 1e-15);
    EXPECT_NEAR(LatGrid.grid_length[2], 2.8284271247461903, 1e-15);
    // minimal value of x, y, z
    EXPECT_DOUBLE_EQ(LatGrid.d_minX, Atom_inp.minX());
    EXPECT_DOUBLE_EQ(LatGrid.d_minY, Atom_inp.minY());
    EXPECT_DOUBLE_EQ(LatGrid.d_minZ, Atom_inp.minZ());
    EXPECT_DOUBLE_EQ(LatGrid.d_minX, -2);
    EXPECT_DOUBLE_EQ(LatGrid.d_minY, -2);
    EXPECT_DOUBLE_EQ(LatGrid.d_minZ, -2);
    // length of cells in x, y, z
    EXPECT_DOUBLE_EQ(LatGrid.cell_x_length, Atom_inp.getCellXLength());
    EXPECT_DOUBLE_EQ(LatGrid.cell_y_length, Atom_inp.getCellYLength());
    EXPECT_DOUBLE_EQ(LatGrid.cell_z_length, Atom_inp.getCellZLength());
    EXPECT_DOUBLE_EQ(LatGrid.cell_x_length, 1.0);
    EXPECT_DOUBLE_EQ(LatGrid.cell_y_length, 1.0);
    EXPECT_DOUBLE_EQ(LatGrid.cell_z_length, 1.0);
    // number of cells in x, y, z
    EXPECT_EQ(LatGrid.dx, Atom_inp.getCellX());
    EXPECT_EQ(LatGrid.dy, Atom_inp.getCellY());
    EXPECT_EQ(LatGrid.dz, Atom_inp.getCellZ());
    EXPECT_EQ(LatGrid.dx, 4);
    EXPECT_EQ(LatGrid.dy, 4);
    EXPECT_EQ(LatGrid.dz, 4);
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
    GlobalV::test_grid = 1;
    double radius = 1e-1000;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid(GlobalV::test_grid);
    LatGrid.init(ofs, *ucell, Atom_inp);
    ofs.close();
}
*/

#undef private

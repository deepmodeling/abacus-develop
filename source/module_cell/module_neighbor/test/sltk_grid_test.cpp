#include "../sltk_grid.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
 *   - Getters
 *     - get the dx, dy, dz, dx_min, dy_min, dz_min
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

TEST_F(SltkGridTest, Foo)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    GlobalV::test_grid = 1;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    Grid LatGrid;
    LatGrid.init(ofs, *ucell, Atom_inp);
    ofs.close();
}

/*
TEST_F(SltkGridTest,Getters)
{
    Grid test(1);
    test.dx=1;
    test.dy=2;
    test.dz=3;
    test.d_minX=4;
    test.d_minY=5;
    test.d_minZ=6;
    double testx=test.getCellX();
    double testy=test.getCellY();
    double testz=test.getCellZ();
    double testminx=test.getD_minX();
    double testminy=test.getD_minY();
    double testminz=test.getD_minZ();
    EXPECT_EQ(testx,1);
    EXPECT_EQ(testy,2);
    EXPECT_EQ(testz,3);
    EXPECT_EQ(testminx,4);
    EXPECT_EQ(testminy,5);
    EXPECT_EQ(testminz,6);
}
*/

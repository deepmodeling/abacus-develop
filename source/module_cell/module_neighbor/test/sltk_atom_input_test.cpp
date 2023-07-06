#include "../sltk_atom_input.h"

#include <string>

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
 * unit test of sltk_atom_input.cpp
 * *********************************************/

/**
 * - Tested Functions:
 *  - read_atom()
 */

void SetGlobalV()
{
    GlobalV::test_grid = 0;
}

class SltkAtomInputTest : public ::testing::Test
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

using SltkAtomInputDeathTest = SltkAtomInputTest;

TEST_F(SltkAtomInputDeathTest, ConstructorWarning1)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    radius = -1;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("search radius < 0,forbidden"));
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputDeathTest, ConstructorWarning2)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].x = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dminX<0.0"));
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputDeathTest, ConstructorWarning3)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].y = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dminY<0.0"));
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputDeathTest, ConstructorWarning4)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    ucell->atoms[0].taud[1].z = -0.25;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dminZ<0.0"));
    ofs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputTest, ConstructorSmall)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 1;
    GlobalV::test_grid = 1;
    radius = 0.5;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ntype = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Amount(atom number) = 2"));
    EXPECT_THAT(str, testing::HasSubstr("Periodic_boundary = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Searching radius(lat0) = 0.5"));
    EXPECT_THAT(str, testing::HasSubstr("CellLength(unit: lat0) = [ 0.707107, 0.707107, 0.707107 ]"));
    EXPECT_THAT(str, testing::HasSubstr("min_tau = [ -0.75, 0, 0 ]"));
    EXPECT_THAT(str, testing::HasSubstr("max_tau = [ 0, 0.75, 0.75 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer+ = [ 2, 2, 2 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer- = [ 2, 2, 2 ]"));
    EXPECT_THAT(str, testing::HasSubstr("CellDim = [ 4, 4, 4 ]"));
    ifs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputTest, Constructor)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    GlobalV::test_grid = 1;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    ofs.close();
    ifs.open("test.out");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ntype = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Amount(atom number) = 2"));
    EXPECT_THAT(str, testing::HasSubstr("Periodic_boundary = 1"));
    EXPECT_THAT(str, testing::HasSubstr("Searching radius(lat0) = 2.55"));
    EXPECT_THAT(str, testing::HasSubstr("CellLength(unit: lat0) = [ 0.707107, 0.707107, 0.707107 ]"));
    EXPECT_THAT(str, testing::HasSubstr("min_tau = [ -0.75, 0, 0 ]"));
    EXPECT_THAT(str, testing::HasSubstr("max_tau = [ 0, 0.75, 0.75 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer+ = [ 6, 6, 6 ]"));
    EXPECT_THAT(str, testing::HasSubstr("glayer- = [ 5, 5, 5 ]"));
    EXPECT_THAT(str, testing::HasSubstr("CellDim = [ 11, 11, 11 ]"));
    ifs.close();
    remove("test.out");
}

TEST_F(SltkAtomInputTest, Getters)
{
    ofs.open("test.out");
    ucell->check_dtau();
    test_atom_in = 2;
    Atom_input Atom_inp(ofs, *ucell, ucell->nat, ucell->ntype, pbc, radius, test_atom_in);
    EXPECT_TRUE(Atom_inp.getExpandFlag());
    EXPECT_EQ(Atom_inp.getAmount(), 2662);
    EXPECT_EQ(Atom_inp.getBoundary(), 1);
    EXPECT_DOUBLE_EQ(Atom_inp.getLatNow(), 10.2);
    EXPECT_NEAR(Atom_inp.getRadius(), 2.55196, 1e-5);
    EXPECT_DOUBLE_EQ(Atom_inp.getCellXLength(), 1.0);
    EXPECT_DOUBLE_EQ(Atom_inp.getCellYLength(), 1.0);
    EXPECT_DOUBLE_EQ(Atom_inp.getCellZLength(), 1.0);
    EXPECT_NEAR(Atom_inp.Clength0(), 7.77817, 1e-5);
    EXPECT_NEAR(Atom_inp.Clength1(), 7.77817, 1e-5);
    EXPECT_NEAR(Atom_inp.Clength2(), 7.77817, 1e-5);
    EXPECT_DOUBLE_EQ(Atom_inp.minX(), -5.0);
    EXPECT_DOUBLE_EQ(Atom_inp.minY(), -5.0);
    EXPECT_DOUBLE_EQ(Atom_inp.minZ(), -5.0);
    EXPECT_EQ(Atom_inp.getCellX(), 11);
    EXPECT_EQ(Atom_inp.getCellY(), 11);
    EXPECT_EQ(Atom_inp.getCellZ(), 11);
    EXPECT_EQ(Atom_inp.getGrid_layerX(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerY(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerZ(), 6);
    EXPECT_EQ(Atom_inp.getGrid_layerX_minus(), 5);
    EXPECT_EQ(Atom_inp.getGrid_layerY_minus(), 5);
    EXPECT_EQ(Atom_inp.getGrid_layerZ_minus(), 5);
    ofs.close();
    remove("test.out");
}

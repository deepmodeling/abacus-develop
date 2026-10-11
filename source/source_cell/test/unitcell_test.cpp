#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_cell/read_stru.h"
#include "memory"
#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "prepare_unitcell.h"
#include <streambuf>
#include <type_traits>
#include <valarray>
#include <vector>


Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor:
 *     - UnitCell() and ~UnitCell()
 *   - Setup:
 *     - setup_from_input(): to set latname, ntype, lmaxmax, init_vel, and lc
 *   - Index
 *     - set_iat2iait(): set index relations in two arrays of Unitcell: iat2it[nat], iat2ia[nat]
 *     - iat2iait(): depends on the above function, but can find both ia & it from iat
 *     - ijat2iaitjajt(): find ia, it, ja, jt from ijat (ijat_max = nat*nat)
 *         which collapses it, ia, jt, ja loop into a single loop
 *   - GetAtomCounts
 *     - get_atom_Counts(): get atomCounts, which is a map from atom type to atom number
 *   - GetOrbitalCounts
 *     - get_orbitalCounts(): get orbitalCounts, which is a map from atom type to orbital number
 *   - CheckDTau
 *     - check_dtau(): move all atomic coordinates into the first unitcell, i.e. in between [0,1)
 *   - CheckTau
 *     - check_tau(): check if any "two atoms are too close"
 */

class UcellTest : public ::testing::Test
{
  protected:
    std::unique_ptr<UnitCell> ucell{new UnitCell};
    std::string output;
};

/// Compile-time guard: UnitCell owns the raw pointer 'atoms' tracked by
/// 'set_atom_flag' and exposes reference aliases into its 'lat' member, so
/// neither copy nor move can be implemented correctly until that ownership
/// model is refactored. Re-enabling either operation must fail this TU.
TEST(UnitCellTypeTraits, NotCopyableOrMovable)
{
    static_assert(!std::is_copy_constructible<UnitCell>::value,
                  "UnitCell must not be copy constructible");
    static_assert(!std::is_copy_assignable<UnitCell>::value,
                  "UnitCell must not be copy assignable");
    static_assert(!std::is_move_constructible<UnitCell>::value,
                  "UnitCell must not be move constructible");
    static_assert(!std::is_move_assignable<UnitCell>::value,
                  "UnitCell must not be move assignable");
}

using UcellDeathTest = UcellTest;

TEST_F(UcellTest, Constructor)
{
    EXPECT_EQ(ucell->Coordinate, "Direct");
    EXPECT_EQ(ucell->latName, "user_defined_lattice");
    EXPECT_DOUBLE_EQ(ucell->lat0, 0.0);
    EXPECT_DOUBLE_EQ(ucell->lat0_angstrom, 0.0);
    EXPECT_EQ(ucell->ntype, 0);
    EXPECT_EQ(ucell->nat, 0);
    EXPECT_EQ(ucell->namax, 0);
    EXPECT_EQ(ucell->nwmax, 0);
    EXPECT_TRUE(ucell->iat2it.empty());
    EXPECT_TRUE(ucell->iat2ia.empty());
    EXPECT_TRUE(ucell->iwt2iat.empty());
    EXPECT_TRUE(ucell->iwt2iw.empty());
    EXPECT_DOUBLE_EQ(ucell->tpiba, 0.0);
    EXPECT_DOUBLE_EQ(ucell->tpiba2, 0.0);
    EXPECT_DOUBLE_EQ(ucell->omega, 0.0);
    EXPECT_FALSE(ucell->set_atom_flag);
}

TEST_F(UcellTest, Setup)
{
    std::string latname_in = "bcc";
    int ntype_in = 1;
    int lmaxmax_in = 2;
    bool init_vel_in = false;
    std::vector<std::string> fixed_axes_in = {"None", "volume", "shape", "a", "b", "c", "ab", "ac", "bc", "abc"};
    for (int i = 0; i < fixed_axes_in.size(); ++i)
    {
        ucell->setup_from_input(latname_in, ntype_in, lmaxmax_in, init_vel_in, fixed_axes_in[i]);
        EXPECT_EQ(ucell->latName, latname_in);
        EXPECT_EQ(ucell->ntype, ntype_in);
        EXPECT_EQ(ucell->lmaxmax, lmaxmax_in);
        EXPECT_EQ(ucell->init_vel, init_vel_in);
        if (fixed_axes_in[i] == "None" || fixed_axes_in[i] == "volume" || fixed_axes_in[i] == "shape")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "a")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "b")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "c")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "ab")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 1);
        }
        else if (fixed_axes_in[i] == "ac")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 1);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "bc")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 1);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
        else if (fixed_axes_in[i] == "abc")
        {
            EXPECT_EQ(ucell->lat_axis_free[0], 0);
            EXPECT_EQ(ucell->lat_axis_free[1], 0);
            EXPECT_EQ(ucell->lat_axis_free[2], 0);
        }
    }
}

TEST_F(UcellTest, Index)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    int iat = 0;
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            EXPECT_EQ(ucell->iat2it[iat], it);
            EXPECT_EQ(ucell->iat2ia[iat], ia);
            // test iat2iait
            int ia_beg, it_beg;
            ucell->iat2iait(iat, &ia_beg, &it_beg);
            EXPECT_EQ(it_beg, it);
            EXPECT_EQ(ia_beg, ia);
            ++iat;
        }
    }
    // test iat2iait: case of (iat >= nat)
    int ia_beg2;
    int it_beg2;
    long long iat2 = ucell->nat + 1;
    EXPECT_FALSE(ucell->iat2iait(iat2, &ia_beg2, &it_beg2));
    // test ijat2iaitjajt
    int ia_test;
    int it_test;
    int ja_test;
    int jt_test;
    long long ijat = 0;
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            for (int jt = 0; jt < utp.natom.size(); ++jt)
            {
                for (int ja = 0; ja < utp.natom[jt]; ++ja)
                {
                    ucell->ijat2iaitjajt(ijat, &ia_test, &it_test, &ja_test, &jt_test);
                    EXPECT_EQ(ia_test, ia);
                    EXPECT_EQ(it_test, it);
                    EXPECT_EQ(ja_test, ja);
                    EXPECT_EQ(jt_test, jt);
                    ++ijat;
                }
            }
        }
    }
}

TEST_F(UcellTest, GetAtomCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, int> atomCounts = ucell->get_atom_Counts();
    EXPECT_EQ(atomCounts[0], 1);
    EXPECT_EQ(atomCounts[1], 2);
}

TEST_F(UcellTest, GetOrbitalCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, int> orbitalCounts = ucell->get_orbital_Counts();
    EXPECT_EQ(orbitalCounts[0], 9);
    EXPECT_EQ(orbitalCounts[1], 9);
}

TEST_F(UcellTest, GetLnchiCounts)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    // test set_iat2itia
    ucell->set_iat2itia();
    std::map<int, std::map<int, int>> LnchiCounts = ucell->get_lnchi_Counts();
    EXPECT_EQ(LnchiCounts[0][0], 1);
    EXPECT_EQ(LnchiCounts[0][1], 1);
    EXPECT_EQ(LnchiCounts[0][2], 1);
    EXPECT_EQ(LnchiCounts[1][0], 1);
    EXPECT_EQ(LnchiCounts[1][1], 1);
    EXPECT_EQ(LnchiCounts[1][2], 1);
}

TEST_F(UcellTest, CheckDTau)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckDTau"];
    ucell = utp.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms,ucell->ntype, ucell->lat0, ucell->latvec);
    for (int it = 0; it < utp.natom.size(); ++it)
    {
        for (int ia = 0; ia < utp.natom[it]; ++ia)
        {
            EXPECT_GE(ucell->atoms[it].taud[ia].x, 0);
            EXPECT_GE(ucell->atoms[it].taud[ia].y, 0);
            EXPECT_GE(ucell->atoms[it].taud[ia].z, 0);
            EXPECT_LT(ucell->atoms[it].taud[ia].x, 1);
            EXPECT_LT(ucell->atoms[it].taud[ia].y, 1);
            EXPECT_LT(ucell->atoms[it].taud[ia].z, 1);
        }
    }
}

TEST_F(UcellTest, CheckTauFalse)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckTau"];
    ucell = utp.SetUcellInfo();
    GlobalV::ofs_warning.open("checktau_warning");
    unitcell::check_tau(ucell->atoms ,ucell->ntype, ucell->lat0);
    GlobalV::ofs_warning.close();
    std::ifstream ifs;
    ifs.open("checktau_warning");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("two atoms are too close!"));
    ifs.close();
    remove("checktau_warning");
}

TEST_F(UcellTest, CheckTauTrue)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-CheckTau"];
    ucell = utp.SetUcellInfo();
    GlobalV::ofs_warning.open("checktau_warning");
    int atom=0;
    //cause the ucell->lat0 is 0.5,if the type of the check_tau has 
    //an int type,it will set to zero,and it will not pass the unittest
    ucell->lat0=0.5;
    ucell->nat=3;
    for (int it=0;it<ucell->ntype;it++)
    {
        for(int ia=0; ia<ucell->atoms[it].na; ++ia)
        {
            
            for (int i=0;i<3;i++)
            {
                ucell->atoms[it].tau[ia][i]=((atom+i)/(ucell->nat*3.0));
            }
            atom+=3;
        }
    }
    EXPECT_EQ(unitcell::check_tau(ucell->atoms ,ucell->ntype, ucell->lat0),true);
    GlobalV::ofs_warning.close();
}

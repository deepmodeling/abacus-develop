#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

/************************************************
 *  unit test of the functions in cal_mw_helper.cpp
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::convert
 *    - convert the data structure for calculation of mw
 *  - SpinConstrain::calculate_MW
 *    - calculate mw from AorbMulP matrix
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, CalculateMW)
{
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    std::map<int, std::map<int, int>> lnchiCounts = {
        {0, {{0, 1}}}
    };
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_lnchiCounts(lnchiCounts);
    sc.set_nspin(4);
    sc.set_npol(2);
    ModuleBase::matrix orbMulP;
    int nlocal = sc.get_nw() / 2;
    orbMulP.create(sc.get_nspin(), nlocal, true);
    EXPECT_EQ(orbMulP.nc, 1);
    EXPECT_EQ(orbMulP.nr, 4);
    orbMulP(0, 0) = 1.0;
    orbMulP(1, 0) = 2.0;
    orbMulP(2, 0) = 3.0;
    orbMulP(3, 0) = 4.0;
    std::vector<std::vector<std::vector<double>>> AorbMulP = sc.convert(orbMulP);
    EXPECT_EQ(AorbMulP.size(), 4);
    EXPECT_EQ(AorbMulP[0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0].size(), 1);
    EXPECT_EQ(AorbMulP[0][0][0], 1.0);
    EXPECT_EQ(AorbMulP[1][0][0], 2.0);
    EXPECT_EQ(AorbMulP[2][0][0], 3.0);
    EXPECT_EQ(AorbMulP[3][0][0], 4.0);
    // calculate_MW
    sc.calculate_MW(AorbMulP);
    testing::internal::CaptureStdout();
    sc.print_Mi(true);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Total Magnetism on atom: 0  (2, 3, 4)"));
}
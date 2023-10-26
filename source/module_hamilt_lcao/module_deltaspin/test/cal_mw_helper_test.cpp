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
 *  unit test of the function of cal_mw
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::cal_mw
 *    - this function calculates the atomic magnetic moment
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, Convert)
{
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.clear_atomCounts();
    sc.clear_orbitalCounts();
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
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
}
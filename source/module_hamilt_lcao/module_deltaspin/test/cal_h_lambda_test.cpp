#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

/************************************************
 *  unit test of the function of cal_h_lambda
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::cal_h_lambda
 *    - this function calculates the h_lambda operator
 */

class SpinConstrainTest : public testing::Test
{
    std::vector<int> nw = {13};
    int nlocal;
    // get sc instance
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, Foo)
{
}
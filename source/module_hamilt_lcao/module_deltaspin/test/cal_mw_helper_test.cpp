#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "module_cell/klist.h"
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
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>>& sc
        = SpinConstrain<std::complex<double>>::getScInstance();
};

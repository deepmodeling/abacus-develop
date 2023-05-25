#include "gtest/gtest.h"

#define private public
#include "module_elecstate/module_charge/charge.h"

namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
double get_ucell_omega()
{
    return 500.0;
}
void Set_GlobalV_Default()
{
    GlobalV::NSPIN = 1;
    GlobalV::device_flag = "cpu";
    GlobalV::precision_flag = "double";
}
} // namespace elecstate

/************************************************
 *  unit test of module_charge/charge.cpp
 ***********************************************/

/**
 * - Tested Functions:
 */

class ChargeTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        elecstate::Set_GlobalV_Default();
    }
    void TearDown() override
    {
    }
};

TEST_F(ChargeTest, Constructor)
{
}

#undef private
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_general/module_surchem/surchem.h"

#include "gtest/gtest.h"

namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}

PotBase* Potential::get_pot_type(const std::string& pot_type)
{
    return new PotBase;
}
} // namespace elecstate

/************************************************
 *  unit test of potential_new.cpp
 ***********************************************/

/**
 * - Tested Functions:
 */

class PotentialNewTest : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
        // set up test data
        // set up test data
        // set up test data
    }

    virtual void TearDown()
    {
        // do nothing
    }
};

TEST_F(PotentialNewTest, Constructor)
{
    // test constructor
}
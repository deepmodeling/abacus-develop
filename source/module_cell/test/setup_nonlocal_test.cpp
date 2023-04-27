#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_cell/setup_nonlocal.h"

/************************************************
 *  unit test of class InfoNonlocal
 ***********************************************/

/**
 * - Tested Functions:
 *   - InfoNonlocal::InfoNonlocal()
 *   - InfoNonlocal::~InfoNonlocal()
*/

class NonlocalTest : public testing::Test
{
protected:
                  InfoNonlocal* infoNonlocal;
                  NonlocalTest()
                  {
                                    infoNonlocal = new InfoNonlocal();
                  }
                  ~NonlocalTest()
                  {
                                    delete infoNonlocal;
                  }
};

TEST_F(NonlocalTest, InfoNonlocal)
{
                  EXPECT_EQ(infoNonlocal->nprojmax, 0);
                  EXPECT_EQ(infoNonlocal->rcutmax_Beta, 0.0);
}
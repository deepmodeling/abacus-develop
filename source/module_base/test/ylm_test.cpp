#include "../ylm.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of class ylm
 ***********************************************/

/**
 * - Tested Functions:
 *   - ZEROS
 *     - set all elements of a double float array to zero
 * */

class ylmTest : public testing::Test
{
};

TEST_F(ylmTest,Zeros)
{
    double aaaa[100];
    ModuleBase::Ylm::ZEROS(aaaa,100); 
    for(int i = 0; i < 100; i++)
	{
        EXPECT_EQ(aaaa[i],0.0);
	}
}

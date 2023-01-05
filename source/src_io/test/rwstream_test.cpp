#include "gtest/gtest.h"
#include "gmock/gmock.h"
/************************************************
 *  unit test of rwstream.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Rwstream()
 *     - Open a binary file
 *   - close()
 *     - Close a binary file
 */

#include "../rwstream.h"

class RwstreamTest : public testing::Test
{
protected:
};

TEST_F(RwstreamTest, Close)
{
	int a1,a2;
	Rwstream wfc("wfc","w");
	wfc << 10 << 100;
	wfc.close();
	wfc.open("wfc","r");
	wfc >> a1;
	wfc >> a2;
	EXPECT_EQ(a1,10);
	EXPECT_EQ(a2,100);
	remove("wfc");
}

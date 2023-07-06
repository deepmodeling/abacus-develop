#include "gtest/gtest.h"
#include "gmock/gmock.h"
/************************************************
 *  unit test of binstream.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Binstream()
 *     - Open a binary file
 *   - close()
 *     - Close a binary file
 */

#include "../binstream.h"

class BinstreamTest : public testing::Test
{
protected:
};

TEST_F(BinstreamTest, variable)
{
	int a1,a2,a3;
	Binstream wfc("wfc","w");
	EXPECT_EQ(!wfc, false);
	wfc << 10 ;
	wfc.close();
	EXPECT_EQ(bool(wfc), false);
	wfc.open("wfc","a");
	EXPECT_EQ(bool(wfc), true);
	wfc << 100;
	wfc.close();
	EXPECT_EQ(!wfc, true);

	wfc.open("wfc","r");
	EXPECT_EQ(bool(wfc), true);
	wfc >> a1;
	wfc >> a2;
	testing::internal::CaptureStdout();
    EXPECT_EXIT(wfc >> a3, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Binstream"));
	EXPECT_EQ(bool(wfc), true);
	wfc.close();
	EXPECT_EQ(bool(wfc), false);

	EXPECT_EQ(a1,10);
	EXPECT_EQ(a2,100);

	wfc.open("wfc","w");
	testing::internal::CaptureStdout();
    EXPECT_EXIT(wfc >> a3, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Binstream"));
	remove("wfc");

	wfc.open("wfc", "r");
	EXPECT_EQ(bool(wfc), false); //If file is not open, return false.

	Binstream *p = new Binstream("wfc" , "r");
	delete p;
}

TEST_F(BinstreamTest, array)
{
	
	int a[10], b[11];
	for(int i = 0 ; i < 10 ; ++i) a[i] = i;
	Binstream wwfc("wfc","w");
	wwfc.write(a, 10);
	wwfc.close();

	Binstream rwfc("wfc","r");
	testing::internal::CaptureStdout();
    EXPECT_EXIT(rwfc.read(b,11);, ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output,testing::HasSubstr("Binstream"));
	rwfc.close();
	rwfc.open("wfc","r");
	rwfc.read(b, 10);
	rwfc.close();
	remove("wfc");

	for(int i = 0 ; i < 10 ; ++i)
	{
		EXPECT_EQ(a[i], b[i]);
	}

	wwfc.open("wfc", "w");
}
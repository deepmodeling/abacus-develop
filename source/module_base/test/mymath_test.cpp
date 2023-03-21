#include "../mymath.h"
#include "gtest/gtest.h"
#include <iostream>
#include <algorithm>

/************************************************
 *  unit test of heap sort functions
 ***********************************************/

/**
 * - Tested Functions:
 *   - Heapsort
 *     - heapsort function in mymath.cpp
 *   - Hpsort
 *     - hpsort function in mymath.cpp
 *
 */

class MymathTest : public testing::Test
{
protected:
	int number;
};

TEST_F(MymathTest,Heapsort)
{
	number = 5;
	double rr[number];
	int index[number];
	rr[0] = 10;
	rr[1] = 9;
	rr[2] = 8;
	rr[3] = 7;
	rr[4] = 6;
	index[0] = 0;
	//for (int i=0;i<number;i++)
	//	std::cout << i << " " << rr[i] << std::endl;
	ModuleBase::heapsort(number,rr,index);
	for (int i=0;i<number-1;i++)
	{
		//std::cout << i << " " << rr[i] << std::endl;
		EXPECT_LT(rr[i],rr[i+1]);
	}
}

TEST_F(MymathTest,Hpsort)
{
	number = 5;
	double rr[number];
	int index[number];
	rr[0] = 10;
	rr[1] = 9;
	rr[2] = 8;
	rr[3] = 7;
	rr[4] = 6;
	index[0] = 0;
	//for (int i=0;i<number;i++)
	//	std::cout << i << " " << rr[i] << std::endl;
	ModuleBase::hpsort(number,rr,index);
	//std::sort(rr,rr+number);
	for (int i=0;i<number-1;i++)
	{
		//std::cout << i << " " << rr[i] << std::endl;
		EXPECT_LT(rr[i],rr[i+1]);
	}
}

#include "../intarray.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of class IntArray
 ***********************************************/

/**
 * - Tested Functions:
 *   - Construct
 *     - construct an int array (2 to 6 dimensions)
 *   - Creat 
 *     - create an int array (2 to 6 dimensions)
 *   - GetArrayCount
 *     - get the total number of int array created
 *   - GetSize
 *     - get the total size of an int array
 *   - GetDim
 *     - get the dimension of an int array
 *   - ZeroOut
 *     - set all elements of an int array to zero
 *   - GetBound
 *     - get the size of each dimension of an int array
 *   - ArrayEqReal
 *     - set all value of an array to an int number
 *   - ArrayEqArray
 *     - equal an intarray to another intarray
 *   - Parentheses
 *     - access element by using operator"()"
 *   - ConstParentheses
 *     - access element by using "()" through pointer
 *     - without changing its elements
 *   -IntArrayAlloc
 *     - Warning of integer array allocation error
 */

namespace ModuleBase
{
void IntArrayAlloc();
}

class IntArrayTest : public testing::Test
{
protected:
	ModuleBase::IntArray a2, a3, a4, a5, a6;
	int aa = 11;
	int bb = 1;
	int count0;
	int count1;
	const int zero = 0;
};

TEST_F(IntArrayTest,GetArrayCount)
{
	count0 = ModuleBase::IntArray::getArrayCount();
	ModuleBase::IntArray c3, c4;
	count1 = ModuleBase::IntArray::getArrayCount();
	EXPECT_EQ((count1-count0),2);
}

TEST_F(IntArrayTest,Construct)
{
	ModuleBase::IntArray x2(1,5);
	ModuleBase::IntArray x3(1,5,3);
	ModuleBase::IntArray x4(1,7,3,4);
	ModuleBase::IntArray x5(1,5,3,8,2);
	ModuleBase::IntArray x6(1,7,3,4,3,2);
	EXPECT_EQ(x2.getSize(),5);
	EXPECT_EQ(x3.getSize(),15);
	EXPECT_EQ(x4.getSize(),84);
	EXPECT_EQ(x5.getSize(),240);
	EXPECT_EQ(x6.getSize(),504);
}

TEST_F(IntArrayTest,Create)
{
	a2.create(2,1);
	a3.create(3,2,1);
	a4.create(4,3,2,1);
	a5.create(5,4,3,2,1);
	a6.create(6,5,4,3,2,1);
	EXPECT_EQ(a2.getSize(),2);
	EXPECT_EQ(a3.getSize(),6);
	EXPECT_EQ(a4.getSize(),24);
	EXPECT_EQ(a5.getSize(),120);
	EXPECT_EQ(a6.getSize(),720);
}


TEST_F(IntArrayTest,GetSize)
{
	ModuleBase::IntArray x3(1,5,3);
	ModuleBase::IntArray x4(1,7,3,4);
	EXPECT_EQ(x3.getSize(),15);
	EXPECT_EQ(x4.getSize(),84);
}

TEST_F(IntArrayTest,GetDim)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	EXPECT_EQ(a2.getDim(),2);
	EXPECT_EQ(a3.getDim(),3);
	EXPECT_EQ(a4.getDim(),4);
	EXPECT_EQ(a5.getDim(),5);
	EXPECT_EQ(a6.getDim(),6);
}

TEST_F(IntArrayTest,ZeroOut)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	a2.zero_out();
	a3.zero_out();
	a4.zero_out();
	a5.zero_out();
	a6.zero_out();
	for (int i=0;i<a2.getSize();++i)
	{
		EXPECT_EQ(a2.ptr[i],zero);
    	}
	for (int i=0;i<a3.getSize();++i)
	{
		EXPECT_EQ(a3.ptr[i],zero);
    	}
	for (int i=0;i<a4.getSize();++i)
	{
		EXPECT_EQ(a4.ptr[i],zero);
    	}
	for (int i=0;i<a5.getSize();++i)
	{
		EXPECT_EQ(a5.ptr[i],zero);
    	}
	for (int i=0;i<a6.getSize();++i)
	{
		EXPECT_EQ(a6.ptr[i],zero);
    	}
}

TEST_F(IntArrayTest,GetBound)
{
	// 2d
	a2.create(2,3);
	EXPECT_EQ(a2.getBound1(),2);
	EXPECT_EQ(a2.getBound2(),3);
	// 3d
	a3.create(3,5,1);
	EXPECT_EQ(a3.getBound1(),3);
	EXPECT_EQ(a3.getBound2(),5);
	EXPECT_EQ(a3.getBound3(),1);
	EXPECT_EQ(a3.getBound4(),1);

	ModuleBase::IntArray x3(1,5,3);
	EXPECT_EQ(x3.getBound1(),1);
	EXPECT_EQ(x3.getBound2(),5);
	EXPECT_EQ(x3.getBound3(),3);
	EXPECT_EQ(x3.getBound4(),0);
	EXPECT_EQ(x3.getBound5(),0);
	EXPECT_EQ(x3.getBound6(),0);
	// 4d
	a4.create(4,3,7,1);
	EXPECT_EQ(a4.getBound1(),4);
	EXPECT_EQ(a4.getBound2(),3);
	EXPECT_EQ(a4.getBound3(),7);
	EXPECT_EQ(a4.getBound4(),1);
	// 5d
	a5.create(5,4,1,2,1);
	EXPECT_EQ(a5.getBound1(),5);
	EXPECT_EQ(a5.getBound2(),4);
	EXPECT_EQ(a5.getBound3(),1);
	EXPECT_EQ(a5.getBound4(),2);
	EXPECT_EQ(a5.getBound5(),1);
	// 6d
	a6.create(6,5,9,3,2,1);
	EXPECT_EQ(a6.getBound1(),6);
	EXPECT_EQ(a6.getBound2(),5);
	EXPECT_EQ(a6.getBound3(),9);
	EXPECT_EQ(a6.getBound4(),3);
	EXPECT_EQ(a6.getBound5(),2);
	EXPECT_EQ(a6.getBound6(),1);
}

TEST_F(IntArrayTest,ArrayEqInt)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	a2 = aa;
	a3 = bb;
	a4 = aa+1;
	a5 = bb -1;
	a6 = 10;
	for (int i=0;i<a2.getSize();++i)
	{
		EXPECT_EQ(a2.ptr[i],aa);
    	}
	for (int i=0;i<a3.getSize();++i)
	{
		EXPECT_EQ(a3.ptr[i],bb);
    	}
	for (int i=0;i<a4.getSize();++i)
	{
		EXPECT_EQ(a4.ptr[i],aa+1);
    	}
	for (int i=0;i<a5.getSize();++i)
	{
		EXPECT_EQ(a5.ptr[i],bb-1);
    	}
	for (int i=0;i<a6.getSize();++i)
	{
		EXPECT_EQ(a6.ptr[i],10);
    	}
}

TEST_F(IntArrayTest,ArrayEqArray)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	a2 = aa;
	a3 = bb;
	a4 = aa+1;
	a5 = bb -1;
	a6 = 10;
	ModuleBase::IntArray x2(2,3);
	ModuleBase::IntArray x3(3,5,1);
	ModuleBase::IntArray x4(4,3,7,1);
	ModuleBase::IntArray x5(5,4,1,2,1);
	ModuleBase::IntArray x6(6,5,9,3,2,1);
	x2 = a2; x3 = a3; x4 = a4; x5 = a5; x6 = a6;
	for (int i=0;i<x2.getSize();++i)
	{
		EXPECT_EQ(x2.ptr[i],aa);
    	}
	for (int i=0;i<x3.getSize();++i)
	{
		EXPECT_EQ(x3.ptr[i],bb);
    	}
	for (int i=0;i<x4.getSize();++i)
	{
		EXPECT_EQ(x4.ptr[i],aa+1);
    	}
	for (int i=0;i<x5.getSize();++i)
	{
		EXPECT_EQ(x5.ptr[i],bb-1);
    	}
	for (int i=0;i<x6.getSize();++i)
	{
		EXPECT_EQ(x6.ptr[i],10);
    	}
}

TEST_F(IntArrayTest,Parentheses)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	a2 = aa;
	a3 = bb;
	a4 = aa+1;
	a5 = bb -1;
	a6 = 10;
	EXPECT_EQ(a2(0,0),aa);
	EXPECT_EQ(a2(1,2),aa);
	EXPECT_EQ(a3(0,0,0),bb);
	EXPECT_EQ(a3(2,4,0),bb);
	EXPECT_EQ(a4(0,0,0,0),aa+1);
	EXPECT_EQ(a4(3,2,6,0),aa+1);
	EXPECT_EQ(a5(0,0,0,0,0),bb-1);
	EXPECT_EQ(a5(4,3,0,1,0),bb-1);
	EXPECT_EQ(a6(0,0,0,0,0,0),10);
	EXPECT_EQ(a6(5,4,8,2,1,0),10);
}

TEST_F(IntArrayTest,ConstParentheses)
{
	a2.create(2,3);
	a3.create(3,5,1);
	a4.create(4,3,7,1);
	a5.create(5,4,1,2,1);
	a6.create(6,5,9,3,2,1);
	a2 = aa;
	a3 = bb;
	a4 = aa+1;
	a5 = bb -1;
	a6 = 10;
	const ModuleBase::IntArray *x2(&a2);
	const ModuleBase::IntArray *x3(&a3);
	const ModuleBase::IntArray *x4(&a4);
	const ModuleBase::IntArray *x5(&a5);
	const ModuleBase::IntArray *x6(&a6);
	EXPECT_EQ((*x2)(0,0),aa);
	EXPECT_EQ((*x2)(1,2),aa);
	EXPECT_EQ((*x3)(0,0,0),bb);
	EXPECT_EQ((*x3)(2,4,0),bb);
	EXPECT_EQ((*x4)(0,0,0,0),aa+1);
	EXPECT_EQ((*x4)(3,2,6,0),aa+1);
	EXPECT_EQ((*x5)(0,0,0,0,0),bb-1);
	EXPECT_EQ((*x5)(4,3,0,1,0),bb-1);
	EXPECT_EQ((*x6)(0,0,0,0,0,0),10);
	EXPECT_EQ((*x6)(5,4,8,2,1,0),10);
}

TEST_F(IntArrayTest,Alloc)
{
	std::string output;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::IntArrayAlloc(), ::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Allocation error for IntArray"));
}

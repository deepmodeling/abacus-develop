#include "../realarray.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
/************************************************
 *  unit test of class realArray
 ***********************************************/

/**
 * - Tested Functions:
 *   - GetArrayCount
 *     - get the total number of real array created
 *   - Construct
 *     - construct real arrays (3 or 4 dimensions)
 *   - Creat
 *     - create a real array (3 or 4 dimensions)
 *   - GetSize
 *     - get the total size of a real array
 *   - GetDim
 *     - get the dimension of a real array
 *   - ZeroOut
 *     - set all elements of a real array to zero
 *   - GetBound
 *     - get the size of each dimension of a real array
 *   - ArrayEqReal
 *     - set all value of an array to a double float
 *   - ArrayEqArray
 *     - equal a realarray to another one
 *   - Parentheses
 *     - access element by using operator"()"
 *   - ConstParentheses
 *     - access element by using "()" through pointer
 *     - without changing its elements
 *   - realArrayAlloc
 *     - output the warning message when allocation error for realArray
 *     - occurs
 *   - zeros
 *     - set all elements of a 1_d array to zero
 */

class realArrayTest : public testing::Test
{
protected:
	ModuleBase::realArray a3, a4, b3, b4;
	double aa = 11.0;
	double bb = 1.0;
	int count0 = 0;
	int count1 = 0;
	const double zero = 0.0;
};

namespace ModuleBase
{
void realArrayAlloc();
}

TEST_F(realArrayTest,GetArrayCount)
{
	count0 = ModuleBase::realArray::getArrayCount();
	ModuleBase::realArray c3, c4;
	count1 = ModuleBase::realArray::getArrayCount();
	EXPECT_EQ((count1-count0),2);
}

TEST_F(realArrayTest,Construct)
{
	ModuleBase::realArray x3(1,5,3);
	ModuleBase::realArray xp3(x3);
	ModuleBase::realArray x4(1,7,3,4);
	ModuleBase::realArray xp4(x4);
	EXPECT_EQ(x3.getSize(),15);
	EXPECT_EQ(xp3.getSize(),15);
	EXPECT_EQ(x4.getSize(),84);
	EXPECT_EQ(xp4.getSize(),84);
}

TEST_F(realArrayTest,Create)
{
	a3.create(1,2,3);
	a4.create(1,2,3,4);
	EXPECT_EQ(a3.getSize(),6);
	EXPECT_EQ(a4.getSize(),24);
}

TEST_F(realArrayTest,GetSize)
{
	ModuleBase::realArray a3(1,5,3);
	//std::cout<< &a3 << &(this->a3) <<std::endl;
	ModuleBase::realArray a4(1,7,3,4);
	EXPECT_EQ(a3.getSize(),15);
	EXPECT_EQ(a4.getSize(),84);
}

TEST_F(realArrayTest,GetDim)
{
	ModuleBase::realArray a3(1,5,3);
	ModuleBase::realArray a4(1,7,3,4);
	EXPECT_EQ(a3.getDim(),3);
	EXPECT_EQ(a4.getDim(),4);
}

TEST_F(realArrayTest,ZeroOut)
{
	a3.create(1,2,3);
	a4.create(1,2,3,4);
	a3.zero_out();
	a4.zero_out();
	for (int i=0;i<a3.getSize();++i)
	{
		EXPECT_EQ(a3.ptr[i],zero);
    	}
	for (int i=0;i<a4.getSize();++i)
	{
		EXPECT_EQ(a4.ptr[i],zero);
    	}
}

TEST_F(realArrayTest,GetBound)
{
	ModuleBase::realArray a3(1,5,3);
	ModuleBase::realArray a4(1,7,3,4);
	EXPECT_EQ(a3.getBound1(),1);
	EXPECT_EQ(a3.getBound2(),5);
	EXPECT_EQ(a3.getBound3(),3);
	EXPECT_EQ(a4.getBound1(),1);
	EXPECT_EQ(a4.getBound2(),7);
	EXPECT_EQ(a4.getBound3(),3);
	EXPECT_EQ(a4.getBound4(),4);
}

TEST_F(realArrayTest,ArrayEqReal)
{
	ModuleBase::realArray a3(1,5,3);
	ModuleBase::realArray a4(1,7,3,4);
	a3 = aa;
	a4 = bb;
	for (int i=0;i<a3.getSize();++i)
	{
		EXPECT_EQ(a3.ptr[i],aa);
    	}
	for (int i=0;i<a4.getSize();++i)
	{
		EXPECT_EQ(a4.ptr[i],bb);
    	}
}

TEST_F(realArrayTest,ArrayEqArray)
{
	ModuleBase::realArray a3(1,5,3);
	ModuleBase::realArray a4(1,7,3,4);
	a3 = aa;
	a4 = bb;
	b3.create(1,5,3);
	b4.create(1,7,3,4);
	b3 = a3;
	b4 = a4;
	for (int i=0;i<b3.getSize();++i)
	{
		EXPECT_EQ(b3.ptr[i],aa);
    	}
	for (int i=0;i<b4.getSize();++i)
	{
		EXPECT_EQ(b4.ptr[i],bb);
    	}
}

TEST_F(realArrayTest,Parentheses)
{
	ModuleBase::realArray a3(1,5,3);
	ModuleBase::realArray a4(1,7,3,4);
	a3 = aa;
	a4 = bb;
	EXPECT_EQ(a3(0,0,0),aa);
	EXPECT_EQ(a3(0,4,2),aa);
	EXPECT_NE(a3(0,0,0),bb);
	EXPECT_NE(a3(0,4,2),bb);
	EXPECT_EQ(a4(0,0,0,0),bb);
	EXPECT_EQ(a4(0,6,2,3),bb);
	EXPECT_NE(a4(0,0,0,0),aa);
	EXPECT_NE(a4(0,6,2,3),aa);
}

TEST_F(realArrayTest,ConstParenthesesExpert)
{
	ModuleBase::realArray a3(1,5,3);
	a3 = aa;
	const ModuleBase::realArray b3 (a3);
	EXPECT_EQ(b3(0,0,0),aa);
	EXPECT_EQ(b3(0,4,2),aa);
	ModuleBase::realArray a4(1,5,3,2);
	a4 = bb;
	const ModuleBase::realArray b4 (a4);
	EXPECT_EQ(b4(0,1,2,1),bb);
	EXPECT_EQ(b4(0,0,0,0),bb);
}

TEST_F(realArrayTest,Alloc)
{
	std::string output;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::realArrayAlloc(), ::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Allocation error for realArray"));
}

//use a pointer to initialize a double array and then call zeros
TEST_F(realArrayTest,Zeros)
{
        double *p = new double [100];
        ModuleBase::zeros<double>(p,100);
        for(int i = 0 ; i < 100; i ++)
        {
                EXPECT_DOUBLE_EQ(p[i],0.0);
        }
        delete [] p;

}

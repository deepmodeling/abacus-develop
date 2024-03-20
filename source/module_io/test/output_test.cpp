#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <random>
/************************************************
 *  unit test of output.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - printr1_d()
 *     - Print out one dimensional array
 *   - printV3()
 *     - Print out ModuleBase::Vector3
 *   - printv31_d()
 *     - Print out one dimensional array of ModuleBase::Vector3
 *   - printM3()
 *     - Print out ModuleBase::Matrix3
 *   - printr3_d()
 *     - Print out three dimensional ModuleBase::realArray
 *   - printrm()
 *     - Print out ModuleBase::matrix
 */

template<typename T>
T* get_simple_array(int num)
{
	T* rand_array = new T[num]();
	assert(num>0);
	for (int i=0;i<num;i++)
	{
	  rand_array[i] = static_cast<int>(i);
	}
	return rand_array;
}

#include "../output.h"

class OutputTest : public testing::Test
{
protected:
	std::ofstream ofs;
	std::ifstream ifs;
	std::string output_str;
	int num = 5;
};

TEST_F(OutputTest, Printr1_d)
{
	ofs.open("running.log");
	int *simple_array = NULL;
	simple_array = get_simple_array<int>(num);
	output::printr1_d(ofs,"simple",simple_array,num);
	ofs.close();
	ifs.open("running.log");
	getline(ifs,output_str);
	getline(ifs,output_str);
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("simple  n1 = 5"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("0           1           2           3           4"));
	ifs.close();
	remove("running.log");
}

TEST_F(OutputTest, PrintV3)
{
	ofs.open("running.log");
	ModuleBase::Vector3<double> u (3.0,4.0,5.0);
	output::printV3(ofs,u);
	ofs.close();
	ifs.open("running.log");
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("3                 4                 5"));
	ifs.close();
	remove("running.log");
}

TEST_F(OutputTest, Printv31_d)
{
	ofs.open("running.log");
	ModuleBase::Vector3<double>* V = new ModuleBase::Vector3<double>[3];
	V[0].set(0.0,1.0,2.0);
	V[1].set(3.0,4.0,5.0);
	V[2].set(6.6,7.0,8.0);
	output::printv31_d(ofs,"vector3 array",V,3);
	ofs.close();
	ifs.open("running.log");
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("vector3 array Dimension = 3"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("0                 1                 2"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("3                 4                 5"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("6.6                 7                 8"));
	ifs.close();
	remove("running.log");
}

TEST_F(OutputTest, PrintM3)
{
	ModuleBase::Matrix3 mb(1,2,3,4,5,6,7,8,9);
	testing::internal::CaptureStdout();
	output::printM3("Matrix3",mb);
	output_str = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output_str,testing::HasSubstr("Matrix3"));
	EXPECT_THAT(output_str,testing::HasSubstr("1                   2                   3"));
	EXPECT_THAT(output_str,testing::HasSubstr("4                   5                   6"));
	EXPECT_THAT(output_str,testing::HasSubstr("7                   8                   9"));
}

TEST_F(OutputTest, Printr3_d)
{
	ofs.open("running.log");
	ModuleBase::realArray a3;
	a3.create(1,1,4);
	a3.zero_out();
	a3(0,0,0) = 0.0;
	a3(0,0,1) = 1.0;
	a3(0,0,2) = 2.0;
	a3(0,0,3) = 3.0;
	output::printr3_d(ofs,"realArray3",a3);
	ofs.close();
	ifs.open("running.log");
	getline(ifs,output_str);
	getline(ifs,output_str);
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("realArray3  b1 = 1 b2 = 1 b3 = 4"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("0                 1                 2                 3"));
	ifs.close();
	remove("running.log");
}

TEST_F(OutputTest, Printrm)
{
	ofs.open("running.log");
	ModuleBase::matrix m23a;
        m23a.create(2,3);
        for (int i=1;i<=6;++i) {m23a.c[i-1] = i*1.0;}
	output::printrm(ofs,"Matrix",m23a);
	ofs.close();
	ifs.open("running.log");
	getline(ifs,output_str);
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("Matrix nr=2 nc=3"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("1           2           3"));
	getline(ifs,output_str);
	EXPECT_THAT(output_str,testing::HasSubstr("4           5           6"));
	ifs.close();
	remove("running.log");
	// on screen
	testing::internal::CaptureStdout();
	output::printrm("Matrix",m23a);
	output_str = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output_str,testing::HasSubstr("Matrix nr=2 nc=3"));
	EXPECT_THAT(output_str,testing::HasSubstr("1           2           3"));
	EXPECT_THAT(output_str,testing::HasSubstr("4           5           6"));
}

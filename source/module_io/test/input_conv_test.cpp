#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/input_conv.h"
#include "module_base/global_variable.h"
#include "for_testing_input_conv.h"

/************************************************
 *  unit test of input_conv.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Convert()
 */

#define private public
#include "module_io/input.h"

class InputConvTest : public testing::Test
{
protected:
	std::string output;
};

TEST_F(InputConvTest, Conv)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::stru_file,INPUT.stru_file);
}

TEST_F(InputConvTest,ParseExpressionDouble)
{
	std::vector<double> vec;
	std::string input = "2*3.5 1.2 0.5";
	Input_Conv::parse_expression(input,vec);
	EXPECT_EQ(vec.size(),4);
	EXPECT_DOUBLE_EQ(vec[0],3.5);
	EXPECT_DOUBLE_EQ(vec[1],3.5);
	EXPECT_DOUBLE_EQ(vec[2],1.2);
	EXPECT_DOUBLE_EQ(vec[3],0.5);
}

#undef private

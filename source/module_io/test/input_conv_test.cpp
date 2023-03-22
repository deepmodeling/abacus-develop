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

#undef private

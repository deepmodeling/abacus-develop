#include "../tool_check.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <cstdio>
#include <fstream>

/************************************************
 *  unit test of functions in tool_check.h
 ***********************************************/

/**
 * - Tested Function
 *   - ModuleBase::CHECK_NAME
 *    - check the next input from ifs is std::string
 *
 *   - ModuleBase::CHECK_INT
 *    - check the next input from ifs is int
 *
 *   - ModuleBase::CHECK_DOUBLE
 *    - check the next input from ifs is double
 *
 *   - ModuleBase::CHECK_STRING
 *    - check the next input from ifs is string
 */

class ToolCheckTest : public testing::Test
{
  protected:
	std::ofstream ofs;
	std::ifstream ifs;
	// define std::string, int, double variables
	std::string name = "abaqus";
	int ecut = 100;
	double occupation = 0.23;
	std::string caltype = "nscf";
	// quit is the swith to control performance of function
	bool quit = false;
	// for capturing stdout
	std::string output = "";
	void TearDown()
	{
		remove("tmp");
	}

};

TEST_F(ToolCheckTest, Name)
{
	ofs.open("tmp");
	// double input to check continus check function
	ofs << name << std::endl;
	ofs << name << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
	testing::internal::CaptureStdout();
	ModuleBase::CHECK_NAME(ifs, "abacus", quit);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not match"));
	// quit check: quit = false
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::CHECK_NAME(ifs, "abacus"), ::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("NOTICE"));
	ifs.close();
}

TEST_F(ToolCheckTest, Int)
{
	ofs.open("tmp");
	// double input to check continus check function
	ofs << ecut << std::endl;
	ofs << ecut << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check
	testing::internal::CaptureStdout();
	ModuleBase::CHECK_INT(ifs, 80, quit);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not match"));
	// quit check: quit = false
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::CHECK_INT(ifs, 80), ::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("NOTICE"));
	ifs.close();
}

TEST_F(ToolCheckTest, Double)
{
	ofs.open("tmp");
	// double input to check continus check function
	ofs << occupation << std::endl;
	ofs << occupation << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check: quit = false
	testing::internal::CaptureStdout();
	ModuleBase::CHECK_DOUBLE(ifs, 0.23002, quit);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not match"));
	// quit check
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::CHECK_DOUBLE(ifs, 0.22998), ::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("NOTICE"));
	ifs.close();
}

TEST_F(ToolCheckTest, String)
{
	ofs.open("tmp");
	// double input to check continus check function
	ofs << caltype << std::endl;
	ofs << caltype << std::endl;
	ofs.close();
	ifs.open("tmp");
	// non-quit check: quit=false
	testing::internal::CaptureStdout();
	ModuleBase::CHECK_STRING(ifs, "scf", quit);
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not match"));
	// quit check
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::CHECK_STRING(ifs, "scf"), ::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("NOTICE"));
	ifs.close();
}

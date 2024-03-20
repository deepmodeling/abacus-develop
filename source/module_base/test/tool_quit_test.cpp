#include "../tool_quit.h"
#include "../global_variable.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of functions in tool_quit.h
 ***********************************************/

/**
 * - Tested Function
 *   - ModuleBase::WARNING
 *     - print out warning info. in warning.log
 *
 *   - ModuleBase::QUIT
 *     - close log files and exit
 *
 *   - ModuleBase::WARNING_QUIT
 *     - combine the above 2 functions
 */

class ToolQuitTest : public testing::Test
{
  protected:
	std::ifstream ifs;
	// for capturing output on screen
	std::string output;
	void SetUp()
	{
		GlobalV::ofs_warning.open("warning.log");
		GlobalV::ofs_running.open("running.log");
		GlobalV::global_out_dir = "OUT/";
	}
	void TearDown()
	{
		remove("warning.log");
		remove("running.log");
	}
};


TEST_F(ToolQuitTest,warning)
{
	ModuleBase::WARNING("INPUT","bad input parameter");
	GlobalV::ofs_warning.close();
	ifs.open("warning.log");
	getline(ifs,output);
	// test output in warning.log file
	EXPECT_THAT(output,testing::HasSubstr("warning"));
	ifs.close();
}

TEST_F(ToolQuitTest,quit)
{
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::QUIT(), ::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	// test output on screen
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
}

TEST_F(ToolQuitTest,quit_with_ret)
{
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::QUIT(1), ::testing::ExitedWithCode(1), "");
	output = testing::internal::GetCapturedStdout();
	// test output on screen
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
}

// use EXPECT_EXIT to test exit codes
TEST_F(ToolQuitTest,warningquit)
{
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::WARNING_QUIT("INPUT","bad input parameter"),
			::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	// test output on screening
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	GlobalV::ofs_warning.close();
	GlobalV::ofs_running.close();
	ifs.open("warning.log");
	getline(ifs,output);
	// test output in warning.log file
	EXPECT_THAT(output,testing::HasSubstr("warning"));
	ifs.close();
	ifs.open("running.log");
	getline(ifs,output);
	// test output in running.log file
	EXPECT_THAT(output,testing::HasSubstr("!!!!!!!"));
	ifs.close();
}

// use EXPECT_EXIT to test exit codes
TEST_F(ToolQuitTest,warningquit_with_ret)
{
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ModuleBase::WARNING_QUIT("INPUT","bad input parameter",1),
			::testing::ExitedWithCode(1), "");
	output = testing::internal::GetCapturedStdout();
	// test output on screening
	EXPECT_THAT(output,testing::HasSubstr("TIME STATISTICS"));
	GlobalV::ofs_warning.close();
	GlobalV::ofs_running.close();
	ifs.open("warning.log");
	getline(ifs,output);
	// test output in warning.log file
	EXPECT_THAT(output,testing::HasSubstr("warning"));
	ifs.close();
	ifs.open("running.log");
	getline(ifs,output);
	// test output in running.log file
	EXPECT_THAT(output,testing::HasSubstr("!!!!!!!"));
	ifs.close();
}
// use __MPI to activate parallel environment
#ifdef __MPI
int main(int argc, char **argv)
{

	MPI_Init(&argc,&argv);

	testing::InitGoogleTest(&argc,argv);
	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif


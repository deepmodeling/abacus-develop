#include "module_io/parse_args.h"

#include "gtest/gtest.h"
#include "module_io/read_input.h"
#include "version.h"

bool ModuleIO::ReadInput::check_mode = false;

TEST(ParseArgsTest, OutVersionTest)
{
    // Test case 1: no arguments
    char arg0[] = "test";
    char* argv[] = {arg0};
    int argc = 1;
    testing::internal::CaptureStdout();
    ModuleIO::parse_args(argc, argv);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ("", output);
    // No output expected

#ifdef VERSION
    std::string output_ref = "ABACUS version " + std::string(VERSION) + "\n";
#else
    std::string output_ref = "ABACUS version unknown\n";
#endif

    // Test case 2: --version argument
    char arg1[] = "--version";
    char* argv1[] = {arg0, arg1};
    argc = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::parse_args(argc, argv1), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);

    // Test case 3: -v argument
    char arg2[] = "-v";
    char* argv2[] = {arg0, arg2};
    argc = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::parse_args(argc, argv2), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);

    // Test case 4: -V argument
    char arg3[] = "-V";
    char* argv3[] = {arg0, arg3};
    argc = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::parse_args(argc, argv3), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);
}

TEST(ParseArgsTest, CheckInput)
{
    char arg0[] = "test";
    char arg1[] = "--check-input";
    char* argv[] = {arg0, arg1};
    int argc = 2;
    ModuleIO::parse_args(argc, argv);
    EXPECT_TRUE(ModuleIO::ReadInput::check_mode);
}
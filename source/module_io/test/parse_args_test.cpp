#include "gtest/gtest.h"
#include "module_io/parse_args.h"
#include "version.h"

TEST(ParseArgsTest, OutVersionTest)
{
    // Test case 1: no arguments
    char* argv[] = {"test"};
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
    char* argv1[] = {"test", "--version"};
    argc = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::parse_args(argc, argv1),::testing::ExitedWithCode(0),"");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);

    // Test case 3: -v argument
    char* argv2[] = {"test", "-v"};
    argc = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleIO::parse_args(argc, argv2),::testing::ExitedWithCode(0),"");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);

    // Test case 4: -V argument
    char* argv3[] = {"test", "-V"};
    argc = 2;
    testing::internal::CaptureStdout();
        EXPECT_EXIT(ModuleIO::parse_args(argc, argv3),::testing::ExitedWithCode(0),"");
    output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output_ref, output);
}
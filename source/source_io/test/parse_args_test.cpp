#include "source_io/parse_args.h"
#include "gtest/gtest.h"
#include "source_io/read_input.h"
#include "source_main/version.h"

// Already deal with Testing.cmake
// #include "build_info.h" 

// This file is modified by ZhouXY-PKU at 2025-12-01

// Refresh test status
class ParseArgsTest : public ::testing::Test {
protected:
    void SetUp() override {
        ModuleIO::ReadInput::check_mode = false;
    }
};

// Test for no argument
TEST_F(ParseArgsTest, NoArguments) {
    char arg0[] = "test";
    char* argv[] = {arg0};
    int argc = 1;

    testing::internal::CaptureStdout();
    ModuleIO::parse_args(argc, argv);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.empty()) << "Expected no output for no arguments.";
}

// Test for abacus version
TEST_F(ParseArgsTest, VersionFlags) {
#ifdef VERSION
    std::string output_ref = "ABACUS version " + std::string(VERSION) + "\n";
#else
    std::string output_ref = "ABACUS version unknown\n";
#endif

    std::vector<std::string> version_args = {"--version", "-v", "-V"};

    for (const auto& arg : version_args) {
        char arg0[] = "test";
        std::vector<char*> argv = {arg0, const_cast<char*>(arg.c_str())};
        int argc = argv.size();
        
        testing::internal::CaptureStdout();
        EXPECT_EXIT(
            { ModuleIO::parse_args(argc, argv.data()); },
            ::testing::ExitedWithCode(0),
            ""
        ) << "Failed for argument: " << arg;
        
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_EQ(output_ref, output) << "Output mismatch for argument: " << arg;
    }
}

// Test for abacus info
TEST_F(ParseArgsTest, InfoFlags) {
    std::vector<std::string> info_args = {"--info", "-i", "-I"};

    for (const auto& arg : info_args) {
        char arg0[] = "test";
        std::vector<char*> argv = {arg0, const_cast<char*>(arg.c_str())};
        int argc = argv.size();

        testing::internal::CaptureStdout();
        EXPECT_EXIT(
            { ModuleIO::parse_args(argc, argv.data()); },
            ::testing::ExitedWithCode(0),
            ""
        ) << "Failed for argument: " << arg;
        
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_TRUE(output.find("ABACUS Core") != std::string::npos) 
            << "Output mismatch for argument: " << arg << "\nCaptured output was: " << output;
    }
}

// Test for unavailable arguments
TEST_F(ParseArgsTest, UnknownArgument) {
    char arg0[] = "test";
    char arg1[] = "--nonexistent-option";
    char* argv[] = {arg0, arg1};
    int argc = 2;

   // testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(1),
        ""
    );
  //  std::string stdout_output = testing::internal::GetCapturedStdout();
    std::string stderr_output = testing::internal::GetCapturedStderr();
  //  std::string combined_output = stdout_output + stderr_output;
    EXPECT_TRUE(stderr_output.find("Usage: abacus") != std::string::npos)
        << "Output did not contain usage information.\n"
        << "Captured stderr was: '" << stderr_output << "'";
}

// Test for --check-input
TEST_F(ParseArgsTest, CheckInputFlag) {
    char arg0[] = "test";
    char arg1[] = "--check-input";
    char* argv[] = {arg0, arg1};
    int argc = 2;

    ModuleIO::parse_args(argc, argv);
    
    EXPECT_TRUE(ModuleIO::ReadInput::check_mode);
}

TEST_F(ParseArgsTest, PriorityVersionOverCheckInput) {
    char arg0[] = "test";
    char arg1[] = "--version";
    char arg2[] = "--check-input";
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("ABACUS version") != std::string::npos)
        << "Output did not contain version information.\nCaptured output was: " << output;
}

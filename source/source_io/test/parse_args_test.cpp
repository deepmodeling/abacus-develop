#include "source_io/parse_args.h"
#include "gtest/gtest.h"
#include "source_io/module_parameter/read_input.h"
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

    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(1),
        "Usage: abacus"
    ) << "Failed for unknown argument test.";
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

// Test for -h without argument (general help)
TEST_F(ParseArgsTest, HelpGeneralShort) {
    char arg0[] = "test";
    char arg1[] = "-h";
    char* argv[] = {arg0, arg1};
    int argc = 2;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("ABACUS") != std::string::npos);
    EXPECT_TRUE(output.find("Usage:") != std::string::npos);
}

// Test for --help without argument (general help)
TEST_F(ParseArgsTest, HelpGeneralLong) {
    char arg0[] = "test";
    char arg1[] = "--help";
    char* argv[] = {arg0, arg1};
    int argc = 2;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("ABACUS") != std::string::npos);
    EXPECT_TRUE(output.find("Usage:") != std::string::npos);
}

// Test for -h with known parameter
TEST_F(ParseArgsTest, HelpKnownParameter) {
    char arg0[] = "test";
    char arg1[] = "-h";
    char arg2[] = "calculation";
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("Parameter:") != std::string::npos);
    EXPECT_TRUE(output.find("calculation") != std::string::npos);
    EXPECT_TRUE(output.find("Type:") != std::string::npos);
    EXPECT_TRUE(output.find("Description:") != std::string::npos);
}

// Test for -h with unknown parameter (should show fuzzy suggestions)
TEST_F(ParseArgsTest, HelpUnknownParameter) {
    char arg0[] = "test";
    char arg1[] = "-h";
    char arg2[] = "ecutwf";  // Similar to ecutwfc
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(1),
        "Unknown parameter.*Did you mean"
    );
}

// Test for -s with results
TEST_F(ParseArgsTest, SearchWithResults) {
    char arg0[] = "test";
    char arg1[] = "-s";
    char arg2[] = "ecut";  // Should match ecutwfc, ecutrho, etc.
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("Found") != std::string::npos);
    EXPECT_TRUE(output.find("parameter(s) matching") != std::string::npos);
    EXPECT_TRUE(output.find("ecut") != std::string::npos);
}

// Test for --search with results
TEST_F(ParseArgsTest, SearchLongWithResults) {
    char arg0[] = "test";
    char arg1[] = "--search";
    char arg2[] = "basis";
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("Found") != std::string::npos);
    EXPECT_TRUE(output.find("parameter(s) matching") != std::string::npos);
}

// Test for -s with no results
TEST_F(ParseArgsTest, SearchWithNoResults) {
    char arg0[] = "test";
    char arg1[] = "-s";
    char arg2[] = "xyzabc123notfound";
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(1),
        "No parameters found matching"
    );
}

// Test for -s without query (error)
TEST_F(ParseArgsTest, SearchWithoutQuery) {
    char arg0[] = "test";
    char arg1[] = "-s";
    char* argv[] = {arg0, arg1};
    int argc = 2;

    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(1),
        "requires a search query"
    );
}

// Test for -h followed by a flag (should show general help, not treat flag as parameter)
TEST_F(ParseArgsTest, HelpFollowedByFlag) {
    char arg0[] = "test";
    char arg1[] = "-h";
    char arg2[] = "--version";
    char* argv[] = {arg0, arg1, arg2};
    int argc = 3;

    testing::internal::CaptureStdout();
    EXPECT_EXIT(
        { ModuleIO::parse_args(argc, argv); },
        ::testing::ExitedWithCode(0),
        ""
    );
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(output.find("Usage:") != std::string::npos);
    // Should show general help, not version
    EXPECT_TRUE(output.find("ABACUS version") == std::string::npos ||
                output.find("Usage:") != std::string::npos);
}


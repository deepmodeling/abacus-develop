#include "source_io/parse_command_line.h"

#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
/// Builds an argv array from strings and returns the token count
/// (excluding the terminating nullptr) via the out-parameter.
/// Lifetime note: `storage` and `strings` must outlive the parse call.
int make_argv(std::vector<std::string>& strings,
              std::vector<char*>& storage)
{
    strings.insert(strings.begin(), "abacus");
    storage.clear();
    storage.reserve(strings.size());
    for (size_t i = 0; i < strings.size(); ++i)
    {
        storage.push_back(&strings[i][0]);
    }
    storage.push_back(nullptr);
    return static_cast<int>(strings.size());
}
} // namespace

TEST(ParseCommandLineTest, ShortFormSingleVar)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"-p", (char*)"a", (char*)"1", nullptr};
    const ModuleIO::CommandLineArgs args
        = ModuleIO::parse_command_line(4, argv);
    ASSERT_EQ(args.vars.size(), 1u);
    EXPECT_EQ(args.vars.at("a"), "1");
}

TEST(ParseCommandLineTest, LongFormSingleVar)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"--parameter", (char*)"ecutwfc", (char*)"120",
                    nullptr};
    const ModuleIO::CommandLineArgs args
        = ModuleIO::parse_command_line(4, argv);
    EXPECT_EQ(args.vars.at("ecutwfc"), "120");
}

TEST(ParseCommandLineTest, MultipleVarsLastWins)
{
    std::vector<std::string> tokens
        = {"-p", "a", "1", "-p", "a", "2", "-p", "b", "x"};
    std::vector<char*> storage;
    const int argc = make_argv(tokens, storage);
    const ModuleIO::CommandLineArgs args
        = ModuleIO::parse_command_line(argc, storage.data());
    EXPECT_EQ(args.vars.at("a"), "2");
    EXPECT_EQ(args.vars.at("b"), "x");
}

TEST(ParseCommandLineTest, VarPlusCustomInput)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"-p",  (char*)"suffix", (char*)"TestRun",
                    (char*)"-in", (char*)"myINPUT",
                    nullptr};
    const ModuleIO::CommandLineArgs args
        = ModuleIO::parse_command_line(6, argv);
    EXPECT_EQ(args.vars.at("suffix"), "TestRun");
    EXPECT_EQ(args.input_file, "myINPUT");
}

TEST(ParseCommandLineTest, DefaultInputPath)
{
    char* argv[] = {(char*)"abacus", nullptr};
    const ModuleIO::CommandLineArgs args
        = ModuleIO::parse_command_line(1, argv);
    EXPECT_EQ(args.input_file, "INPUT");
    EXPECT_TRUE(args.vars.empty());
}

TEST(ParseCommandLineTest, MissingValueThrows)
{
    char* argv[] = {(char*)"abacus", (char*)"-p", (char*)"a", nullptr};
    EXPECT_THROW(ModuleIO::parse_command_line(3, argv), std::runtime_error);
}

TEST(ParseCommandLineTest, EmptyNameThrows)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"-p", (char*)"", (char*)"x", nullptr};
    EXPECT_THROW(ModuleIO::parse_command_line(4, argv), std::runtime_error);
}

TEST(ParseCommandLineTest, DigitLeadingNameThrows)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"-p", (char*)"1abc", (char*)"x", nullptr};
    EXPECT_THROW(ModuleIO::parse_command_line(4, argv), std::runtime_error);
}

TEST(ParseCommandLineTest, UnknownOptionThrows)
{
    char* argv[] = {(char*)"abacus", (char*)"--bogus", nullptr};
    EXPECT_THROW(ModuleIO::parse_command_line(2, argv), std::runtime_error);
}

// The exact-argc contract: with argc == 9 the last pair is truncated
// and must throw (documents that argc counts tokens, not nullptr).
TEST(ParseCommandLineTest, TruncatedArgcThrows)
{
    char* argv[] = {(char*)"abacus",
                    (char*)"-p", (char*)"a", (char*)"1",
                    (char*)"-p", (char*)"a", (char*)"2",
                    (char*)"-p", (char*)"b", (char*)"x",
                    nullptr};
    EXPECT_THROW(ModuleIO::parse_command_line(9, argv), std::runtime_error);
}


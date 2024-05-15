#include "module_base/formatter.h"
#include <iostream>
#include <gtest/gtest.h>
#include <vector>

TEST(FormatterTest, FmtCoreStaticFormat) {
    // const char*
    std::string result = FmtCore::format("Hello, %s!", "world");
    // remove the last '\0' character
    EXPECT_EQ(result, "Hello, world!");
    // std::string
    result = FmtCore::format("Hello, %s!", std::string("world"));
    EXPECT_EQ(result, "Hello, world!");
    // int
    result = FmtCore::format("Hello, %d!", 123);
    EXPECT_EQ(result, "Hello, 123!");
    // float
    result = FmtCore::format("Hello, %f!", 123.456);
    EXPECT_EQ(result, "Hello, 123.456000!");
    // char
    result = FmtCore::format("Hello, %c!", 'a');
    EXPECT_EQ(result, "Hello, a!");
    // invalid format
    result = FmtCore::format("Hello, %z!", "world");
    EXPECT_EQ(result, "Hello, %!");
    // varadic template case
    result = FmtCore::format("Hello, %s, %d, %f, %c!", "world", 123, 123.456, 'a');
    EXPECT_EQ(result, "Hello, world, 123, 123.456000, a!");
}

TEST(FormatterTest, FmtCoreDynamic)
{
    FmtCore fmt("Hello, %s!");
    EXPECT_EQ(fmt.fmt(), "Hello, %s!");
    std::string result = fmt.format(std::string("world"));
    EXPECT_EQ(result, "Hello, world!");

    fmt.reset("Hello, %d!");
    EXPECT_EQ(fmt.fmt(), "Hello, %d!");
    result = fmt.format(123);
    EXPECT_EQ(result, "Hello, 123!");

    fmt.reset("Hello, %f!");
    EXPECT_EQ(fmt.fmt(), "Hello, %f!");
    result = fmt.format(123.456);
    EXPECT_EQ(result, "Hello, 123.456000!");

    fmt.reset("Hello, %c!");
    EXPECT_EQ(fmt.fmt(), "Hello, %c!");
    result = fmt.format('a');
    EXPECT_EQ(result, "Hello, a!");

    // varadic template case
    fmt.reset("Hello, %s, %d, %f, %c!");
    EXPECT_EQ(fmt.fmt(), "Hello, %s, %d, %f, %c!");
    result = fmt.format(std::string("world"), 123, 123.456, 'a');
    EXPECT_EQ(result, "Hello, world, 123, 123.456000, a!");
}

TEST(FormatterTest, FmtPyStrFuncSplit)
{
    std::string fmt = "Hello, %s, %d, %f, %c!";
    // default delimiter, whitespace
    std::vector<std::string> result = FmtCore::split(fmt);
    std::vector<std::string> ref = {"Hello,", "%s,", "%d,", "%f,", "%c!"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
    fmt = "Hello, %s, %d, %f, %c";
    // other delimiter
    result = FmtCore::split(fmt, "%");
    ref = {"Hello, ", "s, ", "d, ", "f, ", "c"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
    // really string case, multiple chars
    result = FmtCore::split(fmt, ", %");
    ref = {"Hello", "s", "d", "f", "c"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
    // no such delimiter
    result = FmtCore::split(fmt, "z");
    ref = {"Hello, %s, %d, %f, %c"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
    // multiple delimiters exist
    fmt = "Hello,       %s,  %d,    %f,  %c!";
    result = FmtCore::split(fmt);
    ref = {"Hello,", "%s,", "%d,", "%f,", "%c!"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
    result = FmtCore::split(fmt, " ");
    ref = {"Hello,", "", "", "", "", "", "", "%s,", "", "%d,", "", "", "", "%f,", "", "%c!"};
    for(int i = 0; i < result.size(); i++)
    {
        EXPECT_EQ(result[i], ref[i]);
    }
}

TEST(FormatterTest, FmtPyStrFuncStartswith)
{
    const std::string fmt = "Hello, %s, %d, %f, %c!";
    EXPECT_TRUE(FmtCore::startswith(fmt, "Hello"));
    EXPECT_FALSE(FmtCore::startswith(fmt, "world"));
}

TEST(FormatterTest, FmtPyStrFuncEndswith)
{
    const std::string fmt = "Hello, %s, %d, %f, %c!";
    EXPECT_TRUE(FmtCore::endswith(fmt, "!"));
    EXPECT_FALSE(FmtCore::endswith(fmt, "world"));
}

TEST(FormatterTest, FmtPyStrFuncStrip)
{
    std::string fmt = "  Hello, %s, %d, %f, %c!  ";
    std::string result = FmtCore::strip(fmt);
    std::string ref = "Hello, %s, %d, %f, %c!";
    EXPECT_EQ(result, ref);
    fmt = "  Hello, %s, %d, %f, %c!  ";
    result = FmtCore::strip(fmt, " ");
    ref = "Hello, %s, %d, %f, %c!";
    EXPECT_EQ(result, ref);
    fmt = "";
    result = FmtCore::strip(fmt);
    ref = "";
    EXPECT_EQ(result, ref);
    fmt = "   ";
    result = FmtCore::strip(fmt);
    ref = "";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtPyStrFuncCenter)
{
    std::string fmt = "Hello, %s, %d, %f, %c!";
    std::string result = FmtCore::center(fmt, 30);
    std::string ref = "    Hello, %s, %d, %f, %c!    ";
    EXPECT_EQ(result, ref);
    result = FmtCore::center(fmt, 30, '*');
    ref = "****Hello, %s, %d, %f, %c!****";
    fmt = "Hello, %s, %d, %f, %c"; // length 21
    result = FmtCore::center(fmt, 30, '*');
    ref = "****Hello, %s, %d, %f, %c*****";
    EXPECT_EQ(result, ref);
    fmt = "";
    result = FmtCore::center(fmt, 30, '*');
    ref = "******************************";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtPyStrFuncReplace)
{
    const std::string fmt = "Hello, %s, %d, %f, %c!";
    std::string result = FmtCore::replace(fmt, "%s", "world");
    std::string ref = "Hello, world, %d, %f, %c!";
    EXPECT_EQ(result, ref);
    result = FmtCore::replace(fmt, "%d", "world");
    ref = "Hello, %s, world, %f, %c!";
    EXPECT_EQ(result, ref);
    result = FmtCore::replace(fmt, "%f", "world");
    ref = "Hello, %s, %d, world, %c!";
    EXPECT_EQ(result, ref);
    result = FmtCore::replace(fmt, "%c", "world");
    ref = "Hello, %s, %d, %f, world!";
    EXPECT_EQ(result, ref);
    result = FmtCore::replace(fmt, "%z", "world");
    ref = "Hello, %s, %d, %f, %c!";
    EXPECT_EQ(result, ref);
    result = FmtCore::replace(fmt, "%", "world");
    ref = "Hello, worlds, worldd, worldf, worldc!";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtPyStrFuncJoin)
{
    const std::vector<std::string> strs = {"Hello", "world", "!"};
    std::string result = FmtCore::join("", strs);
    std::string ref = "Helloworld!";
    EXPECT_EQ(result, ref);
    result = FmtCore::join(" ", strs);
    ref = "Hello world !";
    EXPECT_EQ(result, ref);
    result = FmtCore::join("__", strs);
    ref = "Hello__world__!";
}

TEST(FormatterTest, FmtTableDefaultArgs)
{
    const std::vector<std::string> titles = {"title1", "t i t l e 2", "t-i-t-l-e-3"};
    const std::vector<std::string> fmts = {"%s", "%d", "%f"};
    FmtTable table(titles, 5, fmts);
    const std::vector<std::string> col1 = {"row1", "row2", "row3", "row4", "row5"};
    const std::vector<int> col2 = {1, 2, 3, 4, 5};
    const std::vector<float> col3 = {1.1, 2.2, 3.3, 4.4, 5.5};
    table << col1 << col2 << col3;
    const std::string result = table.str();
    std::cout << result << std::endl;
    std::string ref = "";
    ref += "--------------------------------\n";
    ref += " title1 t i t l e 2 t-i-t-l-e-3 \n";
    ref += "--------------------------------\n";
    ref += "   row1           1    1.100000 \n";
    ref += "   row2           2    2.200000 \n";
    ref += "   row3           3    3.300000 \n";
    ref += "   row4           4    4.400000 \n";
    ref += "   row5           5    5.500000 \n";
    ref += "--------------------------------\n";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtTableCustomArgsAlign)
{
    // shared data
    std::vector<std::string> titles = {"title1", "t i t l e 2", "t-i-t-l-e-3"};
    std::vector<std::string> fmts = {"%s", "%d", "%f"};
    std::vector<std::string> col1 = {"row1", "row2", "row3", "row4", "row5"};
    std::vector<int> col2 = {1, 2, 3, 4, 5};
    std::vector<float> col3 = {1.1, 2.2, 3.3, 4.4, 5.5};
    // align: l and l
    FmtTable table(titles, 5, fmts, {FmtTable::Align::LEFT, FmtTable::Align::LEFT});
    table << col1 << col2 << col3;
    std::string result = table.str();
    std::cout << result << std::endl;
    std::string ref = "";
    ref += "--------------------------------\n";
    ref += " title1 t i t l e 2 t-i-t-l-e-3 \n";
    ref += "--------------------------------\n";
    ref += " row1   1           1.100000    \n";
    ref += " row2   2           2.200000    \n";
    ref += " row3   3           3.300000    \n";
    ref += " row4   4           4.400000    \n";
    ref += " row5   5           5.500000    \n";
    ref += "--------------------------------\n";
    EXPECT_EQ(result, ref);

    // align: r and r
    FmtTable table2(titles, 5, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::RIGHT});
    table2 << col1 << col2 << col3;
    result = table2.str();
    std::cout << result << std::endl;
    ref = "";
    ref += "--------------------------------\n";
    ref += " title1 t i t l e 2 t-i-t-l-e-3 \n";
    ref += "--------------------------------\n";
    ref += "   row1           1    1.100000 \n";
    ref += "   row2           2    2.200000 \n";
    ref += "   row3           3    3.300000 \n";
    ref += "   row4           4    4.400000 \n";
    ref += "   row5           5    5.500000 \n";
    ref += "--------------------------------\n";
    EXPECT_EQ(result, ref);

    // align: l and r
    FmtTable table3(titles, 5, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
    table3 << col1 << col2 << col3;
    result = table3.str();
    std::cout << result << std::endl;
    ref = "";
    ref += "--------------------------------\n";
    ref += " title1 t i t l e 2 t-i-t-l-e-3 \n";
    ref += "--------------------------------\n";
    ref += "   row1           1    1.100000 \n";
    ref += "   row2           2    2.200000 \n";
    ref += "   row3           3    3.300000 \n";
    ref += "   row4           4    4.400000 \n";
    ref += "   row5           5    5.500000 \n";
    ref += "--------------------------------\n";
    EXPECT_EQ(result, ref);

    // align: r and l
    FmtTable table4(titles, 5, fmts, {FmtTable::Align::LEFT, FmtTable::Align::RIGHT});
    table4 << col1 << col2 << col3;
    result = table4.str();
    std::cout << result << std::endl;
    ref = "";
    ref += "--------------------------------\n";
    ref += " title1 t i t l e 2 t-i-t-l-e-3 \n";
    ref += "--------------------------------\n";
    ref += " row1   1           1.100000    \n";
    ref += " row2   2           2.200000    \n";
    ref += " row3   3           3.300000    \n";
    ref += " row4   4           4.400000    \n";
    ref += " row5   5           5.500000    \n";
    ref += "--------------------------------\n";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtTableCustomArgsAlignFrame)
{
    // shared data
    std::vector<std::string> titles = {"title1", "t i t l e 2", "t-i-t-l-e-3"};
    std::vector<std::string> fmts = {"%s", "%d", "%f"};
    std::vector<std::string> col1 = {"row1", "row2", "row3", "row4", "row5"};
    std::vector<int> col2 = {1, 2, 3, 4, 5};
    std::vector<float> col3 = {1.1, 2.2, 3.3, 4.4, 5.5};

    FmtTable table1(titles, 5, fmts, {FmtTable::Align::LEFT, FmtTable::Align::LEFT}, {'+', '?', '*', '.', '^'});
    table1 << col1 << col2 << col3;
    std::string result = table1.str();
    std::cout << result << std::endl;
    std::string ref = "";
    ref += "++++++++++++++++++++++++++++++++\n";
    ref += ".title1 t i t l e 2 t-i-t-l-e-3^\n";
    ref += "????????????????????????????????\n";
    ref += ".row1   1           1.100000   ^\n";
    ref += ".row2   2           2.200000   ^\n";
    ref += ".row3   3           3.300000   ^\n";
    ref += ".row4   4           4.400000   ^\n";
    ref += ".row5   5           5.500000   ^\n";
    ref += "********************************\n";
    EXPECT_EQ(result, ref);
}

TEST(FormatterTest, FmtTableCustomArgsAlignFrameDelim)
{
    // shared data
    std::vector<std::string> titles = {"title1", "t i t l e 2", "t-i-t-l-e-3"};
    std::vector<std::string> fmts = {"%s", "%d", "%f"};
    std::vector<std::string> col1 = {"row1", "row2", "row3", "row4", "row5"};
    std::vector<int> col2 = {1, 2, 3, 4, 5};
    std::vector<float> col3 = {1.1, 2.2, 3.3, 4.4, 5.5};
    FmtTable table1(titles, 5, fmts, 
                    {FmtTable::Align::LEFT, FmtTable::Align::LEFT}, 
                    {'=', '/', '&', '#', '%'},
                    {'"', ']'});
    table1 << col1 << col2 << col3;
    std::string result = table1.str();
    std::cout << result << std::endl;
    std::string ref = "";
    ref += "================================\n";
    ref += "#title1]t i t l e 2]t-i-t-l-e-3%\n";
    ref += "////////////////////////////////\n";
    ref += "#row1  ]1          ]1.100000   %\n";
    ref += "#row2  ]2          ]2.200000   %\n";
    ref += "#row3  ]3          ]3.300000   %\n";
    ref += "#row4  ]4          ]4.400000   %\n";
    ref += "#row5  ]5          ]5.500000   %\n";
    ref += "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    EXPECT_EQ(result, ref);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
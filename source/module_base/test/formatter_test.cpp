#include "module_base/formatter.h"
#include <iostream>
#include <gtest/gtest.h>

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

TEST(FormatterTest, FmtTableDefaultArgs)
{
    std::vector<std::string> titles = {"title1", "t i t l e 2", "t-i-t-l-e-3"};
    std::vector<std::string> fmts = {"%s", "%d", "%f"};
    FmtTable table(titles, 5, fmts);
    std::vector<std::string> col1 = {"row1", "row2", "row3", "row4", "row5"};
    std::vector<int> col2 = {1, 2, 3, 4, 5};
    std::vector<float> col3 = {1.1, 2.2, 3.3, 4.4, 5.5};
    table << col1 << col2 << col3;
    std::string result = table.str();
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
    FmtTable table(titles, 5, fmts, {'l', 'l'});
    table << col1 << col2 << col3;
    std::string result = table.str();
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
    FmtTable table2(titles, 5, fmts, {'r', 'r'});
    table2 << col1 << col2 << col3;
    result = table2.str();
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
    FmtTable table3(titles, 5, fmts, {'r', 'l'});
    table3 << col1 << col2 << col3;
    result = table3.str();
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
    FmtTable table4(titles, 5, fmts, {'l', 'r'});
    table4 << col1 << col2 << col3;
    result = table4.str();
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

    FmtTable table1(titles, 5, fmts, {'l', 'l'}, {'+', '?', '*', '.', '^'});
    table1 << col1 << col2 << col3;
    std::string result = table1.str();
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
                    {'l', 'l'}, 
                    {'=', '/', '&', '#', '%'},
                    {'"', ']'});
    table1 << col1 << col2 << col3;
    std::string result = table1.str();
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
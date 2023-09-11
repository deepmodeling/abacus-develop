#include "../formatter_table.h"
#include <gtest/gtest.h>

/************************************************
 *  unit test of class Table
 ***********************************************/

/**
 * - Tested Functions:
 *   - Table()
 *     - default constructor, default values see formatter.h
 *   - set_mode(int mode)
 *     - setter of mode
 *   - set_col_delimiter(char col_delimiter)
 *     - setter of col_delimiter
 *   - set_frame_switches(std::vector<int> frame_switches)
 *     - setter of frame_switches
 *   - set_frame_delimiters(std::vector<char> frame_delimiters)
 *     - setter of frame_delimiters
 *   - set_frame_mid_switch(bool frame_mid_switch)
 *     - setter of frame_mid_switch
 *   - set_frame_mid_delimiter(char frame_mid_delimiter)
 *     - setter of frame_mid_delimiter
 *   - enable_flexible_width()
 *     - enable the flexible width mode, which adjusts the column width according to the data
 *   - disable_flexible_width()
 *     - disable the flexible width mode, which sets the column width to the maximum width of the data
 *   - enable_up_frame()
 *     - enable the upper frame output
 *   - disable_up_frame()
 *     - disable the upper frame output
 *   - enable_down_frame()
 *     - enable the lower frame output
 *   - disable_down_frame()
 *     - disable the lower frame output
 *   - enable_left_frame()
 *     - enable the left frame output
 *   - disable_left_frame()
 *     - disable the left frame output
 *   - enable_right_frame()
 *     - enable the right frame output
 *   - disable_right_frame()
 *     - disable the right frame output
 *   - enable_mid_frame()
 *      - enable the frame between title and data output
 *   - disable_mid_frame()
 *     - disable the frame between title and data output
 *   - add_col(std::string title, std::vector<std::string> data)
 *     - add a column to the table
 *   - reset()
 *     - reset the table
 *   - adjust_col_width()
 *     - adjust the column width according to the data
 *   - print_table()
 *     - print the table
 */

class TableTest : public testing::Test
{
};

TEST_F(TableTest, DefaultConstructor) {
    formatter::Table table;
    EXPECT_EQ(table.get_mode(), 0);
    EXPECT_EQ(table.get_ncol(), 0);
    EXPECT_EQ(table.get_col_delimiter(), ' ');
    std::vector<int> frame_switches = table.get_frame_switches();
    EXPECT_EQ(frame_switches.size(), 4);
    EXPECT_EQ(frame_switches[0], 1);
    EXPECT_EQ(frame_switches[1], 1);
    EXPECT_EQ(frame_switches[2], 0);
    EXPECT_EQ(frame_switches[3], 0);
    std::vector<char> frame_delimiters = table.get_frame_delimiters();
    EXPECT_EQ(frame_delimiters.size(), 4);
    EXPECT_EQ(frame_delimiters[0], '-');
    EXPECT_EQ(frame_delimiters[1], '-');
    EXPECT_EQ(frame_delimiters[2], '|');
    EXPECT_EQ(frame_delimiters[3], '|');
    EXPECT_EQ(table.get_frame_mid_switch(), 1);
    EXPECT_EQ(table.get_frame_mid_delimiter(), '-');
    EXPECT_TRUE(table.get_flexible_width());
    EXPECT_EQ(table.get_col_max_width(), 0);
    EXPECT_EQ(table.get_total_width(), 0);
}

TEST_F(TableTest, Setter) {
    formatter::Table table;
    table.set_mode(1);
    EXPECT_EQ(table.get_mode(), 1);
    table.set_col_delimiter(',');
    EXPECT_EQ(table.get_col_delimiter(), ',');
    table.set_frame_switches(std::vector<int>({1, 0, 1, 0}));
    std::vector<int> frame_switches = table.get_frame_switches();
    EXPECT_EQ(frame_switches.size(), 4);
    EXPECT_EQ(frame_switches[0], 1);
    EXPECT_EQ(frame_switches[1], 0);
    EXPECT_EQ(frame_switches[2], 1);
    EXPECT_EQ(frame_switches[3], 0);
    table.set_frame_delimiters(std::vector<char>{'-', '-', '|', '|'});
    std::vector<char> frame_delimiters = table.get_frame_delimiters();
    EXPECT_EQ(frame_delimiters.size(), 4);
    EXPECT_EQ(frame_delimiters[0], '-');
    EXPECT_EQ(frame_delimiters[1], '-');
    EXPECT_EQ(frame_delimiters[2], '|');
    EXPECT_EQ(frame_delimiters[3], '|');
    table.set_frame_mid_switch(0);
    EXPECT_EQ(table.get_frame_mid_switch(), 0);
    table.set_frame_mid_delimiter('+');
    EXPECT_EQ(table.get_frame_mid_delimiter(), '+');

    table.enable_flexible_width();
    EXPECT_TRUE(table.get_flexible_width());
    table.disable_flexible_width();
    EXPECT_FALSE(table.get_flexible_width());
    table.enable_up_frame();
    EXPECT_EQ(table.get_frame_switches()[0], 1);
    table.disable_up_frame();
    EXPECT_EQ(table.get_frame_switches()[0], 0);
    table.enable_down_frame();
    EXPECT_EQ(table.get_frame_switches()[1], 1);
    table.disable_down_frame();
    EXPECT_EQ(table.get_frame_switches()[1], 0);
    table.enable_left_frame();
    EXPECT_EQ(table.get_frame_switches()[2], 1);
    table.disable_left_frame();
    EXPECT_EQ(table.get_frame_switches()[2], 0);
    table.enable_right_frame();
    EXPECT_EQ(table.get_frame_switches()[3], 1);
    table.disable_right_frame();
    EXPECT_EQ(table.get_frame_switches()[3], 0);
    table.enable_mid_frame();
    EXPECT_EQ(table.get_frame_mid_switch(), 1);
    table.disable_mid_frame();
    EXPECT_EQ(table.get_frame_mid_switch(), 0);
}

TEST_F(TableTest, AddCol) {
    formatter::Table table;
    int ncol = table.get_ncol();
    std::vector<std::string> col = {"a", "b", "c"};
    std::string col_name = "col";
    table.add_col(col_name, col);
    EXPECT_EQ(table.get_ncol(), ncol + 1);
    EXPECT_EQ(table.get_titles()[ncol], col_name);
}

TEST_F(TableTest, Reset) {
    formatter::Table table;
    std::string col_name1 = "Column 1";
    std::string col_name2 = "Column 2";
    table.add_col(col_name1, (std::vector<std::string>){"Data 1"});
    table.add_col(col_name2, (std::vector<std::string>){"Data 2", "Longer Data 2"});
    table.adjust_col_width();

    table.reset();

    EXPECT_EQ(table.get_ncol(), 0);
    EXPECT_TRUE(table.get_titles().empty());
    EXPECT_EQ(table.get_col_delimiter(), ' ');
    EXPECT_EQ(table.get_frame_switches(), std::vector<int>({1, 1, 0, 0}));
    EXPECT_EQ(table.get_frame_delimiters(), std::vector<char>({'-', '-', '|', '|'}));
    EXPECT_TRUE(table.get_flexible_width());
    EXPECT_TRUE(table.get_col_widths().empty());
    EXPECT_EQ(table.get_col_max_width(), 0);
    EXPECT_EQ(table.get_total_width(), 0);
    EXPECT_TRUE(table.get_data().empty());
}

TEST_F(TableTest, AdjustColWidth) {
    formatter::Table table;

    std::string col_name1 = "Column 1";
    std::string col_name2 = "Column 2";
    table.enable_flexible_width();
    // Test with data
    table.add_col(col_name1, (std::vector<std::string>){"Data 1"});
    table.add_col(col_name2, (std::vector<std::string>){"Data 2", "Longer Data 2"});
    table.adjust_col_width();
    EXPECT_EQ(table.get_col_widths(), std::vector<int>({8, 13}));
    EXPECT_EQ(table.get_col_max_width(), 13);
    EXPECT_EQ(table.get_total_width(), 24);

    table.reset();
    // Test with no data
    table.add_col(col_name1, (std::vector<std::string>){});
    table.add_col(col_name2, (std::vector<std::string>){});
    table.adjust_col_width();
    EXPECT_EQ(table.get_col_widths(), std::vector<int>({8, 8}));
    EXPECT_EQ(table.get_col_max_width(), 8);
    EXPECT_EQ(table.get_total_width(), 19);

    table.reset();
    // Test with different flexible width setting
    table.disable_flexible_width();
    table.add_col(col_name1, (std::vector<std::string>){"Data 1"});
    table.add_col(col_name2, (std::vector<std::string>){"Data 2", "Longer Data 2"});
    table.adjust_col_width();
    EXPECT_EQ(table.get_col_widths(), std::vector<int>({13, 13}));
    EXPECT_EQ(table.get_col_max_width(), 13);
    EXPECT_EQ(table.get_total_width(), 29);
}

TEST_F(TableTest, PrintTable) {
    formatter::Table table;
    table.enable_flexible_width();
    table.set_overall_title("Table Title");
    table.add_col("Column 1", {"Data 1", "Data 2", "Data 3"});
    table.add_col("Column 2", {"Data 4", "Longer Data 5", "Data 6"});
    std::string expected_output = "Table Title\n"
                                  "------------------------\n"
                                  "Column 1 Column 2     \n"
                                  "------------------------\n"
                                  "Data 1   Data 4       \n"
                                  "Data 2   Longer Data 5\n"
                                  "Data 3   Data 6       \n"
                                  "------------------------\n";
    EXPECT_EQ(table.print_table(), expected_output);

    // Test with no data
    table.reset();
    table.disable_flexible_width();
    table.add_col("Column 1", {});
    table.add_col("Column 2", {});
    expected_output = "-------------------\n"
                      "Column 1 Column 2\n"
                      "-------------------\n"
                      "-------------------\n";
    EXPECT_EQ(table.print_table(), expected_output);

    // Test with different frame settings
    table.reset();
    table.enable_flexible_width();
    table.set_overall_title("Table Title");
    table.set_frame_switches({0, 0, 1, 1});
    table.set_frame_delimiters({'*', '*', '|', '|'});
    table.set_frame_mid_switch(false);
    table.add_col("Column 1", {"Data 1", "Data 2", "Data 3"});
    table.add_col("Column 2", {"Data 4", "Longer Data 5", "Data 6"});
    expected_output = "Table Title\n"
                      "|Column 1 Column 2     |\n"
                      "|Data 1   Data 4       |\n"
                      "|Data 2   Longer Data 5|\n"
                      "|Data 3   Data 6       |\n";
    EXPECT_EQ(table.print_table(), expected_output);
}

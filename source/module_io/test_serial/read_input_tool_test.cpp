#include "../read_input_tool.h"
#include <gtest/gtest.h>

// Test fixture for parse_expression tests
class ReadInputTool : public ::testing::Test
{
  protected:
};

TEST_F(ReadInputTool, parse_expression)
{
    // Test case for empty expressions
    {
        std::vector<std::string> expressions = {};
        std::vector<int> result;

        parse_expression(expressions, result);
        EXPECT_TRUE(result.empty());
    }
    // Test case for ""
    {
        std::vector<std::string> expressions = {""};
        std::vector<int> result;

        parse_expression(expressions, result);
        EXPECT_TRUE(result.empty());
    }
    // Test case for expressions without '*'
    {
        std::vector<std::string> expressions = {"3", "4", "10"};
        std::vector<int> expected = {3, 4, 10};
        std::vector<int> result;

        parse_expression(expressions, result);
        EXPECT_EQ(expected.size(), result.size());
        for (size_t i = 0; i < expected.size(); i++)
        {
            EXPECT_EQ(expected[i], result[i]);
        }
    }
    // Test case for expressions with one '*'
    {
        std::vector<std::string> expressions = {"3", "2*4.2", "1*7"};
        std::vector<double> expected = {3.0, 4.2, 4.2, 7.0};
        std::vector<double> result;

        
        parse_expression(expressions, result);
        EXPECT_EQ(expected.size(), result.size());
        for (size_t i = 0; i < expected.size(); i++)
        {
            EXPECT_NEAR(expected[i], result[i], 1e-5);
        }
    }
    // Test case for expressions with more than one '*'
    {
        std::vector<std::string> expressions = {"2*3*4", "1*2*3*4"};
        std::vector<double> result;

        ASSERT_THROW(parse_expression(expressions, result), std::runtime_error);
    }
}

#include "../formatter_fmt.h"
#include <gtest/gtest.h>

/************************************************
 *  unit test of class Fmt
 ***********************************************/

/**
 * - Tested Functions:
 *   - Fmt()
 *     - default constructor, default values see formatter.h
 *   - Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error)
 *     - parameterized constructor
 *   - set_width(int width)
 *     - setter of width
 *   - set_precision(int precision)
 *     - setter of precision
 *   - set_fillChar(char fillChar)
 *     - setter of fillChar
 *   - set_fixed(bool fixed)
 *     - setter of scientific notation boolean, fixed means not scientific
 *   - set_right(bool right)
 *     - setter of right alignment boolean
 *   - set_error(bool error)
 *     - setter of + sign for positive numbers boolean
 *   - format(T value)
 *     - format a value, according to width, precision, fillChar, fixed, right, error
 */

class FmtTest : public testing::Test
{
};

TEST_F(FmtTest, DefaultConstructor) {
    formatter::Fmt fmt;
    EXPECT_EQ(fmt.width(), 4);
    EXPECT_EQ(fmt.precision(), 2);
    EXPECT_EQ(fmt.fillChar(), ' ');
    EXPECT_EQ(fmt.fixed(), true);
    EXPECT_EQ(fmt.right(), true);
    EXPECT_EQ(fmt.error(), false);
}

TEST_F(FmtTest, ParameterizedConstructor) {
    formatter::Fmt fmt(10, 5, '0', false, false, true);
    EXPECT_EQ(fmt.width(), 10);
    EXPECT_EQ(fmt.precision(), 5);
    EXPECT_EQ(fmt.fillChar(), '0');
    EXPECT_EQ(fmt.fixed(), false);
    EXPECT_EQ(fmt.right(), false);
    EXPECT_EQ(fmt.error(), true);
}

TEST_F(FmtTest, Setters) {
    formatter::Fmt fmt;
    fmt.set_width(10);
    fmt.set_precision(5);
    fmt.set_fillChar('0');
    fmt.set_fixed(false);
    fmt.set_right(false);
    fmt.set_error(true);
    EXPECT_EQ(fmt.width(), 10);
    EXPECT_EQ(fmt.precision(), 5);
    EXPECT_EQ(fmt.fillChar(), '0');
    EXPECT_EQ(fmt.fixed(), false);
    EXPECT_EQ(fmt.right(), false);
    EXPECT_EQ(fmt.error(), true);
}

TEST_F(FmtTest, Format) {
    formatter::Fmt fmt;
    EXPECT_EQ(fmt.format(1), "1.00");
    EXPECT_EQ(fmt.format(1.23456789), "1.23");
    EXPECT_EQ(fmt.format(1.23456789e-10), "0.00");
    EXPECT_EQ(fmt.format(1.23456789e10), "12345678900.00");
    EXPECT_EQ(fmt.format(-1), "-1.00");
    EXPECT_EQ(fmt.format(-1.23456789), "-1.23");
    EXPECT_EQ(fmt.format(-1.23456789e-10), "-0.00");
    EXPECT_EQ(fmt.format(-1.23456789e10), "-12345678900.00");

    EXPECT_EQ(fmt.format((std::string)"hello"), "hello");
}

TEST_F(FmtTest, SpecialFormats) {
    formatter::Fmt fmt;
    fmt.set_width(8); fmt.set_precision(4); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(1), "  1.0000");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(1), " +1.0000");
    EXPECT_EQ(fmt.format(-1), " -1.0000");
}

TEST_F(FmtTest, ComplexFormat) {
    formatter::Fmt fmt;
    fmt.set_width(8); fmt.set_precision(4); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(  1.0000,  2.0000)");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "( +1.0000, +2.0000)");
    fmt.set_width(20); fmt.set_precision(10); fmt.set_fillChar(' '); fmt.set_fixed(true); fmt.set_right(true); fmt.set_error(false);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(        1.0000000000,        2.0000000000)");
    fmt.set_error(true);
    EXPECT_EQ(fmt.format(std::complex<double>(1, 2)), "(       +1.0000000000,       +2.0000000000)");
    EXPECT_EQ(fmt.format(std::complex<double>(1, -2)), "(       +1.0000000000,       -2.0000000000)");

}

TEST_F(FmtTest, ToFormat)
{
    formatter::Fmt fmt;
    int width, precision;
    bool fixed, right;
    fmt.to_format("%20.10f", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 10);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, true);
    fmt.to_format("%-20.10f", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 10);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, false);
    fmt.to_format("%20.10e", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 10);
    EXPECT_EQ(fixed, false);
    EXPECT_EQ(right, true);
    fmt.to_format("%-20.10e", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 10);
    EXPECT_EQ(fixed, false);
    EXPECT_EQ(right, false);
    fmt.to_format("%20d", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 0);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, true);
    fmt.to_format("%-20d", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 0);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, false);
    fmt.to_format("%20s", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 0);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, true);
    fmt.to_format("%-20s", width, precision, fixed, right);
    EXPECT_EQ(width, 20);
    EXPECT_EQ(precision, 0);
    EXPECT_EQ(fixed, true);
    EXPECT_EQ(right, false);
}

TEST_F(FmtTest, Format2Args)
{
    formatter::Fmt fmt;
    EXPECT_EQ(fmt.format("%20.10f", 1), "          1.0000000000");
    EXPECT_EQ(fmt.format("%20.10f", 1.23456789), "          1.2345678900");
    EXPECT_EQ(fmt.format("%20.10e", 1.23456789), "   1.2345678900e+00");
    EXPECT_EQ(fmt.format("%20.10e", 1.23456789e10), "   1.2345678900e+10");
    EXPECT_EQ(fmt.format("%20.10e", 1.23456789e-10), "   1.2345678900e-10");
    EXPECT_EQ(fmt.format("%20d", 1), "                   1");
    EXPECT_EQ(fmt.format("%20s", "hello"), "               hello");
}
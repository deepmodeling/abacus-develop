#include "../formatter_physfmt.h"
#include <gtest/gtest.h>

/************************************************
 *  unit test of class PhysicalFmt
 ***********************************************/

/**
 * - Tested Functions:
 *   - PhysicalFmt()
 *     - default constructor, default values see formatter.h
 *   - PhysicalFmt(std::string context, Fmt* p_formatter)
 *     - parameterized constructor, used to bundle with a formatter, set it according to the context and predefined rules, see function body of adjust_formatter()
 *   - adjust_formatter(bool right = true)
 *     - adjust the formatter according to the context
 *   - set_context(std::string context)
 *     - setter of context
 *   - set_p_formatter(Fmt* p_formatter)
 *     - setter of p_formatter, the pointer to Fmt object which is used to format the value
 *   - set_decorator_mode(bool decorator_mode)
 *     - setter of decorator_mode boolean
 */

class PhysFmtTest : public testing::Test
{
};

TEST_F(PhysFmtTest, DefaultConstructor) {
    formatter::PhysicalFmt physfmt;
    EXPECT_EQ(physfmt.get_context(), "none");
    EXPECT_NE(physfmt.get_p_formatter(), nullptr);
    EXPECT_FALSE(physfmt.get_decorator_mode());
}

TEST_F(PhysFmtTest, ParameterizedConstructor) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt("test", &fmt);
    EXPECT_EQ(physfmt.get_context(), "test");
    EXPECT_EQ(physfmt.get_p_formatter(), &fmt);
    EXPECT_TRUE(physfmt.get_decorator_mode());
}

TEST_F(PhysFmtTest, AdjustFormatter) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt("energy", &fmt);
    physfmt.adjust_formatter();
    EXPECT_EQ(fmt.get_width(), 20);
    EXPECT_EQ(fmt.get_precision(), 10);
    EXPECT_EQ(fmt.get_fillChar(), ' ');
    EXPECT_EQ(fmt.get_fixed(), true);
    EXPECT_EQ(fmt.get_right(), true);
    EXPECT_EQ(fmt.get_error(), false);
    physfmt.adjust_formatter(true);
    EXPECT_EQ(fmt.get_width(), 20);
    EXPECT_EQ(fmt.get_precision(), 10);
    EXPECT_EQ(fmt.get_fillChar(), ' ');
    EXPECT_EQ(fmt.get_fixed(), true);
    EXPECT_EQ(fmt.get_right(), false);
    EXPECT_EQ(fmt.get_error(), false);
}

TEST_F(PhysFmtTest, SetContext) {
    formatter::PhysicalFmt physfmt;
    physfmt.set_context("energy");
    EXPECT_EQ(physfmt.get_context(), "energy");
    EXPECT_EQ(physfmt.get_p_formatter()->get_width(), 20);
    EXPECT_EQ(physfmt.get_p_formatter()->get_precision(), 10);
    EXPECT_EQ(physfmt.get_p_formatter()->get_fillChar(), ' ');
    EXPECT_EQ(physfmt.get_p_formatter()->get_fixed(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_right(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_error(), false);
}

TEST_F(PhysFmtTest, Setpformatter) {
    formatter::Fmt fmt;
    formatter::PhysicalFmt physfmt;
    physfmt.set_p_formatter(&fmt);
    EXPECT_EQ(physfmt.get_p_formatter(), &fmt);
}

TEST_F(PhysFmtTest, SetDecoratorMode) {
    formatter::PhysicalFmt physfmt;
    physfmt.set_decorator_mode(true);
    EXPECT_TRUE(physfmt.get_decorator_mode());
}

TEST_F(PhysFmtTest, AdjustFormatterFlexible)
{
    formatter::PhysicalFmt physfmt;
    physfmt.adjust_formatter_flexible(10, 0.5); // means for double, decimal part is 10, integer part is 5
    EXPECT_EQ(physfmt.get_p_formatter()->get_width(), 16);
    EXPECT_EQ(physfmt.get_p_formatter()->get_precision(), 10);
    EXPECT_EQ(physfmt.get_p_formatter()->get_fillChar(), ' ');
    EXPECT_EQ(physfmt.get_p_formatter()->get_fixed(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_right(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_error(), false);

    physfmt.adjust_formatter_flexible(8); // means for int, width is 8
    EXPECT_EQ(physfmt.get_p_formatter()->get_width(), 8);
    EXPECT_EQ(physfmt.get_p_formatter()->get_precision(), 0);
    EXPECT_EQ(physfmt.get_p_formatter()->get_fillChar(), ' ');
    EXPECT_EQ(physfmt.get_p_formatter()->get_fixed(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_right(), true);
    EXPECT_EQ(physfmt.get_p_formatter()->get_error(), false);
}
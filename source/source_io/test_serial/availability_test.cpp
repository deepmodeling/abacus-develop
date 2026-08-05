/// Unit tests for the INPUT availability parser.
///
/// availability.{h,cpp} is dependency-free (only the C++ standard library), so
/// this test exercises the actual parser used by Input_Item directly.
#include "source_io/module_parameter/availability.h"

#include "gtest/gtest.h"

namespace ModuleIO
{

TEST(AvailabilityParser, LeafEqualityRoundTrip)
{
    const AvailabilityExpr e = parse_availability("basis_type==pw");
    EXPECT_TRUE(e.is_leaf());
    EXPECT_EQ(e.condition.param, "basis_type");
    EXPECT_EQ(e.condition.op, "==");
    EXPECT_EQ(e.condition.values, (std::vector<std::string>{"pw"}));
    EXPECT_EQ(e.to_string(), "basis_type==pw");
}

TEST(AvailabilityParser, InListRoundTrip)
{
    const AvailabilityExpr e = parse_availability("vdw_method in [d2, d3_0]");
    EXPECT_TRUE(e.is_leaf());
    EXPECT_EQ(e.condition.op, "in");
    EXPECT_EQ(e.condition.values, (std::vector<std::string>{"d2", "d3_0"}));
    EXPECT_EQ(e.to_string(), "vdw_method in [d2, d3_0]");
}

TEST(AvailabilityParser, ContainsVectorSemantics)
{
    // td_ttype is a Vector; "contains 2" is containment and must stay distinct
    // from scalar membership "in [2]".
    const AvailabilityExpr e = parse_availability("td_ttype contains 2");
    EXPECT_TRUE(e.is_leaf());
    EXPECT_EQ(e.condition.op, "contains");
    EXPECT_EQ(e.condition.values, (std::vector<std::string>{"2"}));
    EXPECT_EQ(e.to_string(), "td_ttype contains 2");
}

TEST(AvailabilityParser, AndGroup)
{
    const AvailabilityExpr e = parse_availability("basis_type==lcao and esolver_type==tddft");
    EXPECT_FALSE(e.is_leaf());
    EXPECT_EQ(e.op, "and");
    ASSERT_EQ(e.children.size(), 2u);
    EXPECT_TRUE(e.children[0].is_leaf());
    EXPECT_EQ(e.children[0].condition.param, "basis_type");
    EXPECT_TRUE(e.children[1].is_leaf());
    EXPECT_EQ(e.children[1].condition.param, "esolver_type");
    EXPECT_EQ(e.to_string(), "basis_type==lcao and esolver_type==tddft");
}

TEST(AvailabilityParser, ParenthesisedOrNesting)
{
    // "and" binds looser than "or"; the grouped "(A or B)" stays a sub-tree.
    const AvailabilityExpr e = parse_availability(
        "symmetry==1 and (dft_functional in [hse, hf] or rpa==true)");
    EXPECT_FALSE(e.is_leaf());
    EXPECT_EQ(e.op, "and");
    ASSERT_EQ(e.children.size(), 2u);
    EXPECT_TRUE(e.children[0].is_leaf());
    EXPECT_EQ(e.children[0].condition.param, "symmetry");
    const AvailabilityExpr& or_node = e.children[1];
    EXPECT_FALSE(or_node.is_leaf());
    EXPECT_EQ(or_node.op, "or");
    ASSERT_EQ(or_node.children.size(), 2u);
    EXPECT_EQ(or_node.children[0].condition.param, "dft_functional");
    EXPECT_EQ(or_node.children[1].condition.param, "rpa");
    EXPECT_EQ(e.to_string(),
              "symmetry==1 and (dft_functional in [hse, hf] or rpa==true)");
}

TEST(AvailabilityParser, EmptyIsAlwaysAvailable)
{
    const AvailabilityExpr e = parse_availability("");
    EXPECT_TRUE(e.is_leaf());
    EXPECT_TRUE(e.condition.param.empty());
}

TEST(AvailabilityParser, KeywordsNearStringEndDoNotCrash)
{
    // Exercise inputs where a keyword ("and"/"or") sits at/just before the
    // string end, which previously could index past the buffer.
    const char* edge_cases[] = {
        "and",
        " or",
        "basis_type==pw or",
        "a and",
        " and basis_type==pw",
        "basis_type==lcao and esolver_type==tddft or ",
    };
    for (const char* s : edge_cases)
    {
        const AvailabilityExpr e = parse_availability(s);
        // Parser must not crash; result is either empty or a parseable tree.
        EXPECT_NO_THROW(e.to_string());
    }
}

} // namespace ModuleIO

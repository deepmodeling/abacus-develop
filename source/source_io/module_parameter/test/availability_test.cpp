/// Unit tests for the INPUT availability parser and metadata validator.
#include "source_io/module_parameter/availability.h"
#include "source_io/module_parameter/availability_validator.h"
#include "source_io/module_parameter/input_item.h"

#include "gtest/gtest.h"

#include <map>
#include <stdexcept>

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

TEST(AvailabilityParser, SlashRemainsPartOfSingleEqualityValue)
{
    const AvailabilityExpr e = parse_availability("basis_type==pw/lcao");
    ASSERT_EQ(e.condition.values.size(), 1u);
    EXPECT_EQ(e.condition.values[0], "pw/lcao");
    EXPECT_EQ(e.to_string(), "basis_type==pw/lcao");
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
    EXPECT_EQ(e.to_string(), "basis_type==lcao and esolver_type==tddft");
}

TEST(AvailabilityParser, ParenthesisedOrNesting)
{
    const AvailabilityExpr e = parse_availability(
        "symmetry==1 and (dft_functional in [hse, hf] or rpa==true)");
    EXPECT_FALSE(e.is_leaf());
    EXPECT_EQ(e.op, "and");
    ASSERT_EQ(e.children.size(), 2u);
    EXPECT_EQ(e.to_string(),
              "symmetry==1 and (dft_functional in [hse, hf] or rpa==true)");
}

TEST(AvailabilityParser, EmptyIsAlwaysAvailable)
{
    const AvailabilityExpr e = parse_availability("");
    EXPECT_TRUE(e.is_leaf());
    EXPECT_TRUE(e.condition.param.empty());
}

TEST(AvailabilityParser, InvalidExpressionsAreRejected)
{
    const char* invalid_expressions[] = {
        "Only used for plane wave basis.",
        "   ",
        "and",
        " or",
        "basis_type==pw or",
        "a and",
        " and basis_type==pw",
        "basis_type==lcao and esolver_type==tddft or ",
        "basis_type==pw==lcao",
        "basis_type==pw)",
        "(basis_type==pw",
    };
    for (const char* expression : invalid_expressions)
    {
        EXPECT_THROW(parse_availability(expression), std::invalid_argument) << expression;
    }
}

TEST(AvailabilityParser, SetterRequiresCanonicalInputAndPreservesStateOnFailure)
{
    Input_Item item("example");
    item.set_availability("basis_type==pw");
    EXPECT_EQ(item.get_availability(), "basis_type==pw");

    EXPECT_THROW(item.set_availability("basis_type == pw"), std::invalid_argument);
    EXPECT_THROW(item.set_availability("basis_type==pw and"), std::invalid_argument);
    EXPECT_EQ(item.get_availability(), "basis_type==pw");
}

TEST(AvailabilityValidator, AcceptsKnownTypedReferences)
{
    const std::map<std::string, AvailabilityValueKind> types = {
        {"basis_type", AvailabilityValueKind::String},
        {"mixing_restart", AvailabilityValueKind::Real},
        {"td_ttype", AvailabilityValueKind::IntegerVector},
        {"relax_method", AvailabilityValueKind::StringVector},
    };
    EXPECT_NO_THROW(validate_availability_expr(
        "example",
        parse_availability("basis_type==pw and mixing_restart>0"),
        types));
    EXPECT_NO_THROW(validate_availability_expr(
        "example", parse_availability("td_ttype contains 2"), types));
    EXPECT_NO_THROW(validate_availability_expr(
        "example", parse_availability("relax_method in [cg 2]"), types));
}

TEST(AvailabilityValidator, InfersDocumentedTypes)
{
    EXPECT_EQ(availability_value_kind("Vector of string"),
              AvailabilityValueKind::StringVector);
    EXPECT_EQ(availability_value_kind("Integer \\[Integer\\](optional)"),
              AvailabilityValueKind::IntegerVector);
    EXPECT_EQ(availability_value_kind("Float"), AvailabilityValueKind::Real);
    EXPECT_EQ(availability_value_kind("Integer or string"), AvailabilityValueKind::Unknown);
}

TEST(AvailabilityValidator, RejectsUnknownAndIncompatibleReferences)
{
    const std::map<std::string, AvailabilityValueKind> types = {
        {"basis_type", AvailabilityValueKind::String},
        {"mixing_restart", AvailabilityValueKind::Real},
    };
    EXPECT_THROW(validate_availability_expr(
                     "example", parse_availability("missing==1"), types),
                 std::invalid_argument);
    EXPECT_THROW(validate_availability_expr(
                     "example", parse_availability("basis_type contains pw"), types),
                 std::invalid_argument);
    EXPECT_THROW(validate_availability_expr(
                     "example", parse_availability("basis_type>0"), types),
                 std::invalid_argument);
    EXPECT_THROW(validate_availability_expr(
                     "example", parse_availability("mixing_restart>none"), types),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, AcceptsDependenciesOfDependencies)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["esolver_type"] = parse_availability("");
    expressions["basis_type"] = parse_availability("");
    expressions["out_dos"] = parse_availability("");
    expressions["method_sto"] = parse_availability("esolver_type==sdft");
    expressions["cal_cond"] = parse_availability("basis_type==pw");
    // Each referenced parameter carries its own enclosing requirement, both
    // directly (method_sto -> esolver_type==sdft) and transitively
    // (cal_cond -> basis_type==pw) along its own path.
    expressions["npart_sto"] = parse_availability(
        "esolver_type==sdft and ((method_sto==2 and out_dos==1) or (basis_type==pw and cal_cond==true))");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "npart_sto", expressions["npart_sto"], expressions));

    // A requirement present only on a different OR branch does not satisfy the
    // referencing path.
    expressions["broken"] = parse_availability(
        "(method_sto==2 and out_dos==1) or (esolver_type==sdft and cal_cond==true)");
    EXPECT_THROW(validate_availability_self_contained(
                     "broken", expressions["broken"], expressions),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, RejectsMissingTransitiveRequirement)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["esolver_type"] = parse_availability("");
    expressions["out_dos"] = parse_availability("");
    expressions["method_sto"] = parse_availability("esolver_type==sdft");
    expressions["npart_sto"] = parse_availability("method_sto==2 and out_dos==1");
    EXPECT_THROW(validate_availability_self_contained(
                     "npart_sto", expressions["npart_sto"], expressions),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, AcceptsOrPrerequisiteImpliedByPath)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["q"] = parse_availability("");
    expressions["r"] = parse_availability("");
    expressions["p"] = parse_availability("q==2 or r==3");

    expressions["factored"] = parse_availability("p==1 and (r==3 or q==2)");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "factored", expressions["factored"], expressions));

    expressions["all"] = parse_availability("q==2 and r==3");
    expressions["and_order"] = parse_availability("all==1 and r==3 and q==2");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "and_order", expressions["and_order"], expressions));

    expressions["distributed"] = parse_availability(
        "(p==1 and q==2) or (p==1 and r==3)");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "distributed", expressions["distributed"], expressions));

    expressions["basis_type"] = parse_availability("");
    expressions["calculation"] = parse_availability("");
    expressions["gamma_only"] = parse_availability("");
    expressions["out_pchg"] = parse_availability(
        "basis_type==pw or (basis_type==lcao and calculation==get_pchg)");
    expressions["if_separate_k"] = parse_availability(
        "(basis_type==pw and out_pchg!=none) or "
        "(basis_type==lcao and calculation==get_pchg and gamma_only==0)");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "if_separate_k", expressions["if_separate_k"], expressions));

    expressions["missing"] = parse_availability("p==1");
    EXPECT_THROW(validate_availability_self_contained(
                     "missing", expressions["missing"], expressions),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, RequiresExactNonEqualityPrerequisite)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["mode"] = parse_availability("");
    expressions["feature"] = parse_availability("mode in [a, b]");

    expressions["explicit"] = parse_availability(
        "feature==enabled and mode in [a, b]");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "explicit", expressions["explicit"], expressions));

    expressions["inferred"] = parse_availability("feature==enabled and mode==a");
    EXPECT_THROW(validate_availability_self_contained(
                     "inferred", expressions["inferred"], expressions),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, ChecksPrerequisitesTransitivelyThroughOr)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["root"] = parse_availability("");
    expressions["q"] = parse_availability("root==ready");
    expressions["r"] = parse_availability("");
    expressions["p"] = parse_availability("q==2 or r==3");

    expressions["complete"] = parse_availability(
        "root==ready and p==1 and (q==2 or r==3)");
    EXPECT_NO_THROW(validate_availability_self_contained(
        "complete", expressions["complete"], expressions));

    expressions["missing"] = parse_availability("p==1 and (q==2 or r==3)");
    EXPECT_THROW(validate_availability_self_contained(
                     "missing", expressions["missing"], expressions),
                 std::invalid_argument);
}

TEST(AvailabilityValidator, RejectsSelfReference)
{
    std::map<std::string, AvailabilityExpr> expressions;
    expressions["self"] = parse_availability("self==true");
    EXPECT_THROW(validate_availability_self_contained(
                     "self", expressions["self"], expressions),
                 std::invalid_argument);

    expressions["duplicate_self"] = parse_availability(
        "duplicate_self==true and duplicate_self==true");
    EXPECT_THROW(validate_availability_self_contained(
                     "duplicate_self", expressions["duplicate_self"], expressions),
                 std::invalid_argument);
}

} // namespace ModuleIO

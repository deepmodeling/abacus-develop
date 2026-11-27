#include "../input_domain.h"
#include "../input_parameter.h"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(InputDomain, UndeclaredDomainKeepsLegacyValidation)
{
    const ModuleIO::InputDomain domain;
    const ModuleIO::InputValues input;
    EXPECT_TRUE(domain.empty());
    EXPECT_TRUE(domain.format_choices().empty());
    EXPECT_TRUE(domain.validate_explicit("legacy", {"any", "tokens"}).empty());
    EXPECT_TRUE(domain.validate_resolved("legacy", input).empty());
}

TEST(InputDomain, RejectsInvalidSchema)
{
    ModuleIO::InputDomain domain;
    EXPECT_THROW(domain.set_choices({}), std::invalid_argument);
    EXPECT_THROW(domain.set_choices({"1", "2", "1"}), std::invalid_argument);
    domain.set_choices({"1"});
    EXPECT_THROW(domain.validate_schema("unbound", {}), std::invalid_argument);
}

TEST(InputDomain, ParsesCompleteScalarIntegers)
{
    ModuleIO::InputDomain domain;
    domain.bind(ModuleIO::AvailabilityValueKind::Integer, [](const Parameter&) { return ModuleIO::InputValue(1); });
    domain.set_choices({"1", "2", "4"});
    for (const std::string token : {"1", "2", "4", "+2", "02"})
    {
        EXPECT_TRUE(domain.validate_explicit("nspin", {token}).empty()) << token;
    }
    for (const std::vector<std::string> tokens : {
             std::vector<std::string>{}, {""}, {"3"}, {"-1"}, {"2abc"}, {"2.0"}, {"2e0"},
             {"999999999999999999999999"}, {"2", "4"}})
    {
        const std::string error = domain.validate_explicit("nspin", tokens);
        EXPECT_NE(error.find("Invalid INPUT parameter 'nspin'"), std::string::npos);
        EXPECT_NE(error.find("Expected one Integer from {1, 2, 4}"), std::string::npos);
        EXPECT_NE(error.find("input-main.html#nspin"), std::string::npos);
    }
}

TEST(InputDomain, ValuesAreReadFromTypedBindingWithoutMutation)
{
    Input_para input;
    ModuleIO::InputDomain domain;
    domain.bind(ModuleIO::AvailabilityValueKind::Integer, [&input](const Parameter&) { return ModuleIO::InputValue(input.gga_grad); });
    domain.set_choices({"-1", "0", "2"});
    const ModuleIO::InputDomain copy = domain;
    for (const int value : {-1, 0, 2})
    {
        input.gga_grad = value;
        EXPECT_TRUE(copy.validate_resolved("gga_grad", {{"gga_grad", ModuleIO::InputValue(input.gga_grad)}}).empty());
        EXPECT_EQ(input.gga_grad, value);
    }
    input.gga_grad = 1;
    const std::string error = copy.validate_resolved("gga_grad", {{"gga_grad", ModuleIO::InputValue(input.gga_grad)}});
    EXPECT_NE(error.find("resolved value '1'"), std::string::npos);
    EXPECT_NE(error.find("input-main.html#gga-grad"), std::string::npos);
    EXPECT_EQ(input.gga_grad, 1);
    EXPECT_EQ(copy.choices(), (std::vector<std::string>{"-1", "0", "2"}));
    EXPECT_EQ(copy.format_choices(), "-1, 0, 2");
}

TEST(InputDomain, NumericBoundsAndNonfiniteValues)
{
    ModuleIO::InputDomain domain;
    domain.bind(ModuleIO::AvailabilityValueKind::Real, [](const Parameter&) { return ModuleIO::InputValue(0.5); });
    domain.set_rule("value>=0 and value<1");
    for (const std::string token : {"0", "0.5", "5e-1", "-0"})
    {
        EXPECT_TRUE(domain.validate_explicit("weight", {token}).empty());
    }
    for (const std::string token : {"-0.1", "1", "nan", "inf", "-inf", "1e999", "0.5garbage", "null"})
    {
        EXPECT_FALSE(domain.validate_explicit("weight", {token}).empty()) << token;
    }
}

TEST(InputDomain, BooleanAliasesAndExactCardinality)
{
    ModuleIO::InputDomain domain;
    domain.bind(ModuleIO::AvailabilityValueKind::Boolean, [](const Parameter&) { return ModuleIO::InputValue(false); });
    domain.scalar();
    for (const std::string token : {"true", "FALSE", "0", "1", "Y", "n", "on", "off", ".true.", ".FALSE."})
    {
        EXPECT_TRUE(domain.validate_explicit("flag", {token}).empty());
    }
    EXPECT_FALSE(domain.validate_explicit("flag", {"2"}).empty());
    EXPECT_FALSE(domain.validate_explicit("flag", {"true", "false"}).empty());
}

TEST(InputDomain, SharedConditionsAndDependencyDiagnostics)
{
    const ModuleIO::InputValues values = {
        {"nspin", ModuleIO::InputValue(2)}, {"noncolin", ModuleIO::InputValue(true)},
        {"basis_type", ModuleIO::InputValue(std::string("pw"))},
        {"directions", ModuleIO::InputValue(std::vector<int>{0, 2})},
        {"threshold", ModuleIO::InputValue(0.5)}};
    for (const std::string predicate : {"nspin in [1, 2, 4]", "threshold>0 and threshold<=0.5",
                                       "basis_type!=lcao", "directions contains 2",
                                       "nspin==4 or noncolin==true"})
    {
        const auto expression = ModuleIO::parse_availability(predicate);
        EXPECT_TRUE(ModuleIO::evaluate_input_condition(expression, values)) << predicate;
    }
    ModuleIO::InputDomain domain;
    domain.add_requirement("noncolin==true", "nspin==4");
    const auto error = domain.validate_requirements("nspin", values);
    EXPECT_NE(error.find("noncolin=true"), std::string::npos);
    EXPECT_NE(error.find("nspin=2"), std::string::npos);
    EXPECT_NE(error.find("Requires nspin==4"), std::string::npos);
    EXPECT_NE(error.find("input-main.html#nspin"), std::string::npos);
    EXPECT_THROW(domain.validate_schema("nspin", {}), std::invalid_argument);
}

TEST(InputDomain, InternalUnsetValueIsNotAPublicChoice)
{
    ModuleIO::InputDomain domain;
    domain.bind(ModuleIO::AvailabilityValueKind::Integer, [](const Parameter&) { return ModuleIO::InputValue(1); });
    domain.set_choices({"1", "2"});
    EXPECT_FALSE(domain.validate_explicit("criterion", {"-1"}).empty());
    const ModuleIO::InputValues unset = {{"criterion", ModuleIO::InputValue(-1)}};
    const ModuleIO::InputValues resolved = {{"criterion", ModuleIO::InputValue(2)}};
    EXPECT_FALSE(domain.validate_resolved("criterion", unset).empty());
    EXPECT_TRUE(domain.validate_resolved("criterion", resolved).empty());
}

TEST(InputDomain, RejectsNonfiniteRuleLiteralsAndVectorScalarComparisons)
{
    ModuleIO::InputDomain domain;
    domain.add_requirement("threshold>nan", "flag==true");
    const std::map<std::string, ModuleIO::AvailabilityValueKind> types = {
        {"threshold", ModuleIO::AvailabilityValueKind::Real},
        {"flag", ModuleIO::AvailabilityValueKind::Boolean},
        {"axes", ModuleIO::AvailabilityValueKind::IntegerVector}};
    EXPECT_THROW(domain.validate_schema("bad_threshold", types), std::invalid_argument);
    ModuleIO::InputDomain vector_domain;
    vector_domain.add_requirement("axes==1", "flag==true");
    EXPECT_THROW(vector_domain.validate_schema("bad_vector", types), std::invalid_argument);
}

#include "../input_domain.h"
#include "../input_parameter.h"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(InputDomain, UndeclaredDomainKeepsLegacyValidation)
{
    const ModuleIO::InputDomain domain;
    const Input_para input;
    EXPECT_TRUE(domain.empty());
    EXPECT_TRUE(domain.format_choices().empty());
    EXPECT_TRUE(domain.validate_explicit("legacy", {"any", "tokens"}).empty());
    EXPECT_TRUE(domain.validate_resolved("legacy", input).empty());
}

TEST(InputDomain, RejectsInvalidSchema)
{
    EXPECT_THROW(ModuleIO::InputDomain(nullptr, {1}), std::invalid_argument);
    EXPECT_THROW(ModuleIO::InputDomain(&Input_para::nspin, {}), std::invalid_argument);
    EXPECT_THROW(ModuleIO::InputDomain(&Input_para::nspin, {1, 2, 1}), std::invalid_argument);
}

TEST(InputDomain, ParsesCompleteScalarIntegers)
{
    const ModuleIO::InputDomain domain(&Input_para::nspin, {1, 2, 4});
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
    const ModuleIO::InputDomain domain(&Input_para::gga_grad, {-1, 0, 2});
    const ModuleIO::InputDomain copy = domain;
    Input_para input;
    for (const int value : {-1, 0, 2})
    {
        input.gga_grad = value;
        EXPECT_TRUE(copy.validate_resolved("gga_grad", input).empty());
        EXPECT_EQ(input.gga_grad, value);
    }
    input.gga_grad = 1;
    const std::string error = copy.validate_resolved("gga_grad", input);
    EXPECT_NE(error.find("resolved value '1'"), std::string::npos);
    EXPECT_NE(error.find("input-main.html#gga-grad"), std::string::npos);
    EXPECT_EQ(input.gga_grad, 1);
    EXPECT_EQ(copy.integer_choices(), (std::vector<int>{-1, 0, 2}));
    EXPECT_EQ(copy.format_choices(), "-1, 0, 2");
}

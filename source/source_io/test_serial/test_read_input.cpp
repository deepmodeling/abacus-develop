#include "source_io/module_parameter/read_input.h"
#include "source_io/input_help.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace
{

class InputDomainReaderTest : public testing::Test
{
  protected:
    std::string filename;

    void SetUp() override
    {
        filename = testing::TempDir() + "abacus_input_domain_INPUT";
    }

    void TearDown() override
    {
        std::remove(filename.c_str());
    }

    void write_input(const std::string& values)
    {
        std::ofstream out(filename);
        out << "INPUT_PARAMETERS\nbasis_type pw\nntype 1\necutwfc 20\n" << values;
        ASSERT_TRUE(out.good());
    }

    std::string reject_input(const std::string& values)
    {
        write_input(values);
        ModuleIO::ReadInput reader(0);
        reader.check_ntype_flag = false;
        Parameter input;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(reader.read_parameters(input, filename), testing::ExitedWithCode(1), "");
        return testing::internal::GetCapturedStdout();
    }
};

TEST_F(InputDomainReaderTest, KeepsExplicitValuesAndContextualDefaults)
{
    for (const int spin : {1, 2, 4})
    {
        for (const int grad : {0, 1, 2})
        {
            const std::string values = "gga_grad " + std::to_string(grad) + "\nnspin " + std::to_string(spin) + "\n";
            write_input(values);
            ModuleIO::ReadInput reader(0);
            reader.check_ntype_flag = false;
            Parameter input;
            reader.read_parameters(input, filename);
            EXPECT_EQ(input.inp.nspin, spin);
            EXPECT_EQ(input.inp.gga_grad, grad);
            const double expected_beta = spin == 1 ? 0.8 : 0.4;
            EXPECT_DOUBLE_EQ(input.inp.mixing_beta, expected_beta);
        }
    }

    write_input("");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_EQ(input.inp.nspin, 1);
    EXPECT_EQ(input.inp.gga_grad, 0);
}

TEST_F(InputDomainReaderTest, RejectsMalformedAndOutOfDomainExplicitValues)
{
    for (const std::string value : {"3", "0", "-1", "2abc", "2.0", "2 4", "999999999999999999999"})
    {
        const std::string values = "nspin " + value + "\n";
        const std::string output = reject_input(values);
        EXPECT_NE(output.find("supplied '" + value + "'"), std::string::npos);
        EXPECT_NE(output.find("Expected one Integer from {1, 2, 4}"), std::string::npos);
        EXPECT_NE(output.find("input-main.html#nspin"), std::string::npos);
        EXPECT_NE(output.find(filename), std::string::npos);
    }
    const std::string output = reject_input("gga_grad 3\n");
    EXPECT_NE(output.find("Expected one Integer from {0, 1, 2}"), std::string::npos);
    EXPECT_NE(output.find("input-main.html#gga-grad"), std::string::npos);
}

TEST_F(InputDomainReaderTest, ReportsIndependentErrorsTogether)
{
    const std::string output = reject_input("nspin 3\ngga_grad -1\n");
    EXPECT_NE(output.find("Invalid INPUT parameter 'nspin'"), std::string::npos);
    EXPECT_NE(output.find("Invalid INPUT parameter 'gga_grad'"), std::string::npos);
}

TEST_F(InputDomainReaderTest, RejectsMissingValuesIncludingLastLineWithoutNewline)
{
    for (const std::string values : {"nspin\n", "nspin", "gga_grad # missing value\n"})
    {
        const std::string output = reject_input(values);
        EXPECT_NE(output.find("supplied ''"), std::string::npos);
    }
    const std::string output = reject_input("nspin\nnspin 2\n");
    EXPECT_NE(output.find("has been read twice"), std::string::npos);
    const std::string legacy_output = reject_input("calculation");
    EXPECT_NE(legacy_output.find("'calculation': missing value (String)"), std::string::npos);
}

TEST_F(InputDomainReaderTest, RetainsExistingSpinDependencyCheck)
{
    const std::string output = reject_input("nspin 2\nnoncolin 1\n");
    EXPECT_NE(output.find("nspin must be 4 when noncolin or lspinorb is enabled"), std::string::npos);
}

TEST(InputDomainHelp, ExportsTheRuntimeRules)
{
    const auto spin = ModuleIO::ParameterHelp::get_metadata("nspin");
    EXPECT_EQ(spin.allowed_values, "1, 2, 4");
    std::ostringstream help;
    ASSERT_TRUE(ModuleIO::ParameterHelp::show_parameter_help("gga_grad", help));
    EXPECT_NE(help.str().find("Allowed values: 0, 1, 2"), std::string::npos);
    std::ostringstream yaml;
    ModuleIO::ParameterHelp::generate_yaml(yaml);
    const std::string dump = yaml.str();
    EXPECT_NE(dump.find("allowed_values: [1, 2, 4]"), std::string::npos);
    EXPECT_NE(dump.find("allowed_values: [0, 1, 2]"), std::string::npos);
}

} // namespace

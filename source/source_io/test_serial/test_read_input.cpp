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
        const std::string written = filename + ".written";
        std::remove(written.c_str());
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
}

TEST_F(InputDomainReaderTest, PreservesLegacyEmptyValues)
{
    write_input("pseudo_dir\norbital_dir\ndfpt_qfile\nout_pchg\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_TRUE(input.inp.pseudo_dir.empty());
    EXPECT_TRUE(input.inp.orbital_dir.empty());
    EXPECT_TRUE(input.inp.dfpt_qfile.empty());
    EXPECT_TRUE(input.inp.out_pchg.empty());
}

TEST_F(InputDomainReaderTest, RejectsMalformedRpaCutoffTokens)
{
    for (const std::string cut_name : {"rpa_cut_coulomb", "rpa_cut_radius", "rpa_cut_width"})
    {
        for (const std::string cut_tokens : {"\nscf_nmax 20\n", "", " invalid\n", " 1 2\n"})
        {
            const std::string cut_input = cut_name + cut_tokens;
            const std::string cut_output = reject_input(cut_input);
            const std::string cut_error = "Invalid INPUT parameter '" + cut_name + "'";
            EXPECT_NE(cut_output.find(cut_error), std::string::npos) << cut_input;
        }
    }
    for (const std::string cut_input : {"rpa_cut_coulomb log_erfc extra\n",
                                       "rpa_cut_radius inscribed extra\n",
                                       "rpa_cut_width 1\n", "rpa_cut_width 1.2501\n",
                                       "rpa_cut_width nan\n", "rpa_cut_width inf\n"})
    {
        const std::string cut_output = reject_input(cut_input);
        EXPECT_NE(cut_output.find("Invalid INPUT parameter 'rpa_cut_"), std::string::npos) << cut_input;
    }
}

TEST_F(InputDomainReaderTest, AcceptsRpaCutoffChoicesAndWidthUpperBound)
{
    for (const std::string cut_mode : {"hard", "log_erfc"})
    {
        for (const std::string cut_radius : {"legacy", "inscribed"})
        {
            const std::string cut_input = "rpa_cut_coulomb " + cut_mode + "\nrpa_cut_radius "
                + cut_radius + "\nrpa_cut_width 1.25\n";
            write_input(cut_input);
            ModuleIO::ReadInput cut_reader(0);
            cut_reader.check_ntype_flag = false;
            Parameter cut_parameters;
            cut_reader.read_parameters(cut_parameters, filename);
            EXPECT_EQ(cut_parameters.inp.rpa_cut_coulomb, cut_mode);
            EXPECT_EQ(cut_parameters.inp.rpa_cut_radius, cut_radius);
            EXPECT_DOUBLE_EQ(cut_parameters.inp.rpa_cut_width, 1.25);
        }
    }
}

TEST_F(InputDomainReaderTest, RejectsMissingOrdinaryStringValues)
{
    for (const std::string label : {"stru_file", "kpoint_file", "suffix", "dft_functional"})
    {
        for (const std::string ending : {"", "\n", " # missing value\n", " ! missing value\n"})
        {
            const std::string values = label + ending;
            const std::string output = reject_input(values);
            const std::string diagnostic = "Invalid INPUT parameter '" + label + "': supplied ''";
            const std::string manual = ModuleIO::input_manual_link(label);
            EXPECT_NE(output.find(diagnostic), std::string::npos);
            EXPECT_NE(output.find("Expected one String."), std::string::npos);
            EXPECT_NE(output.find(manual), std::string::npos);
            EXPECT_NE(output.find(filename), std::string::npos);
        }
    }
}

TEST_F(InputDomainReaderTest, RejectsExtraOrdinaryStringValues)
{
    const std::string output = reject_input("stru_file STRU extra\n");
    EXPECT_NE(output.find("supplied 'STRU extra'"), std::string::npos);
    EXPECT_NE(output.find("Expected one String."), std::string::npos);
    EXPECT_NE(output.find("input-main.html#stru-file"), std::string::npos);
}

TEST_F(InputDomainReaderTest, KeepsExplicitStringPathsAndOmittedDefaults)
{
    write_input("stru_file /tmp/structure\nkpoint_file //server/KPT # custom path\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_EQ(input.inp.stru_file, "/tmp/structure");
    EXPECT_EQ(input.inp.kpoint_file, "//server/KPT");
    EXPECT_EQ(input.inp.suffix, "ABACUS");

    write_input("");
    ModuleIO::ReadInput defaults_reader(0);
    defaults_reader.check_ntype_flag = false;
    Parameter defaults;
    defaults_reader.read_parameters(defaults, filename);
    EXPECT_EQ(defaults.inp.stru_file, "STRU");
    EXPECT_EQ(defaults.inp.kpoint_file, "KPT");
}

TEST_F(InputDomainReaderTest, KeepsTrailingCommentsAndAbsolutePaths)
{
    write_input("scf_thr 1e-7 // about iteration\n"
                "nspin 2 // magnetic calculation\n"
                "pseudo_dir /tmp//pseudo\norbital_dir //server/orbitals\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_DOUBLE_EQ(input.inp.scf_thr, 1e-7);
    EXPECT_EQ(input.inp.nspin, 2);
    EXPECT_EQ(input.inp.pseudo_dir, "/tmp//pseudo/");
    EXPECT_EQ(input.inp.orbital_dir, "//server/orbitals/");
}

TEST_F(InputDomainReaderTest, ReadsWriterOutputWithEmptyOptionalValues)
{
    write_input("nspin 2\ngga_grad 0\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    const std::string written = filename + ".written";
    reader.write_parameters(input, written);

    ModuleIO::ReadInput roundtrip_reader(0);
    roundtrip_reader.check_ntype_flag = false;
    Parameter roundtrip;
    roundtrip_reader.read_parameters(roundtrip, written);
    EXPECT_EQ(roundtrip.inp.nspin, 2);
    EXPECT_EQ(roundtrip.inp.gga_grad, 0);
    EXPECT_EQ(roundtrip.inp.pseudo_dir, input.inp.pseudo_dir);
    EXPECT_EQ(roundtrip.inp.dfpt_qfile, input.inp.dfpt_qfile);
    EXPECT_EQ(roundtrip.inp.out_pchg, input.inp.out_pchg);
}

TEST_F(InputDomainReaderTest, RejectsMissingStringChoiceAtEndOfFile)
{
    const std::string output = reject_input("calculation");
    EXPECT_NE(output.find("Invalid INPUT parameter 'calculation': supplied ''"), std::string::npos);
}

TEST_F(InputDomainReaderTest, RetainsExistingSpinDependencyCheck)
{
    const std::string output = reject_input("nspin 2\nnoncolin 1\n");
    EXPECT_NE(output.find("Requires nspin==4"), std::string::npos);
}

TEST_F(InputDomainReaderTest, RejectsPrivateSentinelsBeforeDefaultResolution)
{
    for (const std::string value : {"-1", "0", "3", "1abc"})
    {
        const std::string output = reject_input("scf_thr_type " + value + "\n");
        EXPECT_NE(output.find("supplied '" + value + "'"), std::string::npos);
        EXPECT_NE(output.find("from {1, 2}"), std::string::npos);
        EXPECT_NE(output.find("input-main.html#scf-thr-type"), std::string::npos);
    }
    for (const std::string basis : {"pw", "lcao"})
    {
        std::ofstream out(filename);
        out << "INPUT_PARAMETERS\nbasis_type " << basis << "\nntype 1\n";
        out.close();
        ModuleIO::ReadInput reader(0);
        reader.check_ntype_flag = false;
        Parameter input;
        reader.read_parameters(input, filename);
        const int expected = basis == "pw" ? 1 : 2;
        EXPECT_EQ(input.inp.scf_thr_type, expected);
    }
    write_input("scf_thr_type 2\nmixing_beta 0\nmixing_beta_mag 0.2\nmixing_gg0_mag 0.7\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_EQ(input.inp.scf_thr_type, 2);
    EXPECT_EQ(input.inp.mixing_beta, 0);
    EXPECT_EQ(input.inp.mixing_beta_mag, 0.2);
    EXPECT_EQ(input.inp.mixing_gg0_mag, 0.7);
}

TEST_F(InputDomainReaderTest, ChecksTypesRangesAndDependenciesThroughOnePath)
{
    const std::vector<std::pair<std::string, std::string>> invalid = {
        {"scf_nmax 2oops", "scf-nmax"}, {"mixing_beta nan", "mixing-beta"},
        {"mixing_beta -0.5", "mixing-beta"}, {"cal_force maybe", "cal-force"},
        {"symmetry 2", "symmetry"}, {"fft_mode 4", "fft-mode"},
        {"md_thermostat unknown", "md-thermostat"}, {"lin_maxiter 0", "lin-maxiter"},
        {"out_wfc_lcao 1", "out-wfc-lcao"}, {"dip_cor_flag true", "dip-cor-flag"}};
    for (const auto& example : invalid)
    {
        const std::string output = reject_input(example.first + "\n");
        EXPECT_NE(output.find("Invalid INPUT parameter"), std::string::npos) << example.first;
        EXPECT_NE(output.find("input-main.html#" + example.second), std::string::npos) << example.first;
    }
    // The legacy nupdown reset used to erase this conflicting explicit value.
    const std::string output = reject_input("nupdown 1\nnspin 1\n");
    EXPECT_NE(output.find("nupdown=1"), std::string::npos);
    EXPECT_NE(output.find("Requires nupdown==0"), std::string::npos);
}

TEST_F(InputDomainReaderTest, PreservesDisabledModesNegativePhysicsAndInactiveSwitches)
{
    write_input("scf_ene_thr -1\nnbspline -1\nmixing_angle -10\n"
                "md_pfirst -2\nmd_plast -1\nout_wfc_lcao 0\ncal_force false\n");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    EXPECT_EQ(input.inp.scf_ene_thr, -1);
    EXPECT_EQ(input.inp.nbspline, -1);
    EXPECT_EQ(input.inp.mixing_angle, -10);
    EXPECT_EQ(input.mdp.md_pfirst, -2);
    EXPECT_EQ(input.mdp.md_plast, -1);
    EXPECT_EQ(input.inp.out_wfc_lcao, 0);
    EXPECT_FALSE(input.inp.cal_force);
    const std::string written = filename + ".written";
    reader.write_parameters(input, written);
    ModuleIO::ReadInput rereader(0);
    rereader.check_ntype_flag = false;
    Parameter restored;
    rereader.read_parameters(restored, written);
    EXPECT_EQ(restored.inp.nspin, 1);
    EXPECT_EQ(restored.inp.mixing_beta_mag, input.inp.mixing_beta_mag);
    EXPECT_EQ(restored.mdp.md_plast, -1);
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
    EXPECT_NE(dump.find("allowed_values: [\"-1\", \"0\", \"1\"]"), std::string::npos);
    EXPECT_NE(dump.find("value_rule:"), std::string::npos);
    EXPECT_NE(dump.find("requirements:"), std::string::npos);
    const auto criterion = ModuleIO::ParameterHelp::get_metadata("scf_thr_type");
    EXPECT_EQ(criterion.allowed_values, "1, 2");
    EXPECT_NE(spin.requirements.find("nspin==4"), std::string::npos);
}

} // namespace

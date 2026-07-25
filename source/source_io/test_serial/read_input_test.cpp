#include "source_io/module_parameter/read_input.h"

#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>

// mock
namespace GlobalV
{
int NPROC = 4;
int MY_RANK = 0;
std::ofstream ofs_running;
std::ofstream ofs_warning;
} // namespace GlobalV
namespace ModuleBase
{
void TITLE(const std::string& class_name, const std::string& function_name, const bool disable)
{
}
namespace GlobalFunc
{
bool SCAN_BEGIN(std::ifstream& ifs, const std::string& TargetName, const bool restart, const bool ifwarn)
{
    return false;
}
} // namespace GlobalFunc
namespace Global_File
{
void make_dir_out(const std::string& suffix,
                  const std::string& calculation,
                  const bool& out_dir,
                  const bool& out_wfc_dir,
                  const int rank,
                  const bool& restart,
                  const bool out_alllog,
                  const std::string& global_out_dir,
                  const std::string& global_stru_dir,
                  const std::string& global_matrix_dir,
                  const std::string& global_wfc_dir,
                  const std::string& global_mlkedf_descriptor_dir,
                  const std::string& global_deepks_label_elec_dir,
                  const std::string& log_file,
                  const bool of_ml_gene_data,
                  const bool deepks_out_freq_elec)
{
}
} // namespace Global_File
} // namespace ModuleBase

/************************************************
 *  unit test of read_input_test.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Selfconsistent_Read:
 *     - read empty INPUT file and write INPUT.ref back
 *     - read INPUT.ref file again and write INPUT
 *   - Check:
 *     - check_mode = true
 */

class InputTest : public testing::Test
{
  protected:
    void TearDown() override
    {
        ModuleIO::ReadInput::check_mode = false;
    }

    bool compare_two_files(const std::string& filename1, const std::string& filename2)
    {
        std::ifstream file1(filename1.c_str());
        std::ifstream file2(filename2.c_str());
        EXPECT_TRUE(file1.is_open());
        EXPECT_TRUE(file2.is_open());

        std::string line1, line2;
        int lineNumber = 1;
        bool allpass = true;
        while (std::getline(file1, line1) && std::getline(file2, line2))
        {
            std::istringstream iss1(line1);
            std::istringstream iss2(line2);

            std::string col1_file1, col2_file1;
            std::string col1_file2, col2_file2;

            // read two columns from each file
            iss1 >> col1_file1 >> col2_file1;
            iss2 >> col1_file2 >> col2_file2;

            // compare two columns
            // compare two columns
            if (col1_file1 != col1_file2 || col2_file1 != col2_file2)
            {
                std::cout << "Mismatch found at line " << lineNumber << " in files " << filename1 << " and "
                          << filename2 << std::endl;
                std::cout << "File1: " << col1_file1 << " " << col2_file1 << std::endl;
                std::cout << "File2: " << col1_file2 << " " << col2_file2 << std::endl;
                allpass = false;
            }

            lineNumber++;
        }

        file1.close();
        file2.close();
        return allpass;
    }
};

TEST_F(InputTest, Selfconsistent_Read)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    { // PW
        std::ofstream emptyfile("empty_INPUT");
        emptyfile << "INPUT_PARAMETERS";
        emptyfile.close();
        Parameter param;
        // readinput.read_parameters(param, "./empty_INPUT");
        EXPECT_NO_THROW(readinput.read_parameters(param, "./empty_INPUT"));
        readinput.write_parameters(param, "./my_INPUT1");
        readinput.clear();
        // readinput.read_parameters(param, "./my_INPUT1");
        EXPECT_NO_THROW(readinput.read_parameters(param, "./my_INPUT1"));
        readinput.write_parameters(param, "./my_INPUT2");
        EXPECT_TRUE(compare_two_files("./my_INPUT1", "./my_INPUT2"));
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT1") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT2") == 0);
        readinput.clear();
    }
    { // LCAO
        std::ofstream emptyfile("empty_INPUT");
        emptyfile << "INPUT_PARAMETERS\n";
        emptyfile << "basis_type           lcao";
        emptyfile.close();
        Parameter param;
        // readinput.read_parameters(param, "./empty_INPUT");
        EXPECT_NO_THROW(readinput.read_parameters(param, "./empty_INPUT"));
        readinput.write_parameters(param, "./my_INPUT1");
        readinput.clear();
        // readinput.read_parameters(param, "./my_INPUT1");
        EXPECT_NO_THROW(readinput.read_parameters(param, "./my_INPUT1"));
        readinput.write_parameters(param, "./my_INPUT2");
        EXPECT_TRUE(compare_two_files("./my_INPUT1", "./my_INPUT2"));
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT1") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT2") == 0);
        readinput.clear();
    }
}

TEST_F(InputTest, Check)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    {
        std::ofstream emptyfile("empty_INPUT");
        emptyfile << "INPUT_PARAMETERS";
        emptyfile.close();

        Parameter param;
        readinput.read_parameters(param, "./empty_INPUT");
        readinput.write_parameters(param, "./INPUT.ref");
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
        readinput.clear();
    }

    ModuleIO::ReadInput::check_mode = true;
    Parameter param;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(readinput.read_parameters(param, "./INPUT.ref"), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("INPUT parameters have been successfully checked!"));
    EXPECT_TRUE(std::remove("./INPUT.ref") == 0);
}

TEST_F(InputTest, BasisIsFinalizedBeforeDependentDefaults)
{
    ModuleIO::ReadInput::check_mode = false;

    std::ofstream input("dependency_basis_INPUT");
    input << "INPUT_PARAMETERS\n"
          << "calculation nscf\n"
          << "basis_type lcao_in_pw\n"
          << "towannier90 true\n"
          << "device cpu\n";
    input.close();

    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    ASSERT_NO_THROW(readinput.read_parameters(param, "dependency_basis_INPUT"));
    EXPECT_EQ(param.inp.basis_type, "lcao");
    EXPECT_DOUBLE_EQ(param.inp.ecutwfc, 100.0);
    EXPECT_TRUE(std::remove("dependency_basis_INPUT") == 0);
}

TEST_F(InputTest, SpinIsFinalizedBeforeNupdown)
{
    ModuleIO::ReadInput::check_mode = false;

    std::ofstream input("dependency_spin_INPUT");
    input << "INPUT_PARAMETERS\n"
          << "noncolin true\n"
          << "nupdown 1\n";
    input.close();

    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    ASSERT_NO_THROW(readinput.read_parameters(param, "dependency_spin_INPUT"));
    EXPECT_EQ(param.inp.nspin, 4);
    EXPECT_TRUE(param.globalv.two_fermi);
    EXPECT_TRUE(std::remove("dependency_spin_INPUT") == 0);
}

TEST_F(InputTest, ScfThresholdIsFinalizedBeforeMixingRestart)
{
    ModuleIO::ReadInput::check_mode = false;

    std::ofstream input("dependency_mixing_INPUT");
    input << "INPUT_PARAMETERS\n"
          << "nspin 2\n"
          << "sc_mag_switch true\n"
          << "sc_scf_thr 10\n";
    input.close();

    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    ASSERT_NO_THROW(readinput.read_parameters(param, "dependency_mixing_INPUT"));
    EXPECT_DOUBLE_EQ(param.inp.scf_thr, 1.0e-9);
    EXPECT_DOUBLE_EQ(param.inp.mixing_restart, 1.0e-10);
    EXPECT_TRUE(std::remove("dependency_mixing_INPUT") == 0);
}

TEST_F(InputTest, SdftModeIsFinalizedBeforeBndpar)
{
    ModuleIO::ReadInput::check_mode = false;

    std::ofstream input("dependency_sdft_INPUT");
    input << "INPUT_PARAMETERS\n"
          << "esolver_type sdft\n"
          << "ks_solver cg\n"
          << "bndpar 2\n"
          << "nbands_sto 0\n";
    input.close();

    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    ASSERT_NO_THROW(readinput.read_parameters(param, "dependency_sdft_INPUT"));
    EXPECT_EQ(param.inp.esolver_type, "ksdft");
    EXPECT_EQ(param.inp.bndpar, 1);
    EXPECT_TRUE(std::remove("dependency_sdft_INPUT") == 0);
}

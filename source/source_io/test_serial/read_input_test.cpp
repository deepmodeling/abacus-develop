#include "source_io/module_parameter/read_input.h"

#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <stdexcept>

// mock
namespace GlobalV
{
int NPROC = 1;
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
        set_nproc(1);
    }

    void set_nproc(const int nproc)
    {
        GlobalV::NPROC = nproc;
    }

    void write_input(const std::string& filename, const std::string& parameters)
    {
        std::ofstream input(filename.c_str());
        input << "INPUT_PARAMETERS\n" << parameters;
    }

    void read_parameters(const std::string& filename, const std::string& parameters, Parameter& param)
    {
        write_input(filename, parameters);
        ModuleIO::ReadInput readinput(0);
        readinput.check_ntype_flag = false;
        try
        {
            readinput.read_parameters(param, filename);
        }
        catch (...)
        {
            std::remove(filename.c_str());
            throw;
        }
        EXPECT_EQ(std::remove(filename.c_str()), 0);
    }

    void expect_invalid_input(const std::string& filename,
                              const std::string& parameters,
                              const std::string& reason)
    {
        Parameter param;
        testing::internal::CaptureStdout();
        EXPECT_THROW(read_parameters(filename, parameters, param), std::runtime_error);
        const std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr(reason));
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
        EXPECT_EQ(param.inp.device, "cpu");
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

TEST_F(InputTest, RejectAutoDevice)
{
    std::ofstream input("auto_device_INPUT");
    input << "INPUT_PARAMETERS\n"
          << "device auto\n";
    input.close();

    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    EXPECT_THROW(readinput.read_parameters(param, "./auto_device_INPUT"), std::runtime_error);
    EXPECT_TRUE(std::remove("./auto_device_INPUT") == 0);
}

TEST_F(InputTest, ValidateNoncollinearSpin)
{
    Parameter valid_param;
    EXPECT_NO_THROW(read_parameters("noncolin_valid_INPUT", "noncolin 1\nnspin 4\n", valid_param));

    expect_invalid_input("noncolin_missing_nspin_INPUT",
                         "noncolin 1\n",
                         "nspin must be 4 when noncolin or lspinorb is enabled");
    expect_invalid_input("noncolin_invalid_nspin_INPUT",
                         "noncolin 1\nnspin 2\n",
                         "nspin must be 4 when noncolin or lspinorb is enabled");
    expect_invalid_input("soc_missing_nspin_INPUT",
                         "lspinorb 1\n",
                         "nspin must be 4 when noncolin or lspinorb is enabled");
}

TEST_F(InputTest, ValidateSdftStochasticBands)
{
    Parameter all_param;
    EXPECT_NO_THROW(read_parameters("sdft_all_INPUT", "esolver_type sdft\nnbands_sto all\n", all_param));
    EXPECT_EQ(all_param.inp.esolver_type, "sdft");

    expect_invalid_input("sdft_zero_INPUT",
                         "esolver_type sdft\nnbands_sto 00\n",
                         "nbands_sto should be in the range of 1 to 100000 or be all");
    expect_invalid_input("ksdft_zero_INPUT",
                         "esolver_type ksdft\nnbands_sto 0\n",
                         "nbands_sto should be in the range of 1 to 100000 or be all");
}

TEST_F(InputTest, ValidateBandParallelization)
{
    set_nproc(4);
    Parameter valid_param;
    EXPECT_NO_THROW(read_parameters("bndpar_valid_INPUT",
                                    "esolver_type sdft\nnbands_sto all\nbndpar 2\n",
                                    valid_param));
    EXPECT_EQ(valid_param.inp.bndpar, 2);

    Parameter bpcg_param;
    EXPECT_NO_THROW(read_parameters("bndpar_bpcg_INPUT", "ks_solver bpcg\nbndpar 2\n", bpcg_param));
    EXPECT_EQ(bpcg_param.inp.bndpar, 2);

    expect_invalid_input("bndpar_zero_INPUT", "bndpar 0\n", "bndpar must be greater than 0");
    expect_invalid_input("bndpar_wrong_solver_INPUT",
                         "bndpar 2\n",
                         "bndpar > 1 requires esolver_type=sdft or ks_solver=bpcg");

    set_nproc(3);
    expect_invalid_input("bndpar_not_divisible_INPUT",
                         "esolver_type sdft\nnbands_sto all\nbndpar 2\n",
                         "The number of processors can not be divided by bndpar");

    set_nproc(1);
    expect_invalid_input("bndpar_too_large_INPUT",
                         "esolver_type sdft\nnbands_sto all\nbndpar 2\n",
                         "bndpar can not exceed the number of MPI processes");
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

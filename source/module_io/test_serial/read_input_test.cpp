#include "module_io/read_input.h"

#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
/************************************************
 *  unit test of read_input_test.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Default
 *     - read empty INPUT file
 *   - Read:
 *     - read input parameters from INPUT files
 *   - Check:
 *     - check_mode = true
 *   - Item_test:
 *     - read in specific values for some items
 */

class InputTest : public testing::Test
{
  protected:
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

TEST_F(InputTest, Default)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    { // PW
        std::ofstream emptyfile("empty_INPUT");
        emptyfile << "INPUT_PARAMETERS";
        emptyfile.close();
        Parameter param;
        readinput.readin_parameters(param, "./empty_INPUT", "./my_INPUT1");
        EXPECT_TRUE(compare_two_files("./my_INPUT1", "./support/INPUT_pw.ref"));
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT1") == 0);
        readinput.clear();
    }
    { // LCAO
        std::ofstream emptyfile("empty_INPUT");
        emptyfile << "INPUT_PARAMETERS\n";
        emptyfile << "basis_type           lcao";
        emptyfile.close();
        Parameter param;
        readinput.readin_parameters(param, "./empty_INPUT", "./my_INPUT2");
        EXPECT_TRUE(compare_two_files("./my_INPUT2", "./support/INPUT_lcao.ref"));
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
        EXPECT_TRUE(std::remove("./my_INPUT2") == 0);
        readinput.clear();
    }
}

TEST_F(InputTest, Read)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    { // PW
        Parameter param;
        readinput.readin_parameters(param, "./support/INPUT_pw.ref", "./my_INPUT3");
        EXPECT_TRUE(compare_two_files("./my_INPUT3", "./support/INPUT_pw.ref"));
        EXPECT_TRUE(std::remove("./my_INPUT3") == 0);
        readinput.clear();
    }
    { // LCAO
        Parameter param;
        readinput.readin_parameters(param, "./support/INPUT_lcao.ref", "./my_INPUT4");
        EXPECT_TRUE(compare_two_files("./my_INPUT4", "./support/INPUT_lcao.ref"));
        EXPECT_TRUE(std::remove("./my_INPUT4") == 0);
        readinput.clear();
    }
}

TEST_F(InputTest, Check)
{
    ModuleIO::ReadInput::check_mode = true;
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(readinput.readin_parameters(param, "./support/INPUT_pw.ref"), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("INPUT parameters have been successfully checked!"));
}
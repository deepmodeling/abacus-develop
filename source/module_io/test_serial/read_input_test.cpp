#include "module_io/read_input.h"

#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
/************************************************
 *  unit test of input.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Default()
 *     - read empty INPUT file
 *   - Read()
 *     - read input parameters from INPUT files
 *   - Check()
 *     - check_mode = true
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
    Parameter param;
    readinput.readin_parameters(param, "./support/empty_INPUT", "./my_INPUT1");
    EXPECT_TRUE(compare_two_files("./my_INPUT1", "./support/INPUT.ref"));
    // EXPECT_TRUE(std::remove("./my_INPUT") == 0);
}

TEST_F(InputTest, Read)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    readinput.readin_parameters(param, "./support/INPUT.ref", "./my_INPUT2");
    EXPECT_TRUE(compare_two_files("./my_INPUT2", "./support/INPUT.ref"));
}
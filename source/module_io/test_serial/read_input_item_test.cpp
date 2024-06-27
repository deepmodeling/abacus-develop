#include "module_base/tool_quit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
/************************************************
 *  unit test of read_input_test_item.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Item_test:
 *     - read in specific values for some items
 */
#define private public
#include "module_io/input_item.h"
#include "module_io/read_input.h"

class InputTest : public testing::Test
{
  protected:
};

TEST_F(InputTest, Item_test)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    { // esolver_type
        std::ofstream file("testINPUT");
        file << "INPUT_PARAMETERS             \n";
        file << "esolver_type   none          \n";
        file.close();
        Parameter param;

        testing::internal::CaptureStdout();
        EXPECT_EXIT(readinput.read_parameters(param, "testINPUT"), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("esolver_type should be ksdft, sdft, ofdft, tddft, lj or dp."));
        EXPECT_TRUE(std::remove("testINPUT") == 0);
        readinput.clear();
    }
    { // nspin
        std::ofstream file("testINPUT");
        file << "INPUT_PARAMETERS             \n";
        file << "nspin              3         \n";
        file.close();
        Parameter param;

        testing::internal::CaptureStdout();
        EXPECT_EXIT(readinput.read_parameters(param, "testINPUT"), ::testing::ExitedWithCode(0), "");
        std::string output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("nspin should be 1, 2 or 4."));
        EXPECT_TRUE(std::remove("testINPUT") == 0);
        readinput.clear();
    }
    { // kspacing
        Parameter param;
        std::string word = "kspacing";
        auto it = std::find_if(
            readinput.input_lists.begin(),
            readinput.input_lists.end(),
            [&word](const std::pair<std::string, ModuleIO::Input_Item>& item) { return item.first == word; });
        if (it != readinput.input_lists.end())
        {
            param.input.kspacing = {-0.1, 0.1, 0.1};
            testing::internal::CaptureStdout();
            EXPECT_EXIT(it->second.checkvalue(it->second, param), ::testing::ExitedWithCode(0), "");
            std::string output = testing::internal::GetCapturedStdout();
            EXPECT_THAT(output, testing::HasSubstr("kspacing must > 0"));

            param.input.kspacing = {0, 1, 2};
            testing::internal::CaptureStdout();
            EXPECT_EXIT(it->second.checkvalue(it->second, param), ::testing::ExitedWithCode(0), "");
            output = testing::internal::GetCapturedStdout();
            EXPECT_THAT(output, testing::HasSubstr("kspacing must > 0"));
        }
    }
    //... It is too slow to test exit(0)
}
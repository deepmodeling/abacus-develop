#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/tool_quit.h"
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
#undef private

class InputTest : public testing::Test
{
  protected:
    std::vector<std::pair<std::string, ModuleIO::Input_Item>>::iterator find_label(
        const std::string& label,
        std::vector<std::pair<std::string, ModuleIO::Input_Item>>& bcastfuncs)
    {
        auto it = std::find_if(
            bcastfuncs.begin(),
            bcastfuncs.end(),
            [&label](const std::pair<std::string, ModuleIO::Input_Item>& item) { return item.first == label; });
        return it;
    }
};


TEST_F(InputTest, Item_test)
{
    ModuleIO::ReadInput readinput(0);
    Parameter param;
    std::string output = "";
    readinput.check_ntype_flag = false;

    { 
        param.input.suffix = "test";
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.global_out_dir, "OUT.test/");
        EXPECT_EQ(param.sys.global_stru_dir, "OUT.test/STRU/");
        EXPECT_EQ(param.sys.global_matrix_dir, "OUT.test/matrix/");

        param.input.basis_type = "lcao";
        param.input.gamma_only = true;
        param.input.esolver_type = "tddft";
        param.input.nspin = 2;
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.gamma_only_local, 0);

        param.input.nspin = 4;
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.domag, 1);
        EXPECT_EQ(param.sys.domag_z, 1);
        EXPECT_EQ(param.sys.npol, 2);

        param.input.nspin = 1;
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.domag, 0);
        EXPECT_EQ(param.sys.domag_z, 0);
        EXPECT_EQ(param.sys.npol, 1);

        param.input.deepks_scf = true;
        param.input.deepks_out_labels = true;
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.deepks_setorb, 1);
    }


}
TEST_F(InputTest, Item_test2)
{
        ModuleIO::ReadInput readinput(0);
    Parameter param;
    std::string output = "";
    readinput.check_ntype_flag = false;

    { 
        param.input.device = "gpu";
        param.input.basis_type = "lcao_in_pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(readinput.set_globalv(param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("not support"));

        param.input.basis_type = "pw";
        param.input.device = "gpu";
        readinput.set_globalv(param);
        EXPECT_EQ(param.sys.device_flag, "gpu");
    }
}
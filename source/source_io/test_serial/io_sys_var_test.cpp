#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/tool_quit.h"
/************************************************
 *  unit test of read_input_test_item.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Item_test:
 *     - read in specific values for some items
 */
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/test_parameters.h"
#define private public
#include "source_io/module_parameter/input_item.h"
#include "source_io/module_parameter/read_input.h"
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
ModuleIO::ReadInput readinput(0);
Parameter param;
std::string output = "";

TEST_F(InputTest, Item_test)
{
    readinput.check_ntype_flag = false;

    { 
        TestParameters::input(param).suffix = "test";
        readinput.set_global_dir(param.inp, TestParameters::sys(param));

        EXPECT_EQ(TestParameters::sys(param).global_out_dir, "OUT.test/");
        EXPECT_EQ(TestParameters::sys(param).global_stru_dir, "OUT.test/STRU/");
        EXPECT_EQ(TestParameters::sys(param).global_matrix_dir, "OUT.test/matrix/");

        readinput.set_globalv(param.inp, TestParameters::sys(param));
    
        TestParameters::input(param).basis_type = "lcao";
        TestParameters::input(param).gamma_only = true;
        TestParameters::input(param).esolver_type = "tddft";
        TestParameters::input(param).nspin = 2;
        readinput.set_globalv(param.inp, TestParameters::sys(param));
        EXPECT_EQ(TestParameters::sys(param).gamma_only_local, 0);
        
        TestParameters::input(param).deepks_scf = true;
        TestParameters::input(param).deepks_out_labels = true;
        readinput.set_globalv(param.inp, TestParameters::sys(param));
        EXPECT_EQ(TestParameters::sys(param).deepks_setorb, 1);

        TestParameters::input(param).nspin = 4;
        TestParameters::input(param).noncolin = true;
        readinput.set_globalv(param.inp, TestParameters::sys(param));
        EXPECT_EQ(TestParameters::sys(param).domag, 1);
        EXPECT_EQ(TestParameters::sys(param).domag_z, 0);
        EXPECT_EQ(TestParameters::sys(param).npol, 2);

        TestParameters::input(param).nspin = 1;
        TestParameters::input(param).lspinorb = true;
        TestParameters::input(param).noncolin = false;
        readinput.set_globalv(param.inp, TestParameters::sys(param));
        EXPECT_EQ(TestParameters::sys(param).domag, 0);
        EXPECT_EQ(TestParameters::sys(param).domag_z, 0);
        EXPECT_EQ(TestParameters::sys(param).npol, 1);

        
    }
}
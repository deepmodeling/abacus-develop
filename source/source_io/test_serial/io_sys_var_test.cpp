#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/input_parameter.h"
#define private public
#include "source_io/module_parameter/read_input.h"
#undef private
#include "source_io/module_parameter/system_parameter.h"
/************************************************
 *  unit test of read_input_test_item.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Item_test:
 *     - read in specific values for some items
 */

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

/// @brief Friend helper to invoke ReadInput private setters in unit tests.
/// @details ReadInput grants friend access to TestParameters so the test can
/// drive its private setup routines without `#define private public`. The
/// class must stay at global scope to match the friend declaration in
/// read_input.h; an anonymous-namespace class would not be the friend.
class TestParameters
{
  public:
    static void set_global_dir(ModuleIO::ReadInput& readinput, const Input_para& inp, System_para& sys)
    {
        readinput.set_global_dir(inp, sys);
    }
    static void set_globalv(ModuleIO::ReadInput& readinput, const Input_para& inp, System_para& sys)
    {
        readinput.set_globalv(inp, sys);
    }
};
ModuleIO::ReadInput readinput(0);
Input_para input;
System_para sys;
std::string output = "";

TEST_F(InputTest, Item_test)
{
    readinput.check_ntype_flag = false;

    {
        input.suffix = "test";
        TestParameters::set_global_dir(readinput, input, sys);

        EXPECT_EQ(sys.global_out_dir, "OUT.test/");
        EXPECT_EQ(sys.global_stru_dir, "OUT.test/STRU/");
        EXPECT_EQ(sys.global_matrix_dir, "OUT.test/matrix/");

        TestParameters::set_globalv(readinput, input, sys);

        input.basis_type = "lcao";
        input.gamma_only = true;
        input.esolver_type = "tddft";
        input.nspin = 2;
        TestParameters::set_globalv(readinput, input, sys);
        EXPECT_EQ(sys.gamma_only_local, 0);

        input.deepks_scf = true;
        input.deepks_out_labels = true;
        TestParameters::set_globalv(readinput, input, sys);
        EXPECT_EQ(sys.deepks_setorb, 1);

        input.nspin = 4;
        input.noncolin = true;
        TestParameters::set_globalv(readinput, input, sys);
        EXPECT_EQ(sys.domag, 1);
        EXPECT_EQ(sys.domag_z, 0);
        EXPECT_EQ(sys.npol, 2);

        input.nspin = 1;
        input.lspinorb = true;
        input.noncolin = false;
        TestParameters::set_globalv(readinput, input, sys);
        EXPECT_EQ(sys.domag, 0);
        EXPECT_EQ(sys.domag_z, 0);
        EXPECT_EQ(sys.npol, 1);
    }
}

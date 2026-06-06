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

static bool all_ks_run_after_set_globalv(const std::string& ks_solver, const int bndpar)
{
    Parameter local_param;
    local_param.input.ks_solver = ks_solver;
    local_param.input.bndpar = bndpar;
    local_param.sys.all_ks_run = true;
    readinput.set_globalv(local_param.inp, local_param.sys);
    return local_param.sys.all_ks_run;
}

TEST_F(InputTest, Item_test)
{
    readinput.check_ntype_flag = false;

    { 
        param.input.suffix = "test";
        readinput.set_global_dir(param.inp, param.sys);

        EXPECT_EQ(param.sys.global_out_dir, "OUT.test/");
        EXPECT_EQ(param.sys.global_stru_dir, "OUT.test/STRU/");
        EXPECT_EQ(param.sys.global_matrix_dir, "OUT.test/matrix/");

        readinput.set_globalv(param.inp, param.sys);
    
        param.input.basis_type = "lcao";
        param.input.gamma_only = true;
        param.input.esolver_type = "tddft";
        param.input.nspin = 2;
        readinput.set_globalv(param.inp, param.sys);
        EXPECT_EQ(param.sys.gamma_only_local, 0);
        
        param.input.deepks_scf = true;
        param.input.deepks_out_labels = true;
        readinput.set_globalv(param.inp, param.sys);
        EXPECT_EQ(param.sys.deepks_setorb, 1);

        param.input.nspin = 4;
        param.input.noncolin = true;
        readinput.set_globalv(param.inp, param.sys);
        EXPECT_EQ(param.sys.domag, 1);
        EXPECT_EQ(param.sys.domag_z, 0);
        EXPECT_EQ(param.sys.npol, 2);

        param.input.nspin = 1;
        param.input.lspinorb = true;
        param.input.noncolin = false;
        readinput.set_globalv(param.inp, param.sys);
        EXPECT_EQ(param.sys.domag, 0);
        EXPECT_EQ(param.sys.domag_z, 0);
        EXPECT_EQ(param.sys.npol, 1);

        
    }
}

TEST_F(InputTest, BandParallelAllKsRun)
{
    EXPECT_FALSE(all_ks_run_after_set_globalv("cg", 2));
    EXPECT_FALSE(all_ks_run_after_set_globalv("dav", 2));
    EXPECT_TRUE(all_ks_run_after_set_globalv("bpcg", 2));
    EXPECT_TRUE(all_ks_run_after_set_globalv("lobpcg", 2));
}

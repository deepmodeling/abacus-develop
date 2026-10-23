#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/fs_compat.h"
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
        // out_app_flag defaults to 1, so stru dir maps to global out dir
        EXPECT_EQ(sys.global_stru_dir, "OUT.test/");
        // out_app_flag defaults to 1, so matrix dir maps to global out dir
        EXPECT_EQ(sys.global_matrix_dir, "OUT.test/");
        EXPECT_EQ(sys.global_wfc_dir, "OUT.test/");
        EXPECT_EQ(sys.global_psi_dir, "OUT.test/");
        EXPECT_EQ(sys.global_chg_dir, "OUT.test/");
        EXPECT_EQ(sys.global_pot_dir, "OUT.test/");
        EXPECT_EQ(sys.global_elf_dir, "OUT.test/");
        EXPECT_EQ(sys.global_dos_dir, "OUT.test/");
        EXPECT_EQ(sys.global_dip_dir, "OUT.test/");
        EXPECT_EQ(sys.global_mat_dir, "OUT.test/");

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

TEST_F(InputTest, SetGlobalDirMdRestart)
{
    readinput.check_ntype_flag = false;

    const std::string case_dir = testing::TempDir() + "/abacus_md_restart_io/";
    ModuleBase::make_directory(case_dir);

    // Restart from a user-provided read_file_dir holding the legacy
    // "STRU_MD_{istep}" naming: the structure file must be picked up.
    const std::string legacy_dir = case_dir + "legacy/";
    ModuleBase::make_directory(legacy_dir);
    const std::string restart_txt = legacy_dir + "Restart_md.txt";
    std::ofstream restart_file(restart_txt.c_str());
    restart_file << 5 << std::endl;
    restart_file << 0.03 << std::endl;
    restart_file.close();
    const std::string stru_legacy = legacy_dir + "STRU_MD_5";
    std::ofstream stru_file_legacy(stru_legacy.c_str());
    stru_file_legacy << "ATOMIC_SPECIES" << std::endl;
    stru_file_legacy.close();

    input.suffix = "autotest";
    input.calculation = "md";
    input.mdp.md_restart = true;
    input.read_file_dir = legacy_dir;
    TestParameters::set_global_dir(readinput, input, sys);
    EXPECT_EQ(sys.global_in_stru, stru_legacy);

    // Restart from a read_file_dir holding the current "STRU{istep + 1}"
    // naming written by run_md (1-based, unified with relax).
    const std::string newnaming_dir = case_dir + "newnaming/";
    ModuleBase::make_directory(newnaming_dir);
    const std::string restart_txt_new = newnaming_dir + "Restart_md.txt";
    std::ofstream restart_file_new(restart_txt_new.c_str());
    restart_file_new << 5 << std::endl;
    restart_file_new << 0.03 << std::endl;
    restart_file_new.close();
    const std::string stru_new = newnaming_dir + "STRU6";
    std::ofstream stru_file_new(stru_new.c_str());
    stru_file_new << "ATOMIC_SPECIES" << std::endl;
    stru_file_new.close();

    input.read_file_dir = newnaming_dir;
    TestParameters::set_global_dir(readinput, input, sys);
    EXPECT_EQ(sys.global_in_stru, stru_new);

    // restore the shared input state so other tests are not affected
    input.calculation = "scf";
    input.mdp.md_restart = false;
    input.read_file_dir = "OUT.$suffix";
}

TEST_F(InputTest, SetGlobalDirMode2)
{
    readinput.check_ntype_flag = false;

    // out_app_flag = 2 routes every output class to its unified
    // subdirectory below OUT.<suffix>/
    input.suffix = "modetwo";
    input.out_app_flag = 2;
    TestParameters::set_global_dir(readinput, input, sys);

    EXPECT_EQ(sys.global_out_dir, "OUT.modetwo/");
    EXPECT_EQ(sys.global_stru_dir, "OUT.modetwo/STRU/");
    EXPECT_EQ(sys.global_psi_dir, "OUT.modetwo/PSI/");
    EXPECT_EQ(sys.global_chg_dir, "OUT.modetwo/CHG/");
    EXPECT_EQ(sys.global_pot_dir, "OUT.modetwo/POT/");
    EXPECT_EQ(sys.global_elf_dir, "OUT.modetwo/ELF/");
    EXPECT_EQ(sys.global_dos_dir, "OUT.modetwo/DOS/");
    EXPECT_EQ(sys.global_dip_dir, "OUT.modetwo/DIP/");
    EXPECT_EQ(sys.global_mat_dir, "OUT.modetwo/MAT/");
    // legacy aliases follow the unified directories
    EXPECT_EQ(sys.global_wfc_dir, sys.global_psi_dir);
    EXPECT_EQ(sys.global_matrix_dir, sys.global_mat_dir);

    // restore the shared input state so other tests are not affected
    input.suffix = "test";
    input.out_app_flag = 1;
}

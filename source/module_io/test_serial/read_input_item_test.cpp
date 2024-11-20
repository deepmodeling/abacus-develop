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
        std::vector<std::pair<std::string, ModuleIO::Input_Item>>& input_lists)
    {
        auto it = std::find_if(
            input_lists.begin(),
            input_lists.end(),
            [&label](const std::pair<std::string, ModuleIO::Input_Item>& item) { return item.first == label; });
        return it;
    }
};

TEST_F(InputTest, Item_test)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;

    std::string output = "";

    { // calculation
        auto it = find_label("calculation", readinput.input_lists);
        param.input.calculation = "get_pchg";
        param.input.basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.calculation = "gen_bessel";
        param.input.basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.calculation = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }

    { // esolver_type
        auto it = find_label("esolver_type", readinput.input_lists);
        param.input.esolver_type = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.esolver_type = "dp";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.esolver_type = "lr";
        param.input.calculation = "scf";
        it = find_label("esolver_type", readinput.input_lists);
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.calculation, "nscf");
    }
    { // nspin
        auto it = find_label("nspin", readinput.input_lists);
        param.input.nspin = 0;
        param.input.noncolin = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nspin, 4);

        param.input.nspin = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // kspacing
        auto it = find_label("kspacing", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.kspacing[0], 1);
        EXPECT_EQ(param.input.kspacing[1], 1);
        EXPECT_EQ(param.input.kspacing[2], 1);
        it->second.str_values = {"1", "2"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.kspacing = {0, -1, 1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.kspacing = {0, 1, 2};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nbands
        auto it = find_label("nbands", readinput.input_lists);
        param.input.nbands = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // symmetry
        auto it = find_label("symmetry", readinput.input_lists);
        param.input.symmetry = "default";
        param.input.gamma_only = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "0");

        param.input.symmetry = "default";
        param.input.gamma_only = false;
        param.input.calculation = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "1");

        param.input.calculation = "md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "0");

        param.input.symmetry = "default";
        param.input.efield_flag = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "0");

        param.input.qo_switch = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "-1");

        param.input.symmetry = "default";
        param.input.berry_phase = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.symmetry, "-1");
    }
    { // nelec
        auto it = find_label("nelec", readinput.input_lists);
        param.input.nelec = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.nelec = 100;
        param.input.nbands = 5;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bndpar
        auto it = find_label("bndpar", readinput.input_lists);
        param.input.esolver_type = "ksdft";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.bndpar, 1);

        param.input.esolver_type = "sdft";
        param.input.bndpar = 2;
        GlobalV::NPROC = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.bndpar, 1);
    }
    { // dft_plus_dmft
        auto it = find_label("dft_plus_dmft", readinput.input_lists);
        param.input.basis_type = "pw";
        param.input.dft_plus_dmft = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // mem_saver
        auto it = find_label("mem_saver", readinput.input_lists);
        param.input.mem_saver = 1;
        param.input.calculation = "scf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mem_saver, 0);

        param.input.mem_saver = 1;
        param.input.calculation = "relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mem_saver, 0);
    }
    { // diag_proc
        auto it = find_label("diago_proc", readinput.input_lists);
        param.input.diago_proc = 0;
        GlobalV::NPROC = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.diago_proc, 1);
    }
    { // cal_force
        auto it = find_label("cal_force", readinput.input_lists);
        param.input.calculation = "cell-relax";
        param.input.cal_force = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.cal_force, true);

        param.input.calculation = "get_wf";
        param.input.cal_force = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.cal_force, false);
    }
    { // ecutrho
        auto it = find_label("ecutrho", readinput.input_lists);
        param.input.ecutwfc = 1;
        param.input.ecutrho = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ecutrho, 4);
        param.input.nx = 0;
        param.input.ecutrho = 5;
        param.input.ecutwfc = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.sys.double_grid, true);

        param.input.ecutwfc = 1;
        param.input.ecutrho = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // pw_diag_thr
        auto it = find_label("pw_diag_thr", readinput.input_lists);
        param.input.pw_diag_thr = 1.0e-2;
        param.input.calculation = "get_S";
        param.input.basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.pw_diag_thr, 1.0e-5);
    }
    { // nb2d
        auto it = find_label("nb2d", readinput.input_lists);
        param.input.nb2d = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // scf_thr
        auto it = find_label("scf_thr", readinput.input_lists);
        param.input.scf_thr = -1;
        param.input.basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.scf_thr, 1.0e-7);

        param.input.scf_thr = -1;
        param.input.basis_type = "pw";
        param.input.calculation = "nscf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.scf_thr, 1.0e-6);

        param.input.scf_thr = -1;
        param.input.basis_type = "pw";
        param.input.calculation = "scf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.scf_thr, 1.0e-9);
    }
    { // scf_thr_type
        auto it = find_label("scf_thr_type", readinput.input_lists);
        param.input.scf_thr_type = -1;
        param.input.basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.scf_thr_type, 2);

        param.input.scf_thr_type = -1;
        param.input.basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.scf_thr_type, 1);
    }
    { // init_wfc
        auto it = find_label("init_wfc", readinput.input_lists);
        param.input.init_wfc = "atomic";
        param.input.calculation = "get_pchg";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.init_wfc, "file");

        param.input.init_wfc = "atomic";
        param.input.basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.init_wfc, "nao");
    }
    { // psi_initializer
        auto it = find_label("psi_initializer", readinput.input_lists);
        param.input.psi_initializer = false;
        param.input.basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.psi_initializer, true);
    }
    { // init_chg
        auto it = find_label("init_chg", readinput.input_lists);
        param.input.init_chg = "get_pchg";
        param.input.init_chg = "";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.init_chg, "atomic");

        param.input.init_chg = "";
        param.input.calculation = "nscf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.init_chg, "file");

        param.input.init_chg = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // chg_extrap
        auto it = find_label("chg_extrap", readinput.input_lists);
        param.input.chg_extrap = "default";
        param.input.calculation = "md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.chg_extrap, "second-order");

        param.input.chg_extrap = "default";
        param.input.calculation = "relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.chg_extrap, "first-order");

        param.input.chg_extrap = "default";
        param.input.calculation = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.chg_extrap, "atomic");

        param.input.chg_extrap = "none";
        param.input.calculation = "get_wf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.chg_extrap, "atomic");
    }
    { // out_chg
        auto it = find_label("out_chg", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_chg = {0};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_chg[0], 1);
    }
    { // out_pot
        auto it = find_label("out_pot", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_pot = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_pot, 0);
    }
    { // out_dos
        auto it = find_label("out_dos", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_dos = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_dos, 0);

        param.input.out_dos = 3;
        param.input.symmetry = "1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.basis_type = "pw";
        param.input.out_dos = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // out_band
        auto it = find_label("out_band", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.out_band[0], 1);
        EXPECT_EQ(param.input.out_band[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.out_band[0], 1);
        EXPECT_EQ(param.input.out_band[1], 2);

        it->second.str_values = {"1", "2", "3"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.calculation = "get_wf";
        param.input.out_band = {1, 2};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_band[0], 0);
    }
    { // out_proj_band
        auto it = find_label("out_proj_band", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_proj_band = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_proj_band, false);

        param.input.basis_type = "pw";
        param.input.out_proj_band = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // read_file_dir
        auto it = find_label("read_file_dir", readinput.input_lists);
        param.input.read_file_dir = "auto";
        param.input.suffix = "test";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.read_file_dir, "OUT.test/");

        param.input.read_file_dir = "test";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.read_file_dir, "test/");
    }
    { // nx
        auto it = find_label("nx", readinput.input_lists);
        param.input.nx = 1;
        param.input.ny = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ny
        auto it = find_label("ny", readinput.input_lists);
        param.input.ny = 1;
        param.input.nz = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nz
        auto it = find_label("nz", readinput.input_lists);
        param.input.nz = 1;
        param.input.nx = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndx
        auto it = find_label("ndx", readinput.input_lists);
        param.input.ndx = 2;
        param.input.nx = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.sys.double_grid, true);

        param.input.ndx = 1;
        param.input.ndy = 0;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ndx = 1;
        param.input.nx = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndy
        auto it = find_label("ndy", readinput.input_lists);
        param.input.ndy = 2;
        param.input.ny = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.sys.double_grid, true);

        param.input.ndy = 1;
        param.input.ndz = 0;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ndy = 1;
        param.input.ny = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndz
        auto it = find_label("ndz", readinput.input_lists);
        param.input.ndz = 2;
        param.input.nz = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.sys.double_grid, true);

        param.input.ndz = 1;
        param.input.nz = 2;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ndz = 1;
        param.input.nz = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // cell_factor
        auto it = find_label("cell_factor", readinput.input_lists);
        param.input.calculation = "cell-relax";
        param.input.cell_factor = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.cell_factor, 2.0);
    }
    { // ks_sovler
        auto it = find_label("ks_solver", readinput.input_lists);
        param.input.ks_solver = "default";
        param.input.basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ks_solver, "cg");

        param.input.ks_solver = "default";
        param.input.basis_type = "lcao";
        param.input.device = "gpu";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ks_solver, "cusolver");
#ifdef __ELPA
        param.input.towannier90 = true;
        param.input.basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ks_solver, "genelpa");
#else
#ifdef __MPI
        param.input.towannier90 = true;
        param.input.basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ks_solver, "scalapack_gvx");
#else
        param.input.towannier90 = true;
        param.input.basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.ks_solver, "lapack");
#endif
#endif
        param.input.ks_solver = "default";
        param.input.basis_type = "lcao";
        param.input.device = "cpu";
        param.input.ks_solver = "lapack";

        param.input.ks_solver = "genelpa";
        param.input.basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ks_solver = "cg";
        param.input.basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ks_solver = "scalapack_gvx";
        param.input.basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ks_solver = "cg";
        param.input.basis_type = "lcao_in_pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // relax_nmax
        auto it = find_label("relax_nmax", readinput.input_lists);
        param.input.calculation = "scf";
        param.input.relax_nmax = 5;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.relax_nmax, 1);

        param.input.calculation = "relax";
        param.input.relax_nmax = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.relax_nmax, 50);
    }
    { // out_stru
        auto it = find_label("out_stru", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_stru = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_stru, false);
    }
    { // cal_stress
        auto it = find_label("cal_stress", readinput.input_lists);
        param.input.calculation = "md";
        param.input.cal_stress = false;
        param.input.esolver_type = "lj";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.cal_stress, true);
    }
    { // fixed_axes
        auto it = find_label("fixed_axes", readinput.input_lists);
        param.input.fixed_axes = "shape";
        param.input.relax_new = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.fixed_axes = "volume";
        param.input.relax_new = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // fixed_ibrav
        auto it = find_label("fixed_ibrav", readinput.input_lists);
        param.input.fixed_ibrav = true;
        param.input.relax_new = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.fixed_ibrav = true;
        param.input.latname = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // fixed_atoms
        auto it = find_label("fixed_atoms", readinput.input_lists);
        param.input.fixed_atoms = true;
        param.input.calculation = "relax";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // relax_method
        auto it = find_label("relax_method", readinput.input_lists);
        param.input.relax_method = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { //relax_new
        auto it = find_label("relax_new", readinput.input_lists);
        param.input.relax_new = true;
        param.input.relax_method = "cg";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.relax_new, true);

        param.input.relax_new = true;
        param.input.relax_method = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.relax_new, false);
    }
    { // force_thr
        auto it = find_label("force_thr", readinput.input_lists);
        param.input.force_thr = -1;
        param.input.force_thr_ev = -1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.force_thr, 1.0e-3);
        EXPECT_EQ(param.input.force_thr_ev, 1.0e-3 * 13.6058 / 0.529177);

        param.input.force_thr = -1;
        param.input.force_thr_ev = 1.0e-3 * 13.6058 / 0.529177;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.force_thr, 1.0e-3);

        param.input.force_thr = 1.0e-3;
        param.input.force_thr_ev = 1.0e-3 * 13.6058 / 0.529177;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.force_thr, 1.0e-3);
        EXPECT_EQ(param.input.force_thr_ev, 1.0e-3 * 13.6058 / 0.529177);
    }
    { // out_level
        auto it = find_label("out_level", readinput.input_lists);
        param.input.out_level = "0";
        param.input.calculation = "md";
        param.sys.out_md_control = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_level, "m");
    }
    { // out_dm
        auto it = find_label("out_dm", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_dm = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_dm, false);

        param.sys.gamma_only_local = false;
        param.input.out_dm = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // out_dm1
        auto it = find_label("out_dm1", readinput.input_lists);
        param.input.calculation = "get_wf";
        param.input.out_dm1 = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_dm1, false);

        param.sys.gamma_only_local = true;
        param.input.out_dm1 = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // use_paw
        auto it = find_label("use_paw", readinput.input_lists);
        param.input.use_paw = true;
        param.input.basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
        param.input.use_paw = true;
        param.input.dft_functional = "default";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // method_sto
        auto it = find_label("method_sto", readinput.input_lists);
        param.input.method_sto = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nbands_sto
        auto it = find_label("nbands_sto", readinput.input_lists);
        param.input.esolver_type = "sdft";

        it->second.str_values = {"all"};
        it->second.read_value(it->second, param);
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nbands_sto, 0);
        EXPECT_EQ(param.input.esolver_type, "sdft");

        it->second.str_values = {"8"};
        it->second.read_value(it->second, param);
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nbands_sto, 8);
        EXPECT_EQ(param.input.esolver_type, "sdft");

        it->second.str_values = {"0"};
        it->second.read_value(it->second, param);
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nbands_sto, 0);
        EXPECT_EQ(param.input.esolver_type, "ksdft");

        param.input.nbands_sto = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        it->second.str_values = {"all"};
        it->second.get_final_value(it->second, param);
        EXPECT_EQ(it->second.final_value.str(), "all");
        it->second.final_value.str("");

        it->second.str_values = {};
        param.input.nbands_sto = 256;
        it->second.get_final_value(it->second, param);
        EXPECT_EQ(it->second.final_value.str(), "256");
    }
    { // basis_type
        auto it = find_label("basis_type", readinput.input_lists);
        param.input.basis_type = "lcao_in_pw";
        param.input.towannier90 = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.basis_type, "lcao");

        param.input.basis_type = "gauss";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // gamma_only
        auto it = find_label("gamma_only", readinput.input_lists);
        param.input.basis_type = "pw";
        param.input.gamma_only = true;
        std::string filename = param.input.kpoint_file;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.gamma_only, false);
        std::ifstream ifs(filename);
        std::string line;
        std::getline(ifs, line);
        ifs.close();
        EXPECT_EQ(line, "K_POINTS");

        param.input.basis_type = "lcao";
        param.input.gamma_only = true;
        param.input.nspin = 4;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

    }
    { // out_mat_r
        auto it = find_label("out_mat_r", readinput.input_lists);
        param.input.esolver_type = "lcao";
        param.input.out_mat_r = true;
        param.sys.gamma_only_local = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("available"));
    }
    { // lcao_ecut
        auto it = find_label("lcao_ecut", readinput.input_lists);
        param.input.lcao_ecut = 0;
        param.input.ecutwfc = 1;
        param.input.basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.lcao_ecut, 1);
    }
    { // out_mat_hs
        auto it = find_label("out_mat_hs", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.out_mat_hs[0], 1);
        EXPECT_EQ(param.input.out_mat_hs[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.out_mat_hs[0], 1);
        EXPECT_EQ(param.input.out_mat_hs[1], 2);

        it->second.str_values = {"1", "2", "3"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.out_mat_hs = {0};
        param.input.qo_switch = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_mat_hs[0], 1);
    }
}
TEST_F(InputTest, Item_test2)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;

    std::string output = "";
    { // out_mat_dh
        auto it = find_label("out_mat_dh", readinput.input_lists);
        param.input.out_mat_dh = true;
        param.input.nspin = 4;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
#ifndef __USECNPY
    { // out_hr_npz
        auto it = find_label("out_hr_npz", readinput.input_lists);
        param.input.out_hr_npz = true;

        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // out_dm_npz
        auto it = find_label("out_dm_npz", readinput.input_lists);
        param.input.out_dm_npz = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
#endif
    { // out_interval
        auto it = find_label("out_interval", readinput.input_lists);
        param.input.out_interval = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // dm_to_rho
        auto it = find_label("dm_to_rho", readinput.input_lists);
        param.input.dm_to_rho = true;
        GlobalV::NPROC = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.dm_to_rho = true;
        param.input.gamma_only = true;
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

#ifndef __USECNPY
        param.input.dm_to_rho = true;
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
#endif
    }
    { // out_wfc_lcao
        auto it = find_label("out_wfc_lcao", readinput.input_lists);
        param.input.out_wfc_lcao = false;
        param.input.qo_switch = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.out_wfc_lcao, true);

        param.input.out_wfc_lcao = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.out_wfc_lcao = 1;
        param.input.basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bx
        auto it = find_label("bx", readinput.input_lists);
        param.input.bx = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.bx = 2;
        param.input.basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.bx, 1);
        EXPECT_EQ(param.input.by, 1);
        EXPECT_EQ(param.input.bz, 1);
    }
    { // by
        auto it = find_label("by", readinput.input_lists);
        param.input.by = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bz
        auto it = find_label("bz", readinput.input_lists);
        param.input.bz = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // mixing_beta
        auto it = find_label("mixing_beta", readinput.input_lists);
        param.input.mixing_beta = -1;
        param.input.nspin = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mixing_beta, 0.8);

        param.input.mixing_beta = -1;
        param.input.nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mixing_beta, 0.4);
        EXPECT_EQ(param.input.mixing_beta_mag, 1.6);
        EXPECT_EQ(param.input.mixing_gg0_mag, 0.0);

        param.input.mixing_beta = -1;
        param.input.nspin = 4;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mixing_beta, 0.4);
        EXPECT_EQ(param.input.mixing_beta_mag, 1.6);
        EXPECT_EQ(param.input.mixing_gg0_mag, 0.0);
    }
    { // mixing_beta_mag
        auto it = find_label("mixing_beta_mag", readinput.input_lists);
        param.input.mixing_beta = 0.3;
        param.input.mixing_beta_mag = -1;
        param.input.nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mixing_beta_mag, 1.2);

        param.input.mixing_beta = 0.5;
        param.input.mixing_beta_mag = -1;
        param.input.nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mixing_beta_mag, 1.6);
    }
    { // dip_cor_flag
        auto it = find_label("dip_cor_flag", readinput.input_lists);
        param.input.dip_cor_flag = true;
        param.input.efield_flag = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // efield_dir
        auto it = find_label("efield_dir", readinput.input_lists);
        param.input.gate_flag = true;
        param.input.efield_flag = true;
        param.input.dip_cor_flag = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_s6
        auto it = find_label("vdw_s6", readinput.input_lists);
        param.input.vdw_s6 = "default";
        param.input.vdw_method = "d2";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.vdw_s6, "0.75");

        // dftd3 parameter will not get its value here
        param.input.vdw_s6 = "default";
        param.input.vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s6, "1.0");
        EXPECT_EQ(param.input.vdw_s6, "default");
    }
    { // vdw_s8
        auto it = find_label("vdw_s8", readinput.input_lists);
        param.input.vdw_s8 = "default";
        param.input.vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s8, "0.722");
        EXPECT_EQ(param.input.vdw_s8, "default");

        param.input.vdw_s8 = "default";
        param.input.vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s8, "0.7875");
        EXPECT_EQ(param.input.vdw_s8, "default");
    }
    { // vdw_a1
        auto it = find_label("vdw_a1", readinput.input_lists);
        param.input.vdw_a1 = "default";
        param.input.vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a1, "1.217");
        EXPECT_EQ(param.input.vdw_a1, "default");

        param.input.vdw_a1 = "default";
        param.input.vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a1, "0.4289");
        EXPECT_EQ(param.input.vdw_a1, "default");
    }
    { // vdw_a2
        auto it = find_label("vdw_a2", readinput.input_lists);
        param.input.vdw_a2 = "default";
        param.input.vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a2, "1.0");
        EXPECT_EQ(param.input.vdw_a2, "default");

        param.input.vdw_a2 = "default";
        param.input.vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a2, "4.4407");
        EXPECT_EQ(param.input.vdw_a2, "default");
    }
    { // vdw_c6_unit
        auto it = find_label("vdw_c6_unit", readinput.input_lists);
        param.input.vdw_C6_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_r0_unit
        auto it = find_label("vdw_r0_unit", readinput.input_lists);
        param.input.vdw_R0_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_type
        auto it = find_label("vdw_cutoff_type", readinput.input_lists);
        param.input.vdw_cutoff_type = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_radius
        auto it = find_label("vdw_cutoff_radius", readinput.input_lists);
        param.input.vdw_cutoff_radius = "default";
        param.input.vdw_method = "d2";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.vdw_cutoff_radius, "56.6918");

        param.input.vdw_cutoff_radius = "default";
        param.input.vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.vdw_cutoff_radius, "95");

        param.input.vdw_cutoff_radius = "default";
        param.input.vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.vdw_cutoff_radius, "95");

        param.input.vdw_cutoff_radius = "default";
        param.input.vdw_method = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.vdw_cutoff_radius, "0");

        param.input.vdw_cutoff_radius = "-1";
        param.input.vdw_method = "d2";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_radius_unit
        auto it = find_label("vdw_radius_unit", readinput.input_lists);
        param.input.vdw_radius_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cn_thr
        auto it = find_label("vdw_cn_thr", readinput.input_lists);
        param.input.vdw_cn_thr = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cn_thr_unit
        auto it = find_label("vdw_cn_thr_unit", readinput.input_lists);
        param.input.vdw_cn_thr_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_period
        auto it = find_label("vdw_cutoff_period", readinput.input_lists);
        it->second.str_values = {"1", "1", "1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(param.input.vdw_cutoff_period[0], 1);
        EXPECT_EQ(param.input.vdw_cutoff_period[1], 1);
        EXPECT_EQ(param.input.vdw_cutoff_period[2], 1);

        it->second.str_values = {"1", "1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.vdw_cutoff_period = {-1, 1, 1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_hybrid_alpha
        auto it = find_label("exx_hybrid_alpha", readinput.input_lists);
        param.input.exx_hybrid_alpha = "default";
        param.input.dft_functional = "HF";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_hybrid_alpha, "1");

        param.input.exx_hybrid_alpha = "default";
        param.input.dft_functional = "PBE0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_hybrid_alpha, "0.25");

        param.input.exx_hybrid_alpha = "default";
        param.input.dft_functional = "SCAN0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_hybrid_alpha, "0.25");

        param.input.exx_hybrid_alpha = "default";
        param.input.dft_functional = "HSE";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_hybrid_alpha, "0.25");

        param.input.exx_hybrid_alpha = "default";
        param.input.dft_functional = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_hybrid_alpha, "0");

        param.input.exx_hybrid_alpha = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_hybrid_step
        auto it = find_label("exx_hybrid_step", readinput.input_lists);
        param.input.exx_hybrid_step = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_real_number
        auto it = find_label("exx_real_number", readinput.input_lists);
        param.input.exx_real_number = "default";
        param.input.gamma_only = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_real_number, "1");

        param.input.exx_real_number = "default";
        param.input.gamma_only = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_real_number, "0");
    }
    { // exx_ccp_rmesh_times
        auto it = find_label("exx_ccp_rmesh_times", readinput.input_lists);
        param.input.exx_ccp_rmesh_times = "default";
        param.input.dft_functional = "HF";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_ccp_rmesh_times, "5");

        param.input.exx_ccp_rmesh_times = "default";
        param.input.dft_functional = "PBE0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_ccp_rmesh_times, "5");

        param.input.exx_ccp_rmesh_times = "default";
        param.input.dft_functional = "SCAN0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_ccp_rmesh_times, "5");

        param.input.exx_ccp_rmesh_times = "default";
        param.input.dft_functional = "HSE";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_ccp_rmesh_times, "1.5");

        param.input.exx_ccp_rmesh_times = "default";
        param.input.dft_functional = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_ccp_rmesh_times, "1");

        param.input.exx_ccp_rmesh_times = "0";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_distribute_type
        auto it = find_label("exx_distribute_type", readinput.input_lists);
        param.input.exx_distribute_type = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_opt_orb_lmax
        auto it = find_label("exx_opt_orb_lmax", readinput.input_lists);
        param.input.exx_opt_orb_lmax = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_opt_orb_ecut
        auto it = find_label("exx_opt_orb_ecut", readinput.input_lists);
        param.input.exx_opt_orb_ecut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_opt_orb_tolerence
        auto it = find_label("exx_opt_orb_tolerence", readinput.input_lists);
        param.input.exx_opt_orb_tolerence = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_symmetry_realspace
        auto it = find_label("exx_symmetry_realspace", readinput.input_lists);
        param.input.exx_symmetry_realspace = true;
        param.input.symmetry="0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.exx_symmetry_realspace, false);
    }
    { // rpa_ccp_rmesh_times
        auto it = find_label("rpa_ccp_rmesh_times", readinput.input_lists);
        param.input.rpa_ccp_rmesh_times = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // berry_phase
        auto it = find_label("berry_phase", readinput.input_lists);
        param.input.berry_phase = true;
        param.input.basis_type = "pw";
        param.input.calculation = "nscf";
        param.input.gdir = 1;

        param.input.gdir = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.calculation = "scf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.basis_type = "lcao_in_pw";
        param.input.calculation = "nscf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
    }
    { // towannier90
        auto it = find_label("towannier90", readinput.input_lists);
        param.input.towannier90 = true;
        param.input.calculation = "nscf";
        param.input.nspin = 2;
        param.input.wannier_spin = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("towannier90"));

        param.input.towannier90 = true;
        param.input.calculation = "scf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("towannier90"));
    }
    { // wannier_method
        auto it = find_label("wannier_method", readinput.input_lists);
        param.input.towannier90 = true;
        param.input.basis_type = "lcao_in_pw";
        param.input.wannier_method = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.wannier_method, 1);
    }
    { // of_hold_rho0
        auto it = find_label("of_hold_rho0", readinput.input_lists);
        param.input.of_wt_rho0 = 1;
        param.input.of_hold_rho0 = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.of_hold_rho0, true);
    }
    { // of_full_pw_dim
        auto it = find_label("of_full_pw_dim", readinput.input_lists);
        param.input.of_full_pw = false;
        param.input.of_full_pw_dim = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.of_full_pw_dim, 0);
    }
    { // of_read_kernel
        auto it = find_label("of_read_kernel", readinput.input_lists);
        param.input.of_read_kernel = true;
        param.input.of_kinetic = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.of_read_kernel, false);
    }
    { // dft_plus_u
        auto it = find_label("dft_plus_u", readinput.input_lists);
        param.input.dft_plus_u = 1;
        param.input.orbital_corr = {-1, -1};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.dft_plus_u, 0);

        param.input.dft_plus_u = 1;
        param.input.basis_type = "pw";
        param.input.ks_solver = "genelpa";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.basis_type = "lcao";
        param.input.ks_solver = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // uramping
        auto it = find_label("uramping", readinput.input_lists);
        param.sys.uramping = 1;
        param.input.orbital_corr = {-1, -1};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.sys.uramping, 0);
    }
    { // onsite_radius
        auto it = find_label("onsite_radius", readinput.input_lists);
        param.input.onsite_radius = 0.0;
        param.input.dft_plus_u = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.onsite_radius, 3.0);
    }
    { // hubbard_u
        auto it = find_label("hubbard_u", readinput.input_lists);
        param.input.ntype = 2;
        it->second.str_values = {"1.0", "2.0"};
        param.sys.hubbard_u = {1.0, 2.0};
        it->second.check_value(it->second, param);
        param.input.ntype = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ntype = 2;
        param.sys.hubbard_u = {1.0, -1.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // orbital_corr
        auto it = find_label("orbital_corr", readinput.input_lists);
        param.input.ntype = 2;
        it->second.str_values = {"1", "2"};
        param.input.orbital_corr = {1, 2};
        it->second.check_value(it->second, param);
        param.input.ntype = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.ntype = 2;
        param.input.orbital_corr = {1, 4};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_nao_ecut
        auto it = find_label("bessel_nao_ecut", readinput.input_lists);
        param.input.bessel_nao_ecut = "default";
        param.input.ecutwfc = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(std::stod(param.input.bessel_nao_ecut), 1.0);

        param.input.bessel_nao_ecut = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_nao_rcut
        auto it = find_label("bessel_nao_rcut", readinput.input_lists);
        param.input.bessel_nao_rcuts = {-1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_descriptor_ecut
        auto it = find_label("bessel_descriptor_ecut", readinput.input_lists);
        param.input.bessel_descriptor_ecut = "default";
        param.input.ecutwfc = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(std::stod(param.input.bessel_descriptor_ecut), 1.0);

        param.input.bessel_descriptor_ecut = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_descriptor_rcut
        auto it = find_label("bessel_descriptor_rcut", readinput.input_lists);
        param.input.bessel_descriptor_rcut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_mag_switch
        auto it = find_label("sc_mag_switch", readinput.input_lists);
        param.input.sc_mag_switch = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_thr
        auto it = find_label("sc_thr", readinput.input_lists);
        param.input.sc_thr = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nsc
        auto it = find_label("nsc", readinput.input_lists);
        param.input.nsc = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nsc_min
        auto it = find_label("nsc_min", readinput.input_lists);
        param.input.nsc_min = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_scf_nmin
        auto it = find_label("sc_scf_nmin", readinput.input_lists);
        param.input.sc_scf_nmin = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // alpha_trial
        auto it = find_label("alpha_trial", readinput.input_lists);
        param.input.alpha_trial = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sccut
        auto it = find_label("sccut", readinput.input_lists);
        param.input.sccut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_scf_thr
        auto it = find_label("sc_scf_thr", readinput.input_lists);
        param.input.sc_scf_thr = -1e-3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // qo_thr
        auto it = find_label("qo_thr", readinput.input_lists);
        param.input.qo_thr = 1e-5;
        it->second.check_value(it->second, param);
    }
    { // qo_strategy
        auto it = find_label("qo_strategy", readinput.input_lists);
        param.input.ntype = 2;

        param.input.qo_basis = "hydrogen";
        param.input.qo_strategy = {"all"};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.qo_strategy.size(), 2);
        EXPECT_EQ(param.input.qo_strategy[0], "all");
        EXPECT_EQ(param.input.qo_strategy[1], "all");

        param.input.qo_strategy = {};
        param.input.qo_basis = "hydrogen";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.qo_strategy.size(), 2);
        EXPECT_EQ(param.input.qo_strategy[0], "energy-valence");
        EXPECT_EQ(param.input.qo_strategy[1], "energy-valence");

        param.input.qo_strategy = {};
        param.input.qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.qo_strategy.size(), 2);
        EXPECT_EQ(param.input.qo_strategy[0], "all");
        EXPECT_EQ(param.input.qo_strategy[1], "all");

        param.input.qo_basis = "test";
        param.input.qo_strategy = {};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // qo_screening_coeff
        auto it = find_label("qo_screening_coeff", readinput.input_lists);
        param.input.ntype = 2;

        it->second.str_values = {"0.2"};
        param.input.qo_screening_coeff = {0.2};
        param.input.qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.qo_screening_coeff.size(), 2);
        EXPECT_EQ(param.input.qo_screening_coeff[0], 0.2);
        EXPECT_EQ(param.input.qo_screening_coeff[1], 0.2);

        param.input.qo_screening_coeff = {};
        param.input.qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.qo_screening_coeff.size(), 2);
        EXPECT_EQ(param.input.qo_screening_coeff[0], 0.1);
        EXPECT_EQ(param.input.qo_screening_coeff[1], 0.1);

        param.input.qo_screening_coeff = {};
        param.input.qo_basis = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.qo_screening_coeff = {0.2, -0.1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.qo_screening_coeff = {0.2, 1e-8};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // md_nstep
        auto it = find_label("md_nstep", readinput.input_lists);
        param.input.mdp.md_nstep = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mdp.md_nstep, 50);
    }
    { // md_prec_level
        auto it = find_label("md_prec_level", readinput.input_lists);
        param.input.mdp.md_prec_level = 1;
        param.input.calculation = "vc-relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mdp.md_prec_level, 0);

        param.input.calculation = "md";
        param.input.mdp.md_prec_level = 1;
        param.input.mdp.md_type = "vc-md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mdp.md_prec_level, 0);
    }
    { // md_tfreq
        auto it = find_label("md_tfreq", readinput.input_lists);
        param.input.mdp.md_tfreq = 0;
        param.input.calculation = "md";
        param.input.mdp.md_dt = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mdp.md_tfreq, 1.0 / 40 / 1.0);
    }
    { // md_pfreq
        auto it = find_label("md_pfreq", readinput.input_lists);
        param.input.mdp.md_pfreq = 0;
        param.input.calculation = "md";
        param.input.mdp.md_dt = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.mdp.md_pfreq, 1.0 / 400 / 1.0);
    }
    { // lj_rule
        auto it = find_label("lj_rule", readinput.input_lists);
        param.input.esolver_type = "lj";
        param.input.mdp.lj_rule = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_rcut
        auto it = find_label("lj_rcut", readinput.input_lists);
        param.input.ntype = 2;
        param.input.esolver_type = "lj";
        it->second.str_values = {"1.0", "2.0"};
        param.input.mdp.lj_rcut = {1.0, 2.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        param.input.mdp.lj_rcut = {1.0, 2.0, -1.0};
        it->second.str_values = {"1.0", "2.0", "-1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_epsilon
        auto it = find_label("lj_epsilon", readinput.input_lists);
        param.input.ntype = 2;
        param.input.esolver_type = "lj";
        param.input.mdp.lj_epsilon = {1.0};
        it->second.str_values = {"1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_sigma
        auto it = find_label("lj_sigma", readinput.input_lists);
        param.input.ntype = 2;
        param.input.esolver_type = "lj";
        param.input.mdp.lj_sigma = {1.0};
        it->second.str_values = {"1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(0), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nocc 
        auto it = find_label("nocc", readinput.input_lists);
        param.input.nocc = 5;
        param.input.nbands = 4;
        param.input.nelec = 0.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nocc, 4);
        param.input.nocc = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(param.input.nocc, 4);
    }
}
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
        std::vector<std::pair<std::string, ModuleIO::Input_Item>>& input_lists)
    {
        auto it = std::find_if(
            input_lists.begin(),
            input_lists.end(),
            [&label](const std::pair<std::string, ModuleIO::Input_Item>& item) { return item.first == label; });
        return it;
    }
};

TEST_F(InputTest, RelaxMethod)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    auto it = find_label("relax_method", readinput.input_lists);

    it->second.str_values = {"cg"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"cg", "2"}));
    EXPECT_TRUE(TestParameters::input(param).uses_simultaneous_relaxation());

    it->second.str_values = {"cg", "1"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"cg", "1"}));
    EXPECT_FALSE(TestParameters::input(param).uses_simultaneous_relaxation());

    it->second.str_values = {"cg", "2"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"cg", "2"}));
    EXPECT_TRUE(TestParameters::input(param).uses_simultaneous_relaxation());

    it->second.str_values = {"bfgs"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"bfgs", "2"}));
    EXPECT_FALSE(TestParameters::input(param).uses_simultaneous_relaxation());

    it->second.str_values = {"bfgs", "1"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"bfgs", "1"}));

    it->second.str_values = {"bfgs", "2"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(TestParameters::input(param).relax_method, (std::vector<std::string>{"bfgs", "2"}));

    for (const std::vector<std::string>& invalid : {
             std::vector<std::string>{"cg", "3"},
             std::vector<std::string>{"bfgs", "3"},
             std::vector<std::string>{"sd", "1"},
             std::vector<std::string>{"none"},
             std::vector<std::string>{"cg", "2", "extra"}})
    {
        it->second.str_values = invalid;
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
    }

    EXPECT_EQ(find_label("relax_new", readinput.input_lists), readinput.input_lists.end());
}

TEST_F(InputTest, Item_test)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;

    std::string output = "";

    { // calculation
        auto it = find_label("calculation", readinput.input_lists);
        TestParameters::input(param).calculation = "get_pchg";
        TestParameters::input(param).basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).calculation = "gen_bessel";
        TestParameters::input(param).basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).calculation = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }

    { // esolver_type
        auto it = find_label("esolver_type", readinput.input_lists);
        TestParameters::input(param).esolver_type = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).esolver_type = "dp";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).esolver_type = "lr";
        TestParameters::input(param).calculation = "scf";
        it = find_label("esolver_type", readinput.input_lists);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("esolver_type=lr requires calculation=nscf"));
    }
    { // nspin
        auto it = find_label("nspin", readinput.input_lists);
        TestParameters::input(param).noncolin = false;
        TestParameters::input(param).nspin = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // smearing_method
        auto it = find_label("smearing_method", readinput.input_lists);
        TestParameters::input(param).smearing_method = "fix";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // kspacing
        auto it = find_label("kspacing", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).kspacing[0], 1);
        EXPECT_EQ(TestParameters::input(param).kspacing[1], 1);
        EXPECT_EQ(TestParameters::input(param).kspacing[2], 1);
        it->second.str_values = {"1", "2"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).kspacing = {0, -1, 1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).kspacing = {0, 1, 2};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nbands
        auto it = find_label("nbands", readinput.input_lists);
        TestParameters::input(param).nbands = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // symmetry
        auto it = find_label("symmetry", readinput.input_lists);
        TestParameters::input(param).symmetry = "default";
        TestParameters::input(param).gamma_only = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "0");

        TestParameters::input(param).symmetry = "default";
        TestParameters::input(param).gamma_only = false;
        TestParameters::input(param).calculation = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "1");

        TestParameters::input(param).calculation = "md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "0");

        TestParameters::input(param).symmetry = "default";
        TestParameters::input(param).efield_flag = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "0");

        TestParameters::input(param).qo_switch = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "-1");

        TestParameters::input(param).symmetry = "default";
        TestParameters::input(param).berry_phase = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).symmetry, "-1");
    }
    { // nelec
        auto it = find_label("nelec", readinput.input_lists);
        TestParameters::input(param).nelec = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).nelec = 100;
        TestParameters::input(param).nbands = 5;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bndpar
        auto it = find_label("bndpar", readinput.input_lists);
        TestParameters::input(param).esolver_type = "sdft";
        TestParameters::input(param).bndpar = 1;
        EXPECT_NO_THROW(it->second.check_value(it->second, param));
    }
    { // dft_plus_dmft
        auto it = find_label("dft_plus_dmft", readinput.input_lists);
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).dft_plus_dmft = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // mem_saver
        auto it = find_label("mem_saver", readinput.input_lists);
        TestParameters::input(param).mem_saver = 1;
        TestParameters::input(param).calculation = "scf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mem_saver, 0);

        TestParameters::input(param).mem_saver = 1;
        TestParameters::input(param).calculation = "relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mem_saver, 0);
    }
    { // gint_precision
        auto it = find_label("gint_precision", readinput.input_lists);
        TestParameters::input(param).gint_precision = "invalid";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).gint_precision = "mix";
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).device = "cpu";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // diag_proc
        auto it = find_label("diago_proc", readinput.input_lists);
        TestParameters::input(param).diago_proc = 0;
        GlobalV::NPROC = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).diago_proc, 1);
    }
    { // cal_force
        auto it = find_label("cal_force", readinput.input_lists);
        TestParameters::input(param).calculation = "cell-relax";
        TestParameters::input(param).cal_force = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).cal_force, true);

        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).cal_force = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).cal_force, false);
    }
    { // ecutrho
        auto it = find_label("ecutrho", readinput.input_lists);
        TestParameters::input(param).ecutwfc = 1;
        TestParameters::input(param).ecutrho = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ecutrho, 4);
        TestParameters::input(param).nx = 0;
        TestParameters::input(param).ecutrho = 5;
        TestParameters::input(param).ecutwfc = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::sys(param).double_grid, true);

        TestParameters::input(param).ecutwfc = 1;
        TestParameters::input(param).ecutrho = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // pw_diag_thr
        auto it = find_label("pw_diag_thr", readinput.input_lists);
        TestParameters::input(param).pw_diag_thr = 1.0e-2;
        TestParameters::input(param).calculation = "get_s";
        TestParameters::input(param).basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).pw_diag_thr, 1.0e-5);
    }
    { // nb2d
        auto it = find_label("nb2d", readinput.input_lists);
        TestParameters::input(param).nb2d = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // scf_thr
        auto it = find_label("scf_thr", readinput.input_lists);
        TestParameters::input(param).scf_thr = -1;
        TestParameters::input(param).basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).scf_thr, 1.0e-7);

        TestParameters::input(param).scf_thr = -1;
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).calculation = "nscf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).scf_thr, 1.0e-6);

        TestParameters::input(param).scf_thr = -1;
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).calculation = "scf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).scf_thr, 1.0e-9);
    }
    { // scf_thr_type
        auto it = find_label("scf_thr_type", readinput.input_lists);
        TestParameters::input(param).scf_thr_type = -1;
        TestParameters::input(param).basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).scf_thr_type, 2);

        TestParameters::input(param).scf_thr_type = -1;
        TestParameters::input(param).basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).scf_thr_type, 1);
    }
    { // init_wfc
        auto it = find_label("init_wfc", readinput.input_lists);
        TestParameters::input(param).init_wfc = "atomic";
        TestParameters::input(param).calculation = "get_pchg";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).init_wfc, "file");

        TestParameters::input(param).init_wfc = "atomic";
        TestParameters::input(param).basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).init_wfc, "nao");
    }
    { // init_chg
        auto it = find_label("init_chg", readinput.input_lists);
        TestParameters::input(param).init_chg = "get_pchg";
        TestParameters::input(param).init_chg = "";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).init_chg, "atomic");

        TestParameters::input(param).init_chg = "";
        TestParameters::input(param).calculation = "nscf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).init_chg, "file");

        TestParameters::input(param).init_chg = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // chg_extrap
        auto it = find_label("chg_extrap", readinput.input_lists);
        TestParameters::input(param).chg_extrap = "default";
        TestParameters::input(param).calculation = "md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).chg_extrap, "second-order");

        TestParameters::input(param).chg_extrap = "default";
        TestParameters::input(param).calculation = "relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).chg_extrap, "first-order");

        TestParameters::input(param).chg_extrap = "default";
        TestParameters::input(param).calculation = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).chg_extrap, "atomic");

        TestParameters::input(param).chg_extrap = "none";
        TestParameters::input(param).calculation = "get_wf";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).chg_extrap, "atomic");
    }
    { // out_chg
        auto it = find_label("out_chg", readinput.input_lists);
        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_chg = {0};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_chg[0], 1);
    }
    { // out_pot
        auto it = find_label("out_pot", readinput.input_lists);
        TestParameters::input(param).calculation = "scf";
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_pot[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_pot[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_pot[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_pot[1], 2);
    }
    { // out_dos
        auto it = find_label("out_dos", readinput.input_lists);
        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_dos = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dos, 0);

        TestParameters::input(param).out_dos = 3;
        TestParameters::input(param).symmetry = "1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).out_dos = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // out_ldos
        auto it = find_label("out_ldos", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_ldos[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_ldos[1], 3);

        it->second.str_values = {"2", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_ldos[0], 2);
        EXPECT_EQ(TestParameters::input(param).out_ldos[1], 2);

        it->second.str_values = {"1", "2", "3"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).out_ldos = {5, 5};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // stm_bias
        auto it = find_label("stm_bias", readinput.input_lists);
        it->second.str_values = {"2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).stm_bias[0], 2);
        EXPECT_EQ(TestParameters::input(param).stm_bias[1], 0.1);
        EXPECT_EQ(TestParameters::input(param).stm_bias[2], 1);

        it->second.str_values = {"3", "0.2", "5"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).stm_bias[0], 3);
        EXPECT_EQ(TestParameters::input(param).stm_bias[1], 0.2);
        EXPECT_EQ(TestParameters::input(param).stm_bias[2], 5);

        it->second.str_values = {"1", "2"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).stm_bias = {1.0, 0.1, 0.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).stm_bias = {1.0, 0.0, 2.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ldos_line
        auto it = find_label("ldos_line", readinput.input_lists);
        it->second.str_values = {"1", "2", "3", "4", "5", "6"};
        it->second.read_value(it->second, param);
        for (int i = 0; i < 6; ++i)
        {
            EXPECT_EQ(TestParameters::input(param).ldos_line[i], i + 1);
        }
        EXPECT_EQ(TestParameters::input(param).ldos_line[6], 100);

        it->second.str_values = {"2", "3", "4", "5", "6", "7", "200"};
        it->second.read_value(it->second, param);
        for (int i = 0; i < 6; ++i)
        {
            EXPECT_EQ(TestParameters::input(param).ldos_line[i], i + 2);
        }
        EXPECT_EQ(TestParameters::input(param).ldos_line[6], 200);

        it->second.str_values = {"1", "2", "3", "4", "5", "6", "7", "8"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ldos_line = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // out_band
        auto it = find_label("out_band", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_band[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_band[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_band[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_band[1], 2);

        it->second.str_values = {"1", "2", "3"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_band = {1, 2};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_band[0], 0);
    }
    { // out_proj_band
        auto it = find_label("out_proj_band", readinput.input_lists);
        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_proj_band = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_proj_band, false);

        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).out_proj_band = true;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // read_file_dir
        auto it = find_label("read_file_dir", readinput.input_lists);
        TestParameters::input(param).read_file_dir = "auto";
        TestParameters::input(param).suffix = "test";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).read_file_dir, "OUT.test/");

        TestParameters::input(param).read_file_dir = "test";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).read_file_dir, "test/");
    }
    { // nx
        auto it = find_label("nx", readinput.input_lists);
        TestParameters::input(param).nx = 1;
        TestParameters::input(param).ny = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ny
        auto it = find_label("ny", readinput.input_lists);
        TestParameters::input(param).ny = 1;
        TestParameters::input(param).nz = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nz
        auto it = find_label("nz", readinput.input_lists);
        TestParameters::input(param).nz = 1;
        TestParameters::input(param).nx = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndx
        auto it = find_label("ndx", readinput.input_lists);
        TestParameters::sys(param).double_grid = false;
        TestParameters::input(param).ndx = 2;
        TestParameters::input(param).nx = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::sys(param).double_grid, true);

        TestParameters::input(param).ndx = 1;
        TestParameters::input(param).ndy = 0;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ndx = 1;
        TestParameters::input(param).nx = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndy
        auto it = find_label("ndy", readinput.input_lists);
        TestParameters::sys(param).double_grid = false;
        TestParameters::input(param).ndy = 2;
        TestParameters::input(param).ny = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::sys(param).double_grid, true);

        TestParameters::input(param).ndy = 1;
        TestParameters::input(param).ndz = 0;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ndy = 1;
        TestParameters::input(param).ny = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // ndz
        auto it = find_label("ndz", readinput.input_lists);
        TestParameters::sys(param).double_grid = false;
        TestParameters::input(param).ndz = 2;
        TestParameters::input(param).nz = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::sys(param).double_grid, true);

        TestParameters::input(param).ndz = 1;
        TestParameters::input(param).nz = 2;
        it->second.str_values = {"1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ndz = 1;
        TestParameters::input(param).nz = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // cell_factor
        auto it = find_label("cell_factor", readinput.input_lists);
        TestParameters::input(param).calculation = "cell-relax";
        TestParameters::input(param).cell_factor = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).cell_factor, 2.0);
    }
    { // ks_sovler
        auto it = find_label("ks_solver", readinput.input_lists);
        TestParameters::input(param).ks_solver = "default";
        TestParameters::input(param).basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ks_solver, "cg");

        TestParameters::input(param).ks_solver = "default";
        TestParameters::input(param).basis_type = "lcao";
        TestParameters::input(param).device = "gpu";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ks_solver, "cusolver");

        TestParameters::input(param).ks_solver = "genelpa";
        TestParameters::input(param).basis_type = "lcao";
        TestParameters::input(param).device = "gpu";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("Please use ks_solver = elpa with device = gpu"));
#ifdef __ELPA
        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ks_solver, "genelpa");
#else
#ifdef __MPI
        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ks_solver, "scalapack_gvx");
#else
        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).basis_type = "lcao_in_pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).ks_solver, "lapack");
#endif
#endif
        TestParameters::input(param).ks_solver = "default";
        TestParameters::input(param).basis_type = "lcao";
        TestParameters::input(param).device = "cpu";
        TestParameters::input(param).ks_solver = "lapack";

        TestParameters::input(param).ks_solver = "genelpa";
        TestParameters::input(param).basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ks_solver = "cg";
        TestParameters::input(param).basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ks_solver = "scalapack_gvx";
        TestParameters::input(param).basis_type = "lcao";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ks_solver = "cg";
        TestParameters::input(param).basis_type = "lcao_in_pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // relax_nmax
        auto it = find_label("relax_nmax", readinput.input_lists);
        TestParameters::input(param).calculation = "scf";
        TestParameters::input(param).relax_nmax = 5;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).relax_nmax, 1);

        TestParameters::input(param).relax_nmax = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).relax_nmax, 0);

        TestParameters::input(param).calculation = "relax";
        TestParameters::input(param).relax_nmax = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).relax_nmax, 0);

        TestParameters::input(param).relax_nmax = -1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).relax_nmax, 50);
    }
    { // out_stru
        auto it = find_label("out_stru", readinput.input_lists);
        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_stru = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }
    { // cal_stress
        auto it = find_label("cal_stress", readinput.input_lists);
        TestParameters::input(param).calculation = "md";
        TestParameters::input(param).cal_stress = false;
        TestParameters::input(param).esolver_type = "lj";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).cal_stress, true);
    }
    { // fixed_axes
        auto it = find_label("fixed_axes", readinput.input_lists);
        TestParameters::input(param).fixed_axes = "shape";
        TestParameters::input(param).relax_method = {"cg", "1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).fixed_axes = "volume";
        TestParameters::input(param).relax_method = {"cg", "1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // fixed_ibrav
        auto it = find_label("fixed_ibrav", readinput.input_lists);
        TestParameters::input(param).fixed_ibrav = true;
        TestParameters::input(param).relax_method = {"cg", "1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).fixed_ibrav = true;
        TestParameters::input(param).latname = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // fixed_atoms
        auto it = find_label("fixed_atoms", readinput.input_lists);
        TestParameters::input(param).fixed_atoms = true;
        TestParameters::input(param).calculation = "relax";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // force_thr
        auto it = find_label("force_thr", readinput.input_lists);
        TestParameters::input(param).force_thr = -1;
        TestParameters::input(param).force_thr_ev = -1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).force_thr, 1.0e-3);
        EXPECT_EQ(TestParameters::input(param).force_thr_ev, 1.0e-3 * 13.6058 / 0.529177);

        TestParameters::input(param).force_thr = -1;
        TestParameters::input(param).force_thr_ev = 1.0e-3 * 13.6058 / 0.529177;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).force_thr, 1.0e-3);

        TestParameters::input(param).force_thr = 1.0e-3;
        TestParameters::input(param).force_thr_ev = 1.0e-3 * 13.6058 / 0.529177;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).force_thr, 1.0e-3);
        EXPECT_EQ(TestParameters::input(param).force_thr_ev, 1.0e-3 * 13.6058 / 0.529177);
    }
    { // out_level
        auto it = find_label("out_level", readinput.input_lists);
        TestParameters::input(param).out_level = "0";
        TestParameters::input(param).calculation = "md";
        TestParameters::sys(param).out_md_control = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_level, "m");
    }
    { // out_dmk
        auto it = find_label("out_dmk", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dmk[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_dmk[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dmk[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_dmk[1], 2);
    }
    { // out_dmr
        auto it = find_label("out_dmr", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dmr[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_dmr[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dmr[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_dmr[1], 2);
    }
    { // method_sto
        auto it = find_label("method_sto", readinput.input_lists);
        TestParameters::input(param).method_sto = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nbands_sto
        auto it = find_label("nbands_sto", readinput.input_lists);
        TestParameters::input(param).esolver_type = "sdft";

        it->second.str_values = {"all"};
        it->second.read_value(it->second, param);
        EXPECT_NO_THROW(it->second.check_value(it->second, param));
        EXPECT_EQ(TestParameters::input(param).nbands_sto, 0);
        EXPECT_EQ(TestParameters::input(param).esolver_type, "sdft");

        it->second.str_values = {"8"};
        it->second.read_value(it->second, param);
        EXPECT_NO_THROW(it->second.check_value(it->second, param));
        EXPECT_EQ(TestParameters::input(param).nbands_sto, 8);
        EXPECT_EQ(TestParameters::input(param).esolver_type, "sdft");

        it->second.str_values = {"1000000"};
        it->second.read_value(it->second, param);
        EXPECT_NO_THROW(it->second.check_value(it->second, param));

        it->second.str_values = {"1000001"};
        it->second.read_value(it->second, param);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        it->second.str_values = {"0"};
        it->second.read_value(it->second, param);
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
        EXPECT_EQ(TestParameters::input(param).nbands_sto, 0);
        EXPECT_EQ(TestParameters::input(param).esolver_type, "sdft");

        it->second.str_values = {"-1"};
        TestParameters::input(param).nbands_sto = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        it->second.str_values = {"all"};
        it->second.get_final_value(it->second, param);
        EXPECT_EQ(it->second.final_value.str(), "all");
        it->second.final_value.str("");

        it->second.str_values = {};
        TestParameters::input(param).nbands_sto = 256;
        it->second.get_final_value(it->second, param);
        EXPECT_EQ(it->second.final_value.str(), "256");
    }
    { // basis_type
        auto it = find_label("basis_type", readinput.input_lists);
        TestParameters::input(param).basis_type = "lcao_in_pw";
        TestParameters::input(param).towannier90 = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).basis_type, "lcao");

        TestParameters::input(param).basis_type = "gauss";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // gamma_only
        auto it = find_label("gamma_only", readinput.input_lists);
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).gamma_only = true;
        std::string filename = TestParameters::input(param).kpoint_file;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).gamma_only, false);
        std::ifstream ifs(filename);
        std::string line;
        std::getline(ifs, line);
        ifs.close();
        EXPECT_EQ(line, "K_POINTS");

        TestParameters::input(param).basis_type = "lcao";
        TestParameters::input(param).gamma_only = true;
        TestParameters::input(param).nspin = 4;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

    }
    { // out_mat_r
        auto it = find_label("out_mat_r", readinput.input_lists);
        TestParameters::input(param).out_hsr[0] = 1;
        TestParameters::sys(param).gamma_only_local = true;
        it->second.check_value(it->second, param);
        TestParameters::input(param).out_hsr[0] = 3;
        it->second.check_value(it->second, param);
        TestParameters::input(param).out_hsr[0] = 0;

        TestParameters::input(param).esolver_type = "lcao";
        TestParameters::input(param).out_mat_r[0] = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("available"));
        TestParameters::input(param).out_mat_r[0] = 0;
        TestParameters::sys(param).gamma_only_local = false;
    }
    { // lcao_ecut
        auto it = find_label("lcao_ecut", readinput.input_lists);
        TestParameters::input(param).lcao_ecut = 0;
        TestParameters::input(param).ecutwfc = 1;
        TestParameters::input(param).basis_type = "lcao";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).lcao_ecut, 1);
    }
    { // out_mat_hs
        auto it = find_label("out_mat_hs", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_hs[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_hs[1], 8);

        it->second.str_values = {"1", "2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_hs[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_hs[1], 2);

    }
    { // out_hsk
        auto it = find_label("out_hsk", readinput.input_lists);
        it->second.str_values = {"1", "12"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hsk[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_hsk[1], 12);

        it->second.str_values = {"2", "12"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hsk[0], 2);
        EXPECT_EQ(TestParameters::input(param).out_hsk[1], 12);
        it->second.check_value(it->second, param);

        TestParameters::input(param).out_hsk[0] = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NPZ output is not implemented"));
    }
    { // out_hsr
        auto it = find_label("out_hsr", readinput.input_lists);
        it->second.str_values = {"1", "10"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hsr[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_hsr[1], 10);

        it->second.str_values = {"2", "12"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hsr[0], 2);
        EXPECT_EQ(TestParameters::input(param).out_hsr[1], 12);
        it->second.check_value(it->second, param);

        TestParameters::input(param).out_hsr[0] = 4;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("format must be 0, 1, 2, or 3"));

#ifndef __CNPY
        TestParameters::input(param).out_hsr[0] = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("ENABLE_CNPY"));
#else
        TestParameters::input(param).out_hsr[0] = 3;
        it->second.check_value(it->second, param);
#endif
    }
    { // out_hr_npz
        auto it = find_label("out_hr_npz", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hr_npz, true);
    }
    { // out_hsr_npz
        auto it = find_label("out_hsr_npz", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_hsr_npz, true);
    }
    { // out_dm_npz
        auto it = find_label("out_dm_npz", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_dm_npz, true);
    }
}

TEST_F(InputTest, HsOutputAliases)
{
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto legacy = find_label("out_mat_hs", readinput.input_lists);
        auto primary = find_label("out_hsk", readinput.input_lists);
        legacy->second.str_values = {"1", "5"};
        primary->second.str_values = {"0"};
        legacy->second.read_value(legacy->second, param);
        primary->second.read_value(primary->second, param);

        readinput.normalize_hs_output_options(param);
        EXPECT_EQ(TestParameters::input(param).out_hsk[0], 0);
        EXPECT_EQ(TestParameters::input(param).out_hsk[1], 8);
    }
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto legacy = find_label("out_mat_hs", readinput.input_lists);
        auto primary = find_label("out_hsk", readinput.input_lists);
        legacy->second.str_values = {"1", "5"};
        primary->second.str_values = {"0"};
        primary->second.read_value(primary->second, param);
        legacy->second.read_value(legacy->second, param);

        readinput.normalize_hs_output_options(param);
        EXPECT_EQ(TestParameters::input(param).out_hsk[0], 0);
        EXPECT_EQ(TestParameters::input(param).out_hsk[1], 8);
    }
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto legacy_text = find_label("out_mat_hs2", readinput.input_lists);
        auto legacy_npz = find_label("out_hsr_npz", readinput.input_lists);
        legacy_text->second.str_values = {"1", "5"};
        legacy_npz->second.str_values = {"1"};
        legacy_text->second.read_value(legacy_text->second, param);
        legacy_npz->second.read_value(legacy_npz->second, param);
        readinput.normalize_hs_output_options(param);

        EXPECT_EQ(TestParameters::input(param).out_hsr[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_hsr[1], 5);
        EXPECT_TRUE(TestParameters::input(param).out_hsr_npz);
        EXPECT_TRUE(TestParameters::input(param).out_hsr_npz_compat);
    }
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto legacy_text = find_label("out_mat_hs2", readinput.input_lists);
        auto legacy_npz = find_label("out_hsr_npz", readinput.input_lists);
        auto primary = find_label("out_hsr", readinput.input_lists);
        legacy_text->second.str_values = {"1", "5"};
        legacy_npz->second.str_values = {"1"};
        primary->second.str_values = {"1", "12"};
        legacy_text->second.read_value(legacy_text->second, param);
        legacy_npz->second.read_value(legacy_npz->second, param);
        primary->second.read_value(primary->second, param);

        readinput.normalize_hs_output_options(param);
        EXPECT_EQ(TestParameters::input(param).out_hsr[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_hsr[1], 12);
        EXPECT_FALSE(TestParameters::input(param).out_hsr_npz);
        EXPECT_FALSE(TestParameters::input(param).out_hsr_npz_compat);
    }
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto legacy_npz = find_label("out_hsr_npz", readinput.input_lists);
        legacy_npz->second.str_values = {"1"};
        legacy_npz->second.read_value(legacy_npz->second, param);
        readinput.normalize_hs_output_options(param);
        EXPECT_EQ(TestParameters::input(param).out_hsr[0], 3);
        EXPECT_EQ(TestParameters::input(param).out_hsr[1], 8);
    }
    {
        ModuleIO::ReadInput readinput(0);
        Parameter param;
        auto primary = find_label("out_hsk", readinput.input_lists);
        primary->second.str_values = {"0"};
        primary->second.read_value(primary->second, param);
        TestParameters::input(param).qo_switch = true;
        readinput.normalize_hs_output_options(param);
        EXPECT_EQ(TestParameters::input(param).out_hsk[0], 1);
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
        TestParameters::input(param).out_mat_dh[0] = 1;
        TestParameters::input(param).nspin = 4;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // dm_to_rho
        auto it = find_label("dm_to_rho", readinput.input_lists);
        TestParameters::input(param).dm_to_rho = true;
        GlobalV::NPROC = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).dm_to_rho = true;
        TestParameters::input(param).gamma_only = true;
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

#ifndef __CNPY
        TestParameters::input(param).dm_to_rho = true;
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
#endif
    }
    { // out_wfc_lcao
        auto it = find_label("out_wfc_lcao", readinput.input_lists);
        TestParameters::input(param).out_wfc_lcao = false;
        TestParameters::input(param).qo_switch = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_wfc_lcao, true);

        TestParameters::input(param).out_wfc_lcao = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).out_wfc_lcao = 1;
        TestParameters::input(param).basis_type = "pw";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bx
        auto it = find_label("bx", readinput.input_lists);
        TestParameters::input(param).bx = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).bx = 2;
        TestParameters::input(param).basis_type = "pw";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).bx, 1);
        EXPECT_EQ(TestParameters::input(param).by, 1);
        EXPECT_EQ(TestParameters::input(param).bz, 1);
    }
    { // by
        auto it = find_label("by", readinput.input_lists);
        TestParameters::input(param).by = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bz
        auto it = find_label("bz", readinput.input_lists);
        TestParameters::input(param).bz = 11;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // mixing_beta
        auto it = find_label("mixing_beta", readinput.input_lists);
        TestParameters::input(param).mixing_beta = -1;
        TestParameters::input(param).nspin = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mixing_beta, 0.8);

        TestParameters::input(param).mixing_beta = -1;
        TestParameters::input(param).nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mixing_beta, 0.4);
        EXPECT_EQ(TestParameters::input(param).mixing_beta_mag, 1.6);
        EXPECT_EQ(TestParameters::input(param).mixing_gg0_mag, 0.0);

        TestParameters::input(param).mixing_beta = -1;
        TestParameters::input(param).nspin = 4;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mixing_beta, 0.4);
        EXPECT_EQ(TestParameters::input(param).mixing_beta_mag, 1.6);
        EXPECT_EQ(TestParameters::input(param).mixing_gg0_mag, 0.0);
    }
    { // mixing_beta_mag
        auto it = find_label("mixing_beta_mag", readinput.input_lists);
        TestParameters::input(param).mixing_beta = 0.3;
        TestParameters::input(param).mixing_beta_mag = -1;
        TestParameters::input(param).nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mixing_beta_mag, 1.2);

        TestParameters::input(param).mixing_beta = 0.5;
        TestParameters::input(param).mixing_beta_mag = -1;
        TestParameters::input(param).nspin = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mixing_beta_mag, 1.6);
    }
    { // dip_cor_flag
        auto it = find_label("dip_cor_flag", readinput.input_lists);
        TestParameters::input(param).dip_cor_flag = true;
        TestParameters::input(param).efield_flag = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // efield_dir
        auto it = find_label("efield_dir", readinput.input_lists);
        TestParameters::input(param).gate_flag = true;
        TestParameters::input(param).efield_flag = true;
        TestParameters::input(param).dip_cor_flag = false;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_s6
        auto it = find_label("vdw_s6", readinput.input_lists);
        TestParameters::input(param).vdw_s6 = "default";
        TestParameters::input(param).vdw_method = "d2";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_s6, "0.75");

        // dftd3 parameter will not get its value here
        TestParameters::input(param).vdw_s6 = "default";
        TestParameters::input(param).vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s6, "1.0");
        EXPECT_EQ(TestParameters::input(param).vdw_s6, "default");
    }
    { // vdw_s8
        auto it = find_label("vdw_s8", readinput.input_lists);
        TestParameters::input(param).vdw_s8 = "default";
        TestParameters::input(param).vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s8, "0.722");
        EXPECT_EQ(TestParameters::input(param).vdw_s8, "default");

        TestParameters::input(param).vdw_s8 = "default";
        TestParameters::input(param).vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_s8, "0.7875");
        EXPECT_EQ(TestParameters::input(param).vdw_s8, "default");
    }
    { // vdw_a1
        auto it = find_label("vdw_a1", readinput.input_lists);
        TestParameters::input(param).vdw_a1 = "default";
        TestParameters::input(param).vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a1, "1.217");
        EXPECT_EQ(TestParameters::input(param).vdw_a1, "default");

        TestParameters::input(param).vdw_a1 = "default";
        TestParameters::input(param).vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a1, "0.4289");
        EXPECT_EQ(TestParameters::input(param).vdw_a1, "default");
    }
    { // vdw_a2
        auto it = find_label("vdw_a2", readinput.input_lists);
        TestParameters::input(param).vdw_a2 = "default";
        TestParameters::input(param).vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a2, "1.0");
        EXPECT_EQ(TestParameters::input(param).vdw_a2, "default");

        TestParameters::input(param).vdw_a2 = "default";
        TestParameters::input(param).vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        // EXPECT_EQ(param.input.vdw_a2, "4.4407");
        EXPECT_EQ(TestParameters::input(param).vdw_a2, "default");
    }
    { // vdw_c6_unit
        auto it = find_label("vdw_c6_unit", readinput.input_lists);
        TestParameters::input(param).vdw_C6_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_r0_unit
        auto it = find_label("vdw_r0_unit", readinput.input_lists);
        TestParameters::input(param).vdw_R0_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_type
        auto it = find_label("vdw_cutoff_type", readinput.input_lists);
        TestParameters::input(param).vdw_cutoff_type = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_radius
        auto it = find_label("vdw_cutoff_radius", readinput.input_lists);
        TestParameters::input(param).vdw_cutoff_radius = "default";
        TestParameters::input(param).vdw_method = "d2";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_radius, "56.6918");

        TestParameters::input(param).vdw_cutoff_radius = "default";
        TestParameters::input(param).vdw_method = "d3_0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_radius, "95");

        TestParameters::input(param).vdw_cutoff_radius = "default";
        TestParameters::input(param).vdw_method = "d3_bj";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_radius, "95");

        TestParameters::input(param).vdw_cutoff_radius = "default";
        TestParameters::input(param).vdw_method = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_radius, "0");

        TestParameters::input(param).vdw_cutoff_radius = "-1";
        TestParameters::input(param).vdw_method = "d2";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_radius_unit
        auto it = find_label("vdw_radius_unit", readinput.input_lists);
        TestParameters::input(param).vdw_radius_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_width2
        auto it = find_label("vdw_cutoff_width2", readinput.input_lists);
        TestParameters::input(param).vdw_cutoff_width2 = -1.0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_width3
        auto it = find_label("vdw_cutoff_width3", readinput.input_lists);
        TestParameters::input(param).vdw_cutoff_width3 = -1.0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cn_thr
        auto it = find_label("vdw_cn_thr", readinput.input_lists);
        TestParameters::input(param).vdw_cn_thr = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cn_thr_unit
        auto it = find_label("vdw_cn_thr_unit", readinput.input_lists);
        TestParameters::input(param).vdw_cn_thr_unit = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // vdw_cutoff_period
        auto it = find_label("vdw_cutoff_period", readinput.input_lists);
        it->second.str_values = {"1", "1", "1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_period[0], 1);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_period[1], 1);
        EXPECT_EQ(TestParameters::input(param).vdw_cutoff_period[2], 1);

        it->second.str_values = {"1", "1"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).vdw_cutoff_period = {-1, 1, 1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_fock_alpha
        auto it = find_label("exx_fock_alpha", readinput.input_lists);
        TestParameters::input(param).exx_fock_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "HF";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_fock_alpha[0], "1");

        TestParameters::input(param).exx_fock_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "PBE0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_fock_alpha[0], "0.25");

        TestParameters::input(param).exx_fock_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "SCAN0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_fock_alpha[0], "0.25");

        TestParameters::input(param).exx_fock_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "cam_pbeh";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_fock_alpha[0], "0.2");
    }
    { // exx_erfc_alpha
        auto it = find_label("exx_erfc_alpha", readinput.input_lists);
        TestParameters::input(param).exx_erfc_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "HSE";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_alpha[0], "0.25");

        TestParameters::input(param).exx_erfc_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "lc_wpbe";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_alpha[0], "-1");

        TestParameters::input(param).exx_erfc_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "lrc_wpbe";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_alpha[0], "-1");

        TestParameters::input(param).exx_erfc_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "lrc_wpbeh";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_alpha[0], "-0.8");

        TestParameters::input(param).exx_erfc_alpha[0] = "default";
        TestParameters::input(param).dft_functional = "cam_pbeh";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_alpha[0], "0.8");
    }
    { // exx_erfc_omega
        auto it = find_label("exx_erfc_omega", readinput.input_lists);
        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "lc_pbe";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.33");

        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "lc_wpbe";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.4");

        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "lrc_wpbe";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.3");

        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "lrc_wpbeh";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.2");

        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "cam_pbeh";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.7");

        TestParameters::input(param).exx_erfc_omega[0] = "default";
        TestParameters::input(param).dft_functional = "hse";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_erfc_omega[0], "0.11");
    }
    { // exx_hybrid_step
        auto it = find_label("exx_hybrid_step", readinput.input_lists);
        TestParameters::input(param).exx_hybrid_step = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_real_number
        auto it = find_label("exx_real_number", readinput.input_lists);
        TestParameters::input(param).exx_real_number = "default";
        TestParameters::input(param).gamma_only = true;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_real_number, "1");

        TestParameters::input(param).exx_real_number = "default";
        TestParameters::input(param).gamma_only = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_real_number, "0");
    }
    { // exx_singluarity_correction
        auto it = find_label("exx_singularity_correction", readinput.input_lists);
        TestParameters::input(param).exx_singularity_correction = "default";
        TestParameters::input(param).dft_functional = "HF";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_singularity_correction, "spencer");

        TestParameters::input(param).exx_singularity_correction = "default";
        TestParameters::input(param).dft_functional = "PBE0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_singularity_correction, "spencer");

        TestParameters::input(param).exx_singularity_correction = "default";
        TestParameters::input(param).dft_functional = "HSE";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_singularity_correction, "limits");
    }
    { // exx_ccp_rmesh_times
        auto it = find_label("exx_ccp_rmesh_times", readinput.input_lists);
        TestParameters::input(param).exx_ccp_rmesh_times = "default";
        TestParameters::input(param).dft_functional = "HF";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_ccp_rmesh_times, "5");

        TestParameters::input(param).exx_ccp_rmesh_times = "default";
        TestParameters::input(param).dft_functional = "PBE0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_ccp_rmesh_times, "5");

        TestParameters::input(param).exx_ccp_rmesh_times = "default";
        TestParameters::input(param).dft_functional = "SCAN0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_ccp_rmesh_times, "5");

        TestParameters::input(param).exx_ccp_rmesh_times = "default";
        TestParameters::input(param).dft_functional = "HSE";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_ccp_rmesh_times, "1.5");

        TestParameters::input(param).exx_ccp_rmesh_times = "default";
        TestParameters::input(param).dft_functional = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_ccp_rmesh_times, "1");

        TestParameters::input(param).exx_ccp_rmesh_times = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_opt_orb_ecut
        auto it = find_label("exx_opt_orb_ecut", readinput.input_lists);
        TestParameters::input(param).exx_opt_orb_ecut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_opt_orb_tolerence
        auto it = find_label("exx_opt_orb_tolerence", readinput.input_lists);
        TestParameters::input(param).exx_opt_orb_tolerence = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // exx_symmetry_realspace
        auto it = find_label("exx_symmetry_realspace", readinput.input_lists);
        TestParameters::input(param).exx_symmetry_realspace = true;
        TestParameters::input(param).symmetry="0";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).exx_symmetry_realspace, false);
    }
    { // rpa_ccp_rmesh_times
        auto it = find_label("rpa_ccp_rmesh_times", readinput.input_lists);
        TestParameters::input(param).rpa_ccp_rmesh_times = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // berry_phase
        auto it = find_label("berry_phase", readinput.input_lists);
        TestParameters::input(param).berry_phase = true;
        TestParameters::input(param).basis_type = "pw";
        TestParameters::input(param).calculation = "nscf";
        TestParameters::input(param).gdir = 1;

        TestParameters::input(param).gdir = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).calculation = "scf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).basis_type = "lcao_in_pw";
        TestParameters::input(param).calculation = "nscf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
    }
    { // towannier90
        auto it = find_label("towannier90", readinput.input_lists);
        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).calculation = "nscf";
        TestParameters::input(param).nspin = 2;
        TestParameters::input(param).wannier_spin = "none";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("towannier90"));

        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).calculation = "scf";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("towannier90"));
    }
    { // wannier_method
        auto it = find_label("wannier_method", readinput.input_lists);
        TestParameters::input(param).towannier90 = true;
        TestParameters::input(param).basis_type = "lcao_in_pw";
        TestParameters::input(param).wannier_method = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).wannier_method, 1);
    }
    { // of_hold_rho0
        auto it = find_label("of_hold_rho0", readinput.input_lists);
        TestParameters::input(param).of_wt_rho0 = 1;
        TestParameters::input(param).of_hold_rho0 = false;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).of_hold_rho0, true);
    }
    { // of_full_pw_dim
        auto it = find_label("of_full_pw_dim", readinput.input_lists);
        TestParameters::input(param).of_full_pw = false;
        TestParameters::input(param).of_full_pw_dim = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).of_full_pw_dim, 0);
    }
    { // of_read_kernel
        auto it = find_label("of_read_kernel", readinput.input_lists);
        TestParameters::input(param).of_read_kernel = true;
        TestParameters::input(param).of_kinetic = "none";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).of_read_kernel, false);
    }
    { // of_ml_gene_data
        auto it = find_label("of_ml_gene_data", readinput.input_lists);

        TestParameters::input(param).of_ml_gene_data = true;
        TestParameters::input(param).esolver_type = "ofdft";
        TestParameters::input(param).basis_type = "pw";
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).of_ml_gene_data = true;
        TestParameters::input(param).esolver_type = "ksdft";
        TestParameters::input(param).basis_type = "lcao";
        GlobalV::NPROC = 1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).of_ml_gene_data = true;
        TestParameters::input(param).esolver_type = "ksdft";
        TestParameters::input(param).basis_type = "pw";
        GlobalV::NPROC = 2;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // dft_plus_u
        auto it = find_label("dft_plus_u", readinput.input_lists);
        TestParameters::input(param).dft_plus_u = 1;
        TestParameters::input(param).orbital_corr = {-1, -1};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).dft_plus_u, 0);
    }
    { // uramping
        auto it = find_label("uramping", readinput.input_lists);
        TestParameters::sys(param).uramping = 1;
        TestParameters::input(param).orbital_corr = {-1, -1};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::sys(param).uramping, 0);
    }
    { // onsite_radius
        auto it = find_label("onsite_radius", readinput.input_lists);
        TestParameters::input(param).onsite_radius = 0.0;
        TestParameters::input(param).dft_plus_u = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).onsite_radius, 3.0);
    }
    { // hubbard_u
        auto it = find_label("hubbard_u", readinput.input_lists);
        TestParameters::input(param).ntype = 2;
        it->second.str_values = {"1.0", "2.0"};
        TestParameters::sys(param).hubbard_u = {1.0, 2.0};
        it->second.check_value(it->second, param);
        TestParameters::input(param).ntype = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ntype = 2;
        TestParameters::sys(param).hubbard_u = {1.0, -1.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // orbital_corr
        auto it = find_label("orbital_corr", readinput.input_lists);
        TestParameters::input(param).ntype = 2;
        it->second.str_values = {"1", "2"};
        TestParameters::input(param).orbital_corr = {1, 2};
        it->second.check_value(it->second, param);
        TestParameters::input(param).ntype = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).ntype = 2;
        TestParameters::input(param).orbital_corr = {1, 4};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_nao_ecut
        auto it = find_label("bessel_nao_ecut", readinput.input_lists);
        TestParameters::input(param).bessel_nao_ecut = "default";
        TestParameters::input(param).ecutwfc = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(std::stod(TestParameters::input(param).bessel_nao_ecut), 1.0);

        TestParameters::input(param).bessel_nao_ecut = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_nao_rcut
        auto it = find_label("bessel_nao_rcut", readinput.input_lists);
        TestParameters::input(param).bessel_nao_rcuts = {-1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_descriptor_ecut
        auto it = find_label("bessel_descriptor_ecut", readinput.input_lists);
        TestParameters::input(param).bessel_descriptor_ecut = "default";
        TestParameters::input(param).ecutwfc = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(std::stod(TestParameters::input(param).bessel_descriptor_ecut), 1.0);

        TestParameters::input(param).bessel_descriptor_ecut = "-1";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // bessel_descriptor_rcut
        auto it = find_label("bessel_descriptor_rcut", readinput.input_lists);
        TestParameters::input(param).bessel_descriptor_rcut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_mag_switch
        auto it = find_label("sc_mag_switch", readinput.input_lists);
        TestParameters::input(param).sc_mag_switch = true;
        // Since sc_mag_switch check is disabled, just call the function without expecting exit
        it->second.check_value(it->second, param);
    }
    { // sc_thr
        auto it = find_label("sc_thr", readinput.input_lists);
        TestParameters::input(param).sc_thr = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nsc
        auto it = find_label("nsc", readinput.input_lists);
        TestParameters::input(param).nsc = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nsc_min
        auto it = find_label("nsc_min", readinput.input_lists);
        TestParameters::input(param).nsc_min = 0;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // alpha_trial
        auto it = find_label("alpha_trial", readinput.input_lists);
        TestParameters::input(param).alpha_trial = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sccut
        auto it = find_label("sccut", readinput.input_lists);
        TestParameters::input(param).sccut = -1;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // sc_scf_thr
        auto it = find_label("sc_scf_thr", readinput.input_lists);
        TestParameters::input(param).sc_scf_thr = -1e-3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // qo_thr
        auto it = find_label("qo_thr", readinput.input_lists);
        TestParameters::input(param).qo_thr = 1e-5;
        it->second.check_value(it->second, param);
    }
    { // qo_strategy
        auto it = find_label("qo_strategy", readinput.input_lists);
        TestParameters::input(param).ntype = 2;

        TestParameters::input(param).qo_basis = "hydrogen";
        TestParameters::input(param).qo_strategy = {"all"};
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).qo_strategy.size(), 2);
        EXPECT_EQ(TestParameters::input(param).qo_strategy[0], "all");
        EXPECT_EQ(TestParameters::input(param).qo_strategy[1], "all");

        TestParameters::input(param).qo_strategy = {};
        TestParameters::input(param).qo_basis = "hydrogen";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).qo_strategy.size(), 2);
        EXPECT_EQ(TestParameters::input(param).qo_strategy[0], "energy-valence");
        EXPECT_EQ(TestParameters::input(param).qo_strategy[1], "energy-valence");

        TestParameters::input(param).qo_strategy = {};
        TestParameters::input(param).qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).qo_strategy.size(), 2);
        EXPECT_EQ(TestParameters::input(param).qo_strategy[0], "all");
        EXPECT_EQ(TestParameters::input(param).qo_strategy[1], "all");

        TestParameters::input(param).qo_basis = "test";
        TestParameters::input(param).qo_strategy = {};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // qo_screening_coeff
        auto it = find_label("qo_screening_coeff", readinput.input_lists);
        TestParameters::input(param).ntype = 2;

        it->second.str_values = {"0.2"};
        TestParameters::input(param).qo_screening_coeff = {0.2};
        TestParameters::input(param).qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff.size(), 2);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff[0], 0.2);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff[1], 0.2);

        TestParameters::input(param).qo_screening_coeff = {};
        TestParameters::input(param).qo_basis = "pswfc";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff.size(), 2);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff[0], 0.1);
        EXPECT_EQ(TestParameters::input(param).qo_screening_coeff[1], 0.1);

        TestParameters::input(param).qo_screening_coeff = {};
        TestParameters::input(param).qo_basis = "test";
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.reset_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).qo_screening_coeff = {0.2, -0.1};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).qo_screening_coeff = {0.2, 1e-8};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // diag_subspace
        auto it = find_label("diag_subspace", readinput.input_lists);
        EXPECT_EQ(it->second.get_availability(), "basis_type==pw and ks_solver==dav_subspace");
    }
    { // md_nstep
        auto it = find_label("md_nstep", readinput.input_lists);
        TestParameters::input(param).mdp.md_nstep = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mdp.md_nstep, 50);
    }
    { // md_prec_level
        auto it = find_label("md_prec_level", readinput.input_lists);
        TestParameters::input(param).mdp.md_prec_level = 1;
        TestParameters::input(param).calculation = "vc-relax";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mdp.md_prec_level, 0);

        TestParameters::input(param).calculation = "md";
        TestParameters::input(param).mdp.md_prec_level = 1;
        TestParameters::input(param).mdp.md_type = "vc-md";
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mdp.md_prec_level, 0);
    }
    { // md_tfreq
        auto it = find_label("md_tfreq", readinput.input_lists);
        TestParameters::input(param).mdp.md_tfreq = 0;
        TestParameters::input(param).calculation = "md";
        TestParameters::input(param).mdp.md_dt = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mdp.md_tfreq, 1.0 / 40 / 1.0);
    }
    { // md_pfreq
        auto it = find_label("md_pfreq", readinput.input_lists);
        TestParameters::input(param).mdp.md_pfreq = 0;
        TestParameters::input(param).calculation = "md";
        TestParameters::input(param).mdp.md_dt = 1.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).mdp.md_pfreq, 1.0 / 400 / 1.0);
    }
    { // lj_rule
        auto it = find_label("lj_rule", readinput.input_lists);
        TestParameters::input(param).esolver_type = "lj";
        TestParameters::input(param).mdp.lj_rule = 3;
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_rcut
        auto it = find_label("lj_rcut", readinput.input_lists);
        TestParameters::input(param).ntype = 2;
        TestParameters::input(param).esolver_type = "lj";
        it->second.str_values = {"1.0", "2.0"};
        TestParameters::input(param).mdp.lj_rcut = {1.0, 2.0};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

        TestParameters::input(param).mdp.lj_rcut = {1.0, 2.0, -1.0};
        it->second.str_values = {"1.0", "2.0", "-1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_epsilon
        auto it = find_label("lj_epsilon", readinput.input_lists);
        TestParameters::input(param).ntype = 2;
        TestParameters::input(param).esolver_type = "lj";
        TestParameters::input(param).mdp.lj_epsilon = {1.0};
        it->second.str_values = {"1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // lj_sigma
        auto it = find_label("lj_sigma", readinput.input_lists);
        TestParameters::input(param).ntype = 2;
        TestParameters::input(param).esolver_type = "lj";
        TestParameters::input(param).mdp.lj_sigma = {1.0};
        it->second.str_values = {"1.0"};
        testing::internal::CaptureStdout();
        EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        output = testing::internal::GetCapturedStdout();
        EXPECT_THAT(output, testing::HasSubstr("NOTICE"));
    }
    { // nocc 
        auto it = find_label("nocc", readinput.input_lists);
        TestParameters::input(param).nocc = 5;
        TestParameters::input(param).nbands = 4;
        TestParameters::input(param).nelec = 0.0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).nocc, 4);
        TestParameters::input(param).nocc = 0;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).nocc, 4);
    }
}

TEST_F(InputTest, Item_test_out_mat_vec)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;

    { // out_mat_dh
        auto it = find_label("out_mat_dh", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_dh[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_dh[1], 8);

        it->second.str_values = {"1", "12"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_dh[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_dh[1], 12);
    }
    { // out_mat_ds
        auto it = find_label("out_mat_ds", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_ds[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_ds[1], 8);

        it->second.str_values = {"1", "10"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_ds[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_ds[1], 10);
    }
    { // out_mat_t
        auto it = find_label("out_mat_t", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_t[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_t[1], 8);

        it->second.str_values = {"1", "15"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_t[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_t[1], 15);
    }
    { // out_mat_r
        auto it = find_label("out_mat_r", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_r[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_r[1], 8);

        it->second.str_values = {"1", "6"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_r[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_r[1], 6);
    }
    { // out_mat_xc2
        auto it = find_label("out_mat_xc2", readinput.input_lists);
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_xc2[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_xc2[1], 8);

        it->second.str_values = {"1", "9"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_mat_xc2[0], 1);
        EXPECT_EQ(TestParameters::input(param).out_mat_xc2[1], 9);
    }
}

TEST_F(InputTest, OutStru)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    auto it = find_label("out_stru", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());

    // --- Valid numeric values ---
    {
        it->second.str_values = {"0"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }
    {
        it->second.str_values = {"1"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }
    {
        it->second.str_values = {"2"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 2);
    }

    // --- Backward-compatible boolean aliases (true -> 1, false -> 0) ---
    {
        it->second.str_values = {"true"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }
    {
        it->second.str_values = {"TRUE"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }
    {
        it->second.str_values = {".true."};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }
    {
        it->second.str_values = {"Yes"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }
    {
        it->second.str_values = {"false"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }
    {
        it->second.str_values = {"FALSE"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }
    {
        it->second.str_values = {".false."};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }
    {
        it->second.str_values = {"No"};
        it->second.read_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);
    }

    // --- Valid value check_value passes ---
    {
        for (const int v : {0, 1, 2})
        {
            TestParameters::input(param).out_stru = v;
            // Expect no exit / no crash; check_value is a void function that only
            // calls WARNING_QUIT on bad input.
            it->second.check_value(it->second, param);
        }
    }

    // --- reset_value: calculation in offlist forces out_stru to 0 ---
    {
        TestParameters::input(param).calculation = "get_wf";
        TestParameters::input(param).out_stru = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);

        TestParameters::input(param).calculation = "nscf";
        TestParameters::input(param).out_stru = 2;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 0);

        // Non-offlist calculation preserves value
        TestParameters::input(param).calculation = "cell-relax";
        TestParameters::input(param).out_stru = 1;
        it->second.reset_value(it->second, param);
        EXPECT_EQ(TestParameters::input(param).out_stru, 1);
    }

    // --- Invalid integer values -> WARNING_QUIT via check_value ---
    {
        for (const std::string& s : {"3", "-1", "-2", "4", "10"})
        {
            it->second.str_values = {s};
            it->second.read_value(it->second, param);
            EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
        }
    }

    // --- Non-numeric / malformed inputs -> WARNING_QUIT via read_value ---
    {
        for (const std::string& s : {"abc", "2.5", "2abc", "-1abc", "xyz", ""})
        {
            it->second.str_values = {s};
            EXPECT_EXIT(it->second.read_value(it->second, param), ::testing::ExitedWithCode(1), "");
        }
    }
}

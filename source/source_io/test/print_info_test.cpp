#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/test_parameters.h"
#define private public
#include "source_cell/klist.h"
#include "source_cell/parallel_kpoints.h"
#include "source_cell/unitcell.h"
#include "source_io/module_unk/berryphase.h"
#include "source_io/module_output/print_info.h"
#include "prepare_unitcell.h"
#undef private
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

bool berryphase::berry_phase_flag=false;


/************************************************
 *  unit test of print_info.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - print_parameters()
 *     - setup calculation parameters
 */

class PrintInfoTest : public testing::Test
{
protected:
	std::string output;
	UnitCell* ucell;
	K_Vectors* kv;
	void SetUp()
	{
		ucell = new UnitCell;
		kv = new K_Vectors;
	}
	void TearDown()
	{
		delete ucell;
		delete kv;
	}
};

TEST_F(PrintInfoTest, SetupParameters)
{
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo();
	std::string k_file = "./support/KPT";
	kv->spin_mult = 1;
	const bool gamma_only_local = false;
	const double kspacing[3] = {0.0, 0.0, 0.0};
	const std::string kmesh_type = "gamma";
	const double koffset[3] = {0.0, 0.0, 0.0};
	kv->read_kpoints(*ucell, k_file, gamma_only_local, kspacing, kmesh_type, koffset, GlobalV::ofs_running, GlobalV::ofs_warning, GlobalV::MY_RANK);
	EXPECT_EQ(kv->get_nkstot(),512);
	std::vector<std::string> cal_type = {"scf","relax","cell-relax","md"};
	std::vector<std::string> md_types = {"fire","nve","nvt","npt","langevin","msst"};
	GlobalV::MY_RANK = 0;
	for(int i=0; i<cal_type.size(); ++i)
	{
		if(cal_type[i] != "md")
		{
			TestParameters::sys().gamma_only_local = false;
			TestParameters::input().calculation = cal_type[i];
			testing::internal::CaptureStdout();
            EXPECT_NO_THROW(ModuleIO::print_parameters(*ucell, *kv, TestParameters::input()));
            output = testing::internal::GetCapturedStdout();
			if(TestParameters::input().calculation == "scf")
			{
				EXPECT_THAT(output,testing::HasSubstr("Self-consistent calculations"));
			}
			else if(TestParameters::input().calculation == "relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Ion relaxation calculations"));
			}
			else if(TestParameters::input().calculation == "cell-relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Cell relaxation calculations"));
			}
		}
		else
		{
			TestParameters::sys().gamma_only_local = true;
            TestParameters::input().calculation = cal_type[i];
			for(int j=0; j<md_types.size(); ++j)
			{
                TestParameters::input().mdp.md_type = md_types[j];
                testing::internal::CaptureStdout();
                EXPECT_NO_THROW(ModuleIO::print_parameters(*ucell, *kv, TestParameters::input()));
                output = testing::internal::GetCapturedStdout();
                EXPECT_THAT(output,testing::HasSubstr("Molecular Dynamics simulations"));
                if (PARAM.mdp.md_type == "fire")
                {
                    EXPECT_THAT(output,testing::HasSubstr("FIRE"));
                }
                else if (PARAM.mdp.md_type == "nve")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NVE"));
                }
                else if (PARAM.mdp.md_type == "nvt")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NVT"));
                }
                else if (PARAM.mdp.md_type == "npt")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NPT"));
                }
                else if (PARAM.mdp.md_type == "langevin")
                {
                    EXPECT_THAT(output,testing::HasSubstr("Langevin"));
                }
                else if (PARAM.mdp.md_type == "msst")
                {
                    EXPECT_THAT(output,testing::HasSubstr("MSST"));
                }
			}
		}
	}
	std::vector<std::string> basis_type = {"lcao","pw","lcao_in_pw"};
	for(int i=0; i<basis_type.size(); ++i)
	{
		TestParameters::input().basis_type = basis_type[i];
		testing::internal::CaptureStdout();
        EXPECT_NO_THROW(ModuleIO::print_parameters(*ucell, *kv, TestParameters::input()));
        output = testing::internal::GetCapturedStdout();
		if(TestParameters::input().basis_type == "lcao")
		{
			EXPECT_THAT(output,testing::HasSubstr("Use Systematically Improvable Atomic bases"));
		}
		else if(TestParameters::input().basis_type == "lcao_in_pw")
		{
			EXPECT_THAT(output,testing::HasSubstr("Expand Atomic bases into plane waves"));
		}
		else if(TestParameters::input().basis_type == "pw")
		{
			EXPECT_THAT(output,testing::HasSubstr("Use plane wave basis"));
		}
	}
}

TEST_F(PrintInfoTest, PrintScreen)
{
	int stress_step = 11;
	int force_step = 101;
	int istep = 1001;
	std::vector<std::string> cal_type = {"scf","nscf","md","relax","cell-relax"};
	for(int i=0; i<cal_type.size(); ++i)
	{
		TestParameters::input().calculation = cal_type[i];
		if(TestParameters::input().calculation=="scf")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("SELF-CONSISTENT"));
		}
		else if(TestParameters::input().calculation=="nscf")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("NONSELF-CONSISTENT"));
		}
		else if(TestParameters::input().calculation=="md")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("STEP OF MOLECULAR DYNAMICS"));
		}
		else
		{
			if(TestParameters::input().calculation=="relax")
			{
				testing::internal::CaptureStdout();
                ModuleIO::print_screen(stress_step, force_step, istep);
                output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("RELAX STEP"));
			}
			else if(TestParameters::input().calculation=="cell-relax")
			{
				testing::internal::CaptureStdout();
                ModuleIO::print_screen(stress_step, force_step, istep);
                output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("RELAX STEP"));
			EXPECT_THAT(output,testing::HasSubstr("CELL#"));
			EXPECT_THAT(output,testing::HasSubstr("IONS#"));
			}
		}
	}
}

TEST_F(PrintInfoTest, PrintTime)
{
	time_t time_start = std::time(nullptr);
	time_t time_finish = std::time(nullptr);
	testing::internal::CaptureStdout();
    EXPECT_NO_THROW(ModuleIO::print_time(time_start, time_finish));
    output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("START  Time"));
	EXPECT_THAT(output,testing::HasSubstr("FINISH Time"));
	EXPECT_THAT(output,testing::HasSubstr("TOTAL  Time"));
}

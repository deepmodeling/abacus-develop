#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#include "module_cell/klist.h"
#include "module_cell/parallel_kpoints.h"
#include "module_cell/unitcell.h"
#include "module_io/berryphase.h"
#include "module_io/print_info.h"
#include "prepare_unitcell.h"
#undef private
#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}
#endif
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

bool berryphase::berry_phase_flag=false;

namespace GlobalC
{
	Parallel_Kpoints Pkpoints;
	UnitCell ucell;
}

/************************************************
 *  unit test of print_info.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - setup_parameters()
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
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->get_nkstot(),512);
	std::vector<std::string> cal_type = {"scf","relax","cell-relax","md"};
	std::vector<std::string> md_types = {"fire","nve","nvt","npt","langevin","msst"};
	GlobalV::MY_RANK = 0;
	for(int i=0; i<cal_type.size(); ++i)
	{
		if(cal_type[i] != "md")
		{
			PARAM.sys.gamma_only_local = false;
			PARAM.input.calculation = cal_type[i];
			testing::internal::CaptureStdout();
            EXPECT_NO_THROW(ModuleIO::setup_parameters(*ucell, *kv));
            output = testing::internal::GetCapturedStdout();
			if(PARAM.input.calculation == "scf")
			{
				EXPECT_THAT(output,testing::HasSubstr("Self-consistent calculations"));
			}
			else if(PARAM.input.calculation == "relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Ion relaxation calculations"));
			}
			else if(PARAM.input.calculation == "cell-relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Cell relaxation calculations"));
			}
		}
		else
		{
			PARAM.sys.gamma_only_local = true;
      PARAM.input.calculation = cal_type[i];
			for(int j=0; j<md_types.size(); ++j)
			{
                PARAM.input.mdp.md_type = md_types[j];
                testing::internal::CaptureStdout();
                EXPECT_NO_THROW(ModuleIO::setup_parameters(*ucell, *kv));
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
		PARAM.input.basis_type = basis_type[i];
		testing::internal::CaptureStdout();
        EXPECT_NO_THROW(ModuleIO::setup_parameters(*ucell, *kv));
        output = testing::internal::GetCapturedStdout();
		if(PARAM.input.basis_type == "lcao")
		{
			EXPECT_THAT(output,testing::HasSubstr("Use Systematically Improvable Atomic bases"));
		}
		else if(PARAM.input.basis_type == "lcao_in_pw")
		{
			EXPECT_THAT(output,testing::HasSubstr("Expand Atomic bases into plane waves"));
		}
		else if(PARAM.input.basis_type == "pw")
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
		PARAM.input.calculation = cal_type[i];
		if(PARAM.input.calculation=="scf")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("SELF-CONSISTENT"));
		}
		else if(PARAM.input.calculation=="nscf")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("NONSELF-CONSISTENT"));
		}
		else if(PARAM.input.calculation=="md")
		{
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("STEP OF MOLECULAR DYNAMICS"));
		}
		else
		{
			PARAM.input.relax_new = false;
			if(PARAM.input.calculation=="relax")
			{
				testing::internal::CaptureStdout();
                ModuleIO::print_screen(stress_step, force_step, istep);
                output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("STEP OF ION RELAXATION"));
			}
			else if(PARAM.input.calculation=="cell-relax")
			{
				testing::internal::CaptureStdout();
                ModuleIO::print_screen(stress_step, force_step, istep);
                output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("RELAX CELL"));
				EXPECT_THAT(output,testing::HasSubstr("RELAX IONS"));
			}
			PARAM.input.relax_new = true;
			testing::internal::CaptureStdout();
            ModuleIO::print_screen(stress_step, force_step, istep);
            output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("STEP OF RELAXATION"));
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

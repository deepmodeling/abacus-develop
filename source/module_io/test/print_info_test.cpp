#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_io/print_info.h"
#include "module_io/input.h"
#include "prepare_unitcell.h"
#include "module_cell/unitcell.h"
#include "module_cell/klist.h"
#include "module_cell/parallel_kpoints.h"
#include "module_io/berryphase.h"

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

bool berryphase::berry_phase_flag=0;

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

TEST_F(PrintInfoTest, Constructor)
{
	EXPECT_NO_THROW(Print_Info pinfo);
}

TEST_F(PrintInfoTest, SetupParameters)
{
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo();
	std::string k_file = "./support/KPT";
	kv->nspin = 1;
	kv->read_kpoints(k_file);
	EXPECT_EQ(kv->nkstot,512);
	std::vector<std::string> cal_type = {"scf","relax","cell-relax","md"};
	std::vector<std::string> md_types = {"fire","nve","nvt","npt","langevin","msst"};
	GlobalV::MY_RANK = 0;
	for(int i=0; i<cal_type.size(); ++i)
	{
		if(cal_type[i] != "md")
		{
			GlobalV::COLOUR = false;
			GlobalV::GAMMA_ONLY_LOCAL = false;
			GlobalV::CALCULATION = cal_type[i];
			testing::internal::CaptureStdout();
			EXPECT_NO_THROW(Print_Info::setup_parameters(*ucell,*kv));
			output = testing::internal::GetCapturedStdout();
			if(GlobalV::CALCULATION == "scf")
			{
				EXPECT_THAT(output,testing::HasSubstr("Self-consistent calculations"));
			}
			else if(GlobalV::CALCULATION == "relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Ion relaxation calculations"));
			}
			else if(GlobalV::CALCULATION == "cell-relax")
			{
				EXPECT_THAT(output,testing::HasSubstr("Cell relaxation calculations"));
			}
		}
		else
		{
			GlobalV::COLOUR = true;
			GlobalV::GAMMA_ONLY_LOCAL = true;
			GlobalV::CALCULATION = cal_type[i];
			for(int j=0; j<md_types.size(); ++j)
			{
				INPUT.mdp.md_type = md_types[j];
                testing::internal::CaptureStdout();
                EXPECT_NO_THROW(Print_Info::setup_parameters(*ucell,*kv));
                output = testing::internal::GetCapturedStdout();
                EXPECT_THAT(output,testing::HasSubstr("Molecular Dynamics simulations"));
                if(INPUT.mdp.md_type == "fire")
                {
                    EXPECT_THAT(output,testing::HasSubstr("FIRE"));
                }
                else if(INPUT.mdp.md_type == "nve")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NVE"));
                }
                else if(INPUT.mdp.md_type == "nvt")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NVT"));
                }
                else if(INPUT.mdp.md_type == "npt")
                {
                    EXPECT_THAT(output,testing::HasSubstr("NPT"));
                }
                else if(INPUT.mdp.md_type == "langevin")
                {
                    EXPECT_THAT(output,testing::HasSubstr("Langevin"));
                }
                else if(INPUT.mdp.md_type == "msst")
                {
                    EXPECT_THAT(output,testing::HasSubstr("MSST"));
                }
			}
		}
	}
	std::vector<std::string> basis_type = {"lcao","pw","lcao_in_pw"};
	for(int i=0; i<basis_type.size(); ++i)
	{
		GlobalV::BASIS_TYPE = basis_type[i];
		testing::internal::CaptureStdout();
		EXPECT_NO_THROW(Print_Info::setup_parameters(*ucell,*kv));
		output = testing::internal::GetCapturedStdout();
		if(GlobalV::BASIS_TYPE == "lcao")
		{
			EXPECT_THAT(output,testing::HasSubstr("Use Systematically Improvable Atomic bases"));
		}
		else if(GlobalV::BASIS_TYPE == "lcao_in_pw")
		{
			EXPECT_THAT(output,testing::HasSubstr("Expand Atomic bases into plane waves"));
		}
		else if(GlobalV::BASIS_TYPE == "pw")
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
		GlobalV::CALCULATION = cal_type[i];
		if(GlobalV::CALCULATION=="scf")
		{
			testing::internal::CaptureStdout();
			Print_Info::print_screen(stress_step,force_step,istep);
			output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("SELF-CONSISTENT"));
		}
		else if(GlobalV::CALCULATION=="nscf")
		{
			testing::internal::CaptureStdout();
			Print_Info::print_screen(stress_step,force_step,istep);
			output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("NONSELF-CONSISTENT"));
		}
		else if(GlobalV::CALCULATION=="md")
		{
			testing::internal::CaptureStdout();
			Print_Info::print_screen(stress_step,force_step,istep);
			output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("STEP OF MOLECULAR DYNAMICS"));
		}
		else
		{
			GlobalV::relax_new = false;
			if(GlobalV::CALCULATION=="relax")
			{
				testing::internal::CaptureStdout();
				Print_Info::print_screen(stress_step,force_step,istep);
				output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("STEP OF ION RELAXATION"));
			}
			else if(GlobalV::CALCULATION=="cell-relax")
			{
				testing::internal::CaptureStdout();
				Print_Info::print_screen(stress_step,force_step,istep);
				output = testing::internal::GetCapturedStdout();
				EXPECT_THAT(output,testing::HasSubstr("RELAX CELL"));
				EXPECT_THAT(output,testing::HasSubstr("RELAX IONS"));
			}
			GlobalV::relax_new = true;
			testing::internal::CaptureStdout();
			Print_Info::print_screen(stress_step,force_step,istep);
			output = testing::internal::GetCapturedStdout();
			EXPECT_THAT(output,testing::HasSubstr("STEP OF RELAXATION"));
		}
	}
}

TEST_F(PrintInfoTest, PrintTime)
{
	time_t time_start = std::time(NULL);
	time_t time_finish = std::time(NULL);
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(Print_Info::print_time(time_start,time_finish));
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("START  Time"));
	EXPECT_THAT(output,testing::HasSubstr("FINISH Time"));
	EXPECT_THAT(output,testing::HasSubstr("TOTAL  Time"));
}
#undef private

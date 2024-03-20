#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include <vector>
#include <valarray>
#ifdef __MPI
#include "mpi.h"
#endif
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
#endif
Magnetism::Magnetism()
{
	this->tot_magnetization = 0.0;
	this->abs_magnetization = 0.0;
	this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
	delete[] this->start_magnetization;
}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadCellPPWarning1
 *     - read_cell_pseudopots(): error when average the pseudopotential: error_ap
 *   - ReadCellPPWarning2
 *     - read_cell_pseudopots(): Couldn't find pseudopotential file:: error == 1
 *   - ReadCellPPWarning3
 *     - read_cell_pseudopots(): Pseudopotential data do not match: error ==2
 *     - error==3 is currently difficult to reach in read_pseudo_vwr
 *   - ReadCellPPWarning4
 *     - read_cell_pseudopots(): dft_functional from INPUT does not match that in pseudopot file
 *   - ReadCellPPWarning5
 *     - read_cell_pseudopots(): Unknown pseudopotential type
 *   - ReadCellPP
 *     - read_cell_pseudopots(): read pp files with flag_empty_element set
 *   - CalMeshx
 *     - cal_meshx(): calculate max mesh info from atomic pseudo potential file
 *   - CalNatomwfc1
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo potential file
 *     - NSPIN != 4
 *     - this corresponds to number_of_wfc, PP_CHI in pp file, and atoms[it].ncpp.lchi[ncpp.nchi]
 *     - setup the total number of PAOs: pseudopotential atomic orbitals
 *   - CalNatomwfc2
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo potential file
 *     - NSPIN ==4, has_so = false
 *   - CalNatomwfc3
 *     - cal_natomwfc(): calculate total number of atomic orbitals in pseudo potential file
 *     - NSPIN ==4, has_so = true
 *   - CalNwfc1
 *     - cal_nwfc(): calcuate the total number of local basis: NSPIN != 4
 *     - this corresponds to number_of_proj, PP_BETA in pp file, and atoms[it].l_nchi[nw], nw from orb file
 *     - setup GlobalV::NLOCAL
 *     - interfaces initialed in this function:
 *       - itia2iat
 *       - iat2iwt
 *       - itiaiw2iwt
 *       - iwt2iat
 *       - iwt2iw
 *   - CalNwfc2
 *     - cal_nwfc(): calcuate the total number of local basis: NSPIN == 4
 *   - CheckStructure
 *     - check_structure(): check if too atoms are two close
 *   - ReadPseudoWarning1
 *     - read_pseudo(): All DFT functional must consistent.
 *   - ReadPseudoWarning2
 *     - read_pseudo(): number valence electrons > corresponding minimum possible of an element
 *   - CalNelec: UnitCell::cal_nelec
 *     - calculate the total number of valence electrons from psp files
 */

//mock function
#ifdef __LCAO
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}
#endif

class UcellTest : public ::testing::Test
{
protected:
	UcellTestPrepare utp = UcellTestLib["C1H2-Read"];
	std::unique_ptr<UnitCell> ucell;
	std::ofstream ofs;
	std::string pp_dir;
	std::string output;
	void SetUp()
	{
		ofs.open("running.log");
		GlobalV::relax_new = utp.relax_new;
		GlobalV::global_out_dir = "./";
		ucell = utp.SetUcellInfo();
		GlobalV::LSPINORB = false;
		pp_dir = "./support/";
		GlobalV::PSEUDORCUT = 15.0;
		GlobalV::DFT_FUNCTIONAL = "default";
		GlobalV::test_unitcell = 1;
		GlobalV::test_pseudo_cell = 1;
		GlobalV::NSPIN = 1;
		GlobalV::BASIS_TYPE = "pw";
	}
	void TearDown()
	{
		ofs.close();
	}
};

using UcellDeathTest = UcellTest;

TEST_F(UcellDeathTest,ReadCellPPWarning1)
{
	GlobalV::LSPINORB = true;
	ucell->pseudo_fn[1] = "H_sr.upf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_cell_pseudopots(pp_dir,ofs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("error when average the pseudopotential."));
}

TEST_F(UcellDeathTest,ReadCellPPWarning2)
{
	pp_dir = "./arbitrary/";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_cell_pseudopots(pp_dir,ofs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Couldn't find pseudopotential file"));
}

TEST_F(UcellDeathTest,ReadCellPPWarning3)
{
	ucell->pseudo_type[0] = "upf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_cell_pseudopots(pp_dir,ofs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Pseudopotential data do not match."));
}

TEST_F(UcellDeathTest,ReadCellPPWarning4)
{
	GlobalV::DFT_FUNCTIONAL = "LDA";
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(ucell->read_cell_pseudopots(pp_dir,ofs));
	output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("dft_functional readin is: LDA"));
    EXPECT_THAT(output, testing::HasSubstr("dft_functional in pseudopot file is: PBE"));
    EXPECT_THAT(output, testing::HasSubstr("Please make sure this is what you need"));
}

TEST_F(UcellDeathTest, ReadCellPPWarning5)
{
    ucell->pseudo_type[0] = "upf0000";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ucell->read_cell_pseudopots(pp_dir, ofs), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("Unknown pseudopotential type."));
}

TEST_F(UcellTest,ReadCellPP)
{
	ucell->atoms[1].flag_empty_element = true;
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_EQ(ucell->atoms[0].ncpp.pp_type,"NC");
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so); //becomes false in average_p
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	EXPECT_EQ(ucell->atoms[0].ncpp.nchi,2); //3=>2 in average_p
	EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	ofs.close();
	std::ifstream ifs;
    	ifs.open("running.log");
    	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    	EXPECT_THAT(str, testing::HasSubstr("Read in pseudopotential file is C.upf"));
    	EXPECT_THAT(str, testing::HasSubstr("pseudopotential type = NC"));
    	EXPECT_THAT(str, testing::HasSubstr("exchange-correlation functional = PBE"));
    	EXPECT_THAT(str, testing::HasSubstr("valence electrons = 4"));
    	EXPECT_THAT(str, testing::HasSubstr("Read in pseudopotential file is H.upf"));
    	EXPECT_THAT(str, testing::HasSubstr("valence electrons = 0"));
}

TEST_F(UcellTest,CalMeshx)
{
	ucell->read_cell_pseudopots(pp_dir,ofs);
	ucell->cal_meshx();
	EXPECT_EQ(ucell->atoms[0].ncpp.msh,1247);
	EXPECT_EQ(ucell->atoms[1].ncpp.msh,1165);
	EXPECT_EQ(ucell->meshx,1247);
}

TEST_F(UcellTest,CalNatomwfc1)
{
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	ucell->cal_natomwfc(ofs);
	EXPECT_EQ(ucell->atoms[0].ncpp.nchi,2);
	EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	EXPECT_EQ(ucell->atoms[0].na,1);
	EXPECT_EQ(ucell->atoms[1].na,2);
	EXPECT_EQ(ucell->natomwfc,(1+3)*1+1*2);
}

TEST_F(UcellTest,CalNatomwfc2)
{
	GlobalV::LSPINORB = false;
	GlobalV::NSPIN = 4;
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	ucell->cal_natomwfc(ofs);
	EXPECT_EQ(ucell->atoms[0].ncpp.nchi,2);
	EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	EXPECT_EQ(ucell->atoms[0].na,1);
	EXPECT_EQ(ucell->atoms[1].na,2);
	EXPECT_EQ(ucell->natomwfc,((1+3)*1+1*2)*2);
}

TEST_F(UcellTest,CalNatomwfc3)
{
	GlobalV::LSPINORB = true;
	GlobalV::NSPIN = 4;
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_TRUE(ucell->atoms[0].ncpp.has_so);
	EXPECT_TRUE(ucell->atoms[1].ncpp.has_so);
	ucell->cal_natomwfc(ofs);
	EXPECT_EQ(ucell->atoms[0].ncpp.nchi,3);
	EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	EXPECT_EQ(ucell->atoms[0].na,1);
	EXPECT_EQ(ucell->atoms[1].na,2);
	EXPECT_EQ(ucell->natomwfc,((2*0+2)+(2*1+2)+(2*1))*1+(2*0+2)*2);
}

TEST_F(UcellTest,CalNwfc1)
{
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	ucell->cal_nwfc(ofs);
	EXPECT_EQ(ucell->atoms[0].iw2l[8],2);
	EXPECT_EQ(ucell->atoms[0].iw2n[8],0);
	EXPECT_EQ(ucell->atoms[0].iw2m[8],4);
	EXPECT_EQ(ucell->atoms[1].iw2l[8],2);
	EXPECT_EQ(ucell->atoms[1].iw2n[8],0);
	EXPECT_EQ(ucell->atoms[1].iw2m[8],4);
	EXPECT_EQ(ucell->atoms[1].iw2_ylm[8],8);
	//here is the default table for pw basis calculation
	// nw = 1*1 + 3*1 + 5*1 = 9
	//    L N m  L*L+m
	// 0  0 0 0    0
	// 1  1 0 0    0
	// 2  1 0 1    2
	// 3  1 0 2    2
	// 4  2 0 0    4
	// 5  2 0 1    5
	// 6  2 0 2    6
	// 7  2 0 3    7
	// 8  2 0 4    8
	EXPECT_EQ(ucell->atoms[0].na,1);
	EXPECT_EQ(ucell->atoms[1].na,2);
	EXPECT_EQ(ucell->namax,2);
	EXPECT_EQ(ucell->atoms[0].nw,9);
	EXPECT_EQ(ucell->atoms[1].nw,9);
	EXPECT_EQ(ucell->nwmax,9);
	EXPECT_EQ(GlobalV::NLOCAL,3*9);
	//check itia2iat
	EXPECT_EQ(ucell->itia2iat.getSize(), 4);
	EXPECT_EQ(ucell->itia2iat(0,0), 0);
	EXPECT_EQ(ucell->itia2iat(0,1), 0);
	EXPECT_EQ(ucell->itia2iat(1,0), 1);
	EXPECT_EQ(ucell->itia2iat(1,1), 2);
	//check iat2iwt
	EXPECT_EQ(ucell->get_npol(), 1);
	EXPECT_EQ(ucell->get_iat2iwt()[0], 0);
	EXPECT_EQ(ucell->get_iat2iwt()[1], 9);
	EXPECT_EQ(ucell->get_iat2iwt()[2], 18);
	//check itiaiw2iwt
	EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 0), 0);
	EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 1), 1);
	EXPECT_EQ(ucell->itiaiw2iwt(0, 0, 8), 8);
	EXPECT_EQ(ucell->itiaiw2iwt(1, 0, 0), 9);
	EXPECT_EQ(ucell->itiaiw2iwt(1, 1, 0), 18);
	//check itia2iat
	EXPECT_EQ(ucell->itia2iat.getSize(), 4);
	EXPECT_EQ(ucell->itia2iat(0, 0), 0);
	EXPECT_EQ(ucell->itia2iat(0, 1), 0);
	EXPECT_EQ(ucell->itia2iat(1, 0), 1);
	EXPECT_EQ(ucell->itia2iat(1, 1), 2);
	//check iwt2iat
	EXPECT_EQ(ucell->iwt2iat[0], 0);
	EXPECT_EQ(ucell->iwt2iat[10], 1);
	EXPECT_EQ(ucell->iwt2iat[20], 2);
	//check iwt2iw
	EXPECT_EQ(ucell->iwt2iw[0], 0);
	EXPECT_EQ(ucell->iwt2iw[10], 1);
	EXPECT_EQ(ucell->iwt2iw[20], 2);
}

TEST_F(UcellTest,CalNwfc2)
{
	GlobalV::NSPIN = 4;
	GlobalV::BASIS_TYPE = "lcao";
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	ucell->cal_nwfc(ofs);
	EXPECT_EQ(GlobalV::NLOCAL,3*9*2);
}

TEST_F(UcellDeathTest,CheckStructure)
{
	ucell->read_cell_pseudopots(pp_dir,ofs);
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	//trial 1
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(ucell->check_structure(0.2));
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("WARNING: Some atoms are too close!!!"));
	//trial 2
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->check_structure(0.4),::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("The structure is unreasonable!"));
	//trial 3
	ucell->atoms[0].ncpp.psd = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(ucell->check_structure(0.2));
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Notice: symbol 'arbitrary' is not an element symbol!!!! set the covalent radius to be 0."));
}

TEST_F(UcellDeathTest,ReadPseudoWarning1)
{
	GlobalV::global_pseudo_dir = pp_dir;
	GlobalV::out_element_info = 1;
	GlobalV::MIN_DIST_COEF = 0.2;
	ucell->pseudo_fn[1] = "H_sr_lda.upf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(ucell->read_pseudo(ofs),::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("All DFT functional must consistent."));
}

TEST_F(UcellDeathTest,ReadPseudoWarning2)
{
	GlobalV::global_pseudo_dir = pp_dir;
	GlobalV::out_element_info = 1;
	GlobalV::MIN_DIST_COEF = 0.2;
	ucell->pseudo_fn[0] = "Al_ONCV_PBE-1.0.upf";
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(ucell->read_pseudo(ofs));
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Warning: the number of valence electrons in pseudopotential > 3 for Al: [Ne] 3s2 3p1"));
}

TEST_F(UcellTest, CalNelec)
{
    ucell->read_cell_pseudopots(pp_dir, ofs);
    EXPECT_EQ(4,ucell->atoms[0].ncpp.zv);
    EXPECT_EQ(1,ucell->atoms[1].ncpp.zv);
    EXPECT_EQ(1,ucell->atoms[0].na);
    EXPECT_EQ(2,ucell->atoms[1].na);
    double nelec = 0;
    ucell->cal_nelec(nelec);
    EXPECT_DOUBLE_EQ(6,nelec);
}

#ifdef __MPI
#include "mpi.h"
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	testing::InitGoogleTest(&argc, argv);

	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

	int result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}
#endif

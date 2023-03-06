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
 *     - upf.functional_error == 1
 *   - ReadCellPP
 *     - read_cell_pseudopots(): read pp files with flag_empty_element set
 *   - UpdatePosTaud
 *     - update_pos_taud(const double* pos)
 *     - bcast_atoms_tau() is also called in the above function, which calls Atom::bcast_atom with many
 *       atomic info in addition to tau
 *   - BcastUnitcell
 *     - bcast basic info of unitcell and basic info of atoms
 *   - BcastUnitcell2
 *     - calls bcast_atoms2() to bcast atomic pseudo info
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
 *   - CalNwfc2
 *     - cal_nwfc(): calcuate the total number of local basis: NSPIN == 4
 *   - CalUx
 *     - cal_ux(): 
 *   - CheckStructure
 *     - check_structure(): check if too atoms are two close
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
	UcellTestPrepare UTP = UcellTestPrepare("C1H2-Index",	//system-name
				"bcc",		//latname
				2,		//lmaxmax
				true,		//init_vel
				true,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				1.8897261254578281, //lat0
				{10.0,0.0,0.0,	//latvec
				 0.0,10.0,0.0,
				 0.0,0.0,10.0},
				{"C","H"},	//elements
				{"C.upf","H.upf"},	//upf file
				{"upf201","upf201"},	//upf types
				{"C.orb","H.orb"},	//orb file
				{1,2},		//number of each elements
				{12.0,1.0},	//atomic mass
				"Direct",	//coordination type
				{0.1,0.1,0.1,	//atomic coordinates
				 0.12,0.12,0.12,
				 0.08,0.08,0.08},
				{1,1,1,	//if atom can move: mbl
				 0,0,0,
				 0,0,1},
				{0.1,0.1,0.1,	//velocity: vel
				 0.1,0.1,0.1,
				 0.1,0.1,0.1});
	std::unique_ptr<UnitCell> ucell;
	std::ofstream ofs;
	std::string pp_dir;
	std::string output;
	void SetUp()
	{
		ofs.open("running.log");
		GlobalV::relax_new = UTP.relax_new;
		GlobalV::global_out_dir = "./";
		ucell = UTP.SetUcellInfo();
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
	EXPECT_THAT(output,testing::HasSubstr("dft_functional from INPUT does not match that in pseudopot file"));
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

TEST_F(UcellTest,UpdatePosTaud)
{
	double* pos_in = new double[ucell->nat*3];
	ModuleBase::Vector3<double>* tmp = new ModuleBase::Vector3<double>[ucell->nat];
	ucell->set_iat2itia();
	for(int iat=0; iat<ucell->nat; ++iat)
	{
		pos_in[iat*3] = 0.01;
		pos_in[iat*3+1] = 0.01;
		pos_in[iat*3+2] = 0.01;
		int it, ia;
		ucell->iat2iait(iat,&ia,&it);
		tmp[iat] = ucell->atoms[it].taud[ia];
	}
	ucell->update_pos_taud(pos_in);
	for(int iat=0; iat<ucell->nat; ++iat)
	{
		int it, ia;
		ucell->iat2iait(iat,&ia,&it);
		EXPECT_DOUBLE_EQ(ucell->atoms[it].taud[ia].x,tmp[iat].x+0.01);
		EXPECT_DOUBLE_EQ(ucell->atoms[it].taud[ia].y,tmp[iat].y+0.01);
		EXPECT_DOUBLE_EQ(ucell->atoms[it].taud[ia].z,tmp[iat].z+0.01);
	}
	delete[] pos_in;
}

TEST_F(UcellTest,BcastUnitcell2)
{
	if(GlobalV::MY_RANK==0)
	{
		ucell->read_cell_pseudopots(pp_dir,ofs);
	}
	ucell->bcast_unitcell2();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(ucell->atoms[0].ncpp.nbeta,6);
		EXPECT_EQ(ucell->atoms[0].ncpp.nchi,3);
		EXPECT_EQ(ucell->atoms[1].ncpp.nbeta,3);
		EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	}
}

TEST_F(UcellTest,BcastUnitcell)
{
	GlobalV::NSPIN=4;
	ucell->bcast_unitcell();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(ucell->Coordinate,"Direct");
		EXPECT_DOUBLE_EQ(ucell->a1.x,10.0);
		EXPECT_EQ(ucell->atoms[0].na,2);
		EXPECT_EQ(ucell->atoms[1].na,1);
	}
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

TEST_F(UcellTest,CalUx1)
{
	ucell->atoms[0].m_loc_[0].set(0,-1,0);
	ucell->atoms[1].m_loc_[0].set(1,1,1);
	ucell->atoms[1].m_loc_[1].set(0,0,0);
	ucell->cal_ux();
	EXPECT_FALSE(ucell->magnet.lsign_);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[0],0);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[1],-1);
	EXPECT_DOUBLE_EQ(ucell->magnet.ux_[2],0);
}

TEST_F(UcellTest,CalUx2)
{
	ucell->atoms[0].m_loc_[0].set(0,0,0);
	ucell->atoms[1].m_loc_[0].set(1,1,1);
	ucell->atoms[1].m_loc_[1].set(0,0,0);
	//(0,0,0) is also parallel to (1,1,1)
	ucell->cal_ux();
	EXPECT_TRUE(ucell->magnet.lsign_);
	EXPECT_NEAR(ucell->magnet.ux_[0],0.57735,1e-5);
	EXPECT_NEAR(ucell->magnet.ux_[1],0.57735,1e-5);
	EXPECT_NEAR(ucell->magnet.ux_[2],0.57735,1e-5);
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

int main(int argc, char **argv)
{
#ifdef __MPI
	MPI_Init(&argc, &argv);
#endif
	testing::InitGoogleTest(&argc, argv);
	int result = RUN_ALL_TESTS();
#ifdef __MPI
	MPI_Finalize();
#endif
	return result;
}

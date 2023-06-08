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
 *   - UpdatePosTaud
 *     - update_pos_taud(const double* pos)
 *     - bcast_atoms_tau() is also called in the above function, which calls Atom::bcast_atom with many
 *       atomic info in addition to tau
 *   - BcastUnitcell
 *     - bcast basic info of unitcell and basic info of atoms
 *   - BcastUnitcell2
 *     - calls bcast_atoms2() to bcast atomic pseudo info
 *   - ReadPseudo
 *     - read_pseudo()
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

#ifdef __MPI
TEST_F(UcellTest,BcastUnitcell2)
{
	ucell->read_cell_pseudopots(pp_dir,ofs);
	ucell->bcast_unitcell2();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(ucell->atoms[0].ncpp.nbeta,4);
		EXPECT_EQ(ucell->atoms[0].ncpp.nchi,2);
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
		EXPECT_EQ(ucell->atoms[0].na,1);
		EXPECT_EQ(ucell->atoms[1].na,2);
	}
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

TEST_F(UcellTest,ReadPseudo)
{
	GlobalV::global_pseudo_dir = pp_dir;
	GlobalV::out_element_info = 1;
	GlobalV::MIN_DIST_COEF = 0.2;
	ucell->read_pseudo(ofs);
	//check_structure will print some warning info
	//output nonlocal file
	if(GlobalV::MY_RANK==0)
	{
		std::ifstream ifs;
		ifs.open("./C/C.NONLOCAL");
		EXPECT_TRUE(ifs.good());
		ifs.close();
		ifs.open("./H/H.NONLOCAL");
		EXPECT_TRUE(ifs.good());
		ifs.close();
		std::string command1 = "test -d C && rm -rf C";
		std::string command2 = "test -d H && rm -rf H";
		auto error1 = std::system( command1.c_str() );
		EXPECT_EQ(error1,0);
		auto error2 = std::system( command2.c_str() );
		EXPECT_EQ(error2,0);
	}
	//read_cell_pseudopots
	//bcast_unitcell2
	EXPECT_FALSE(ucell->atoms[0].ncpp.has_so);
	EXPECT_FALSE(ucell->atoms[1].ncpp.has_so);
	EXPECT_EQ(ucell->atoms[0].ncpp.nbeta,4);
	EXPECT_EQ(ucell->atoms[0].ncpp.nchi,2);
	EXPECT_EQ(ucell->atoms[1].ncpp.nbeta,3);
	EXPECT_EQ(ucell->atoms[1].ncpp.nchi,1);
	//cal_meshx
	EXPECT_EQ(ucell->meshx,1247);
	//cal_natomwfc
	EXPECT_EQ(ucell->natomwfc,(1+3)*1+1*2);
	//cal_nwfc
	EXPECT_EQ(ucell->lmax,2);
	EXPECT_EQ(ucell->lmax_ppwf,1);
}

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

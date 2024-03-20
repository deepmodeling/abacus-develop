#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include<streambuf>
#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of Atom_pseudo
 ***********************************************/

/**
 * - Tested Functions:
 *   - Atom_pseudo
 *     - constructor of class Atom_pseudo
 *   - ~Atom_pseudo
 *     - deconstructor of class Atom_pseudo
 *   - set_d_so
 *     - set spin-orbital info from pseudopotential
 *   - bcast_atom_pseudo
 *     - bcast upf201 pp info to other processes
 */

#define private public
#include "module_cell/read_pp.h"
#include "module_cell/pseudo.h"
#include "module_cell/atom_pseudo.h"

class AtomPseudoTest : public testing::Test
{
protected:
	std::unique_ptr<Pseudopot_upf> upf{new Pseudopot_upf};
	std::unique_ptr<Atom_pseudo> atom_pseudo{new Atom_pseudo};
};

TEST_F(AtomPseudoTest, SetDSo)
{
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
	std::ifstream ifs;
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs);
	atom_pseudo->set_pseudo(*upf);
	ifs.close();
	EXPECT_EQ(atom_pseudo->nh,14);
	EXPECT_TRUE(atom_pseudo->has_so);
	ModuleBase::ComplexMatrix d_so_in(atom_pseudo->nh*2,atom_pseudo->nh*2);
	int nproj = 6;
	int nproj_soc = 4;
	bool has_so = 1;
	GlobalV::NSPIN = 4;
	atom_pseudo->set_d_so(d_so_in,nproj,nproj_soc,has_so);
	EXPECT_NEAR(atom_pseudo->d_so(0,0,0).real(),1e-8,1e-7);
	EXPECT_NEAR(atom_pseudo->d_so(0,0,0).imag(),1e-8,1e-7);
	GlobalV::LSPINORB = 1;
	atom_pseudo->set_d_so(d_so_in,nproj,nproj_soc,has_so);
	EXPECT_NEAR(atom_pseudo->d_so(0,0,0).real(),1e-8,1e-7);
	EXPECT_NEAR(atom_pseudo->d_so(0,0,0).imag(),1e-8,1e-7);
#ifdef __MPI
	}
#endif
}

#ifdef __MPI
TEST_F(AtomPseudoTest, BcastAtomPseudo)
{
	if(GlobalV::MY_RANK==0)
	{
		std::ifstream ifs;
		ifs.open("./support/C.upf");
		GlobalV::PSEUDORCUT = 15.0;
		upf->read_pseudo_upf201(ifs);
		atom_pseudo->set_pseudo(*upf);
		ifs.close();
	}
	atom_pseudo->bcast_atom_pseudo();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(atom_pseudo->nbeta,6);
		EXPECT_EQ(atom_pseudo->nchi,3);
		EXPECT_DOUBLE_EQ(atom_pseudo->rho_atc[0],8.7234550809E-01);
	}
}

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
#undef private

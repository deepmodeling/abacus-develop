#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/write_dos_pw.h"
#include "module_io/input.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "for_testing_klist.h"
#include "dos_test.h"

/************************************************
 *  unit test of write_dos_pw
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_dos_pw()
 *     - the function to calculate and print out
 *     - density of states in pw basis calculation
 */

Input INPUT;

class DosPWTest : public ::testing::Test
{
protected:
	K_Vectors* kv = nullptr;
	ModuleBase::matrix ekb;
	ModuleBase::matrix wg;
	void SetUp()
	{
		kv = new K_Vectors;
	}
	void TearDown()
	{
		delete kv;
	}
};

TEST_F(DosPWTest,Dos1)
{
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	GlobalV::NSPIN = 1;
	INPUT.dos_emax_ev = dosp.emax_ev;
	INPUT.dos_setemax = true;
	INPUT.dos_emin_ev = dosp.emin_ev;
	INPUT.dos_setemin = true;
	kv->nks = dosp.nks;
	kv->nkstot = dosp.nkstot;
	kv->isk.reserve(kv->nks);
	kv->wk.reserve(kv->nks);
	for(int ik=0; ik<kv->nks; ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}
	GlobalV::NBANDS = dosp.nbands;
	ModuleIO::write_dos_pw(dosp.ekb,
			dosp.wg,
			*kv,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff);
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("DOS1_smearing.dat");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("             3200")); // number of electrons is 32
		ifs.close();
		ifs.open("DOS1");
		std::string str1((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str1, testing::HasSubstr("4800")); //number of energy points is (18-(-6))/0.005
		ifs.close();
		remove("DOS1_smearing.dat");
		remove("DOS1");
#ifdef __MPI
	}
#endif
}

TEST_F(DosPWTest,Dos2)
{
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	GlobalV::NSPIN = 1;
	INPUT.dos_emax_ev = dosp.emax_ev;
	INPUT.dos_setemax = false;
	INPUT.dos_emin_ev = dosp.emin_ev;
	INPUT.dos_setemin = false;
	kv->nks = dosp.nks;
	kv->nkstot = dosp.nkstot;
	kv->isk.reserve(kv->nks);
	kv->wk.reserve(kv->nks);
	for(int ik=0; ik<kv->nks; ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}
	GlobalV::NBANDS = dosp.nbands;
	ModuleIO::write_dos_pw(dosp.ekb,
			dosp.wg,
			*kv,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff);
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("DOS1_smearing.dat");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("             3197")); // number of electrons is 32
		ifs.close();
		ifs.open("DOS1");
		std::string str1((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str1, testing::HasSubstr("4531")); //number of energy points
		ifs.close();
		remove("DOS1_smearing.dat");
		remove("DOS1");
#ifdef __MPI
	}
#endif
}

#ifdef __MPI
int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv);

	testing::InitGoogleTest(&argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif

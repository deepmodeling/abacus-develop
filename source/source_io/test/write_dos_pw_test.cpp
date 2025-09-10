#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_io/write_dos_pw.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "for_testing_klist.h"
#include "dos_test.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

/************************************************
 *  unit test of write_dos_pw
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_dos_pw()
 *     - the function to calculate and print out
 *     - density of states in pw basis calculation
 */


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
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	PARAM.input.nspin = 1;
	PARAM.input.dos_emax_ev = dosp.emax_ev;
	PARAM.sys.dos_setemax = true;
	PARAM.input.dos_emin_ev = dosp.emin_ev;
	PARAM.sys.dos_setemin = true;
	kv->set_nks(dosp.nks);
	kv->set_nkstot(dosp.nkstot);
	kv->isk.reserve(kv->get_nks());
	kv->wk.reserve(kv->get_nks());
	for(int ik=0; ik<kv->get_nks(); ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}
	PARAM.input.nbands = dosp.nbands;

    // initialize the Fermi energy
    elecstate::efermi fermi_energy;

	std::ofstream ofs("write_dos_pw.log");
    
    UnitCell ucell;

	ModuleIO::write_dos_pw(
            ucell, // this should be unitcell, 2025-04-12
            dosp.ekb,
			dosp.wg,
			*kv,
			PARAM.inp.nbands,
            -1, // istep_in
			fermi_energy,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff,
			ofs);
    ofs.close();
    remove("write_dos_pw.log");

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("dos.txt");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("4801 # number of points"));
        EXPECT_THAT(str, testing::HasSubstr("           -4.6           0.25        0.28125        1.42515       0.159819"));
        EXPECT_THAT(str, testing::HasSubstr("             18              0             16              0             16"));
		ifs.close();
		remove("dos.txt");
#ifdef __MPI
	}
#endif
}


TEST_F(DosPWTest,Dos2)
{
    //is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	PARAM.input.nspin = 1;
	PARAM.input.dos_emax_ev = dosp.emax_ev;
	PARAM.sys.dos_setemax = false;
	PARAM.input.dos_emin_ev = dosp.emin_ev;
	PARAM.sys.dos_setemin = false;
	kv->set_nks(dosp.nks);
	kv->set_nkstot(dosp.nkstot);
	kv->isk.reserve(kv->get_nks());
	kv->wk.reserve(kv->get_nks());
	for(int ik=0; ik<kv->get_nks(); ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}
	PARAM.input.nbands = dosp.nbands;

    // initialize the Fermi energy
    elecstate::efermi fermi_energy;

	std::ofstream ofs("write_dos_pw.log");

    UnitCell ucell;

	ModuleIO::write_dos_pw(
			ucell,
			dosp.ekb,
			dosp.wg,
			*kv,
			PARAM.inp.nbands,
			-1, // istep_in
            fermi_energy,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff,
            ofs);
    ofs.close();
    remove("write_dos_pw.log");

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("dos.txt");
		std::string str1((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str1, testing::HasSubstr("4532 # number of points"));
        EXPECT_THAT(str1, testing::HasSubstr("       -5.38811        0.03125        0.03125"));
        EXPECT_THAT(str1, testing::HasSubstr("        3.07189         0.1875        5.46875"));
		ifs.close();
		remove("dos.txt");
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
    
    // only test the second one
    // ::testing::GTEST_FLAG(filter) = "DosPWTest.Dos2";

	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif

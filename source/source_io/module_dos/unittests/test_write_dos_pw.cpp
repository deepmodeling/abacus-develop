#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../write_dos_pw.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "./for_testing_klist.h"
#include "./dos_test.h"

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
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,dos_test_data::nks,dos_test_data::nkstot,dos_test_data::nbands);
	dosp.set_isk();
	dosp.set_wk();
	dosp.set_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	kv->set_nks(dosp.nks);
	kv->set_nkstot(dosp.nkstot);
	kv->isk.reserve(kv->get_nks());
	kv->wk.reserve(kv->get_nks());
	for(int ik=0; ik<kv->get_nks(); ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}

    // initialize the Fermi energy
    elecstate::Efermi fermi_energy;

	std::ofstream ofs("write_dos_pw.log");

    UnitCell ucell;

    const int nspin = 1;
    const int out_dos = 1;
    const bool dos_setemax = true;
    const bool dos_setemin = true;
    const bool two_fermi = false;
    const bool out_app_flag = false;
    const int bndpar = 1;
    const std::string global_out_dir = "./";

	ModuleIO::write_dos_pw(
            ucell, // this should be unitcell, 2025-04-12
            dosp.ekb,
			dosp.wg,
			*kv,
			dosp.nbands,
            -1, // istep_in
			fermi_energy,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff,
			nspin,
			out_dos,
			dos_setemax,
			dosp.emax_ev,
			dos_setemin,
			dosp.emin_ev,
			two_fermi,
			out_app_flag,
			bndpar,
			global_out_dir,
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
		EXPECT_THAT(str, testing::HasSubstr("   -4.995000    0.500000    0.500000"));
		EXPECT_THAT(str, testing::HasSubstr("   18.000000    0.000000    8.000000"));
		ifs.close();
		remove("dos.txt");
#ifdef __MPI
	}
#endif
}


TEST_F(DosPWTest,Dos2)
{
    //is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,dos_test_data::nks,dos_test_data::nkstot,dos_test_data::nbands);
	dosp.set_isk();
	dosp.set_wk();
	dosp.set_istate_info();
	EXPECT_EQ(dosp.is,0);
	double dos_scale = 0.01;
	kv->set_nks(dosp.nks);
	kv->set_nkstot(dosp.nkstot);
	kv->isk.reserve(kv->get_nks());
	kv->wk.reserve(kv->get_nks());
	for(int ik=0; ik<kv->get_nks(); ++ik)
	{
		kv->isk[ik] = dosp.isk[ik];
		kv->wk[ik] = dosp.wk[ik];
	}

    // initialize the Fermi energy
    elecstate::Efermi fermi_energy;

	std::ofstream ofs("write_dos_pw.log");

    UnitCell ucell;

    const int nspin = 1;
    const int out_dos = 1;
    const bool two_fermi = false;
    const bool out_app_flag = false;
    const int bndpar = 1;
    const std::string global_out_dir = "./";
    const bool dos_setemax = false;
    const bool dos_setemin = false;

	ModuleIO::write_dos_pw(
			ucell,
			dosp.ekb,
			dosp.wg,
			*kv,
			dosp.nbands,
			-1, // istep_in
            fermi_energy,
			dosp.de_ev,
			dos_scale,
			dosp.bcoeff,
			nspin,
			out_dos,
			dos_setemax,
			dosp.emax_ev,
			dos_setemin,
			dosp.emin_ev,
			two_fermi,
			out_app_flag,
			bndpar,
			global_out_dir,
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
		EXPECT_THAT(str1, testing::HasSubstr("2526 # number of points"));
		EXPECT_THAT(str1, testing::HasSubstr("   -4.997500    0.500000    0.500000"));
		EXPECT_THAT(str1, testing::HasSubstr("    7.567500    0.000000    8.000000"));
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

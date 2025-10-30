#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_io/cal_dos.h"
#include "source_base/global_variable.h"
#include <string>
#ifdef __MPI
#include "mpi.h"
#endif
#include "dos_test.h"

/************************************************
 *  unit test of ca_dos
 ***********************************************/

/**
 * - Tested Functions:
 *   - cal_dos()
 *     - the function to calculate and print out
 *     - density of states
 */


class DosTest : public ::testing::Test
{
protected:
	std::string output;
};

TEST_F(DosTest,Dos)
{
	//is,fa,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);

    const int istep = 1;
	ModuleIO::cal_dos(dosp.is,
			dosp.fa,
			dosp.de_ev,
			dosp.emax_ev,
			dosp.emin_ev,
			dosp.bcoeff,
			dosp.nks,
			dosp.nkstot,
			dosp.wk,
			dosp.isk,
			dosp.nbands,
			dosp.ekb,
			dosp.wg,
			istep);

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open(dosp.fa.c_str());
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("4801 # number of points"));
        EXPECT_THAT(str, testing::HasSubstr("          -5.39        0.03125        0.03125       0.178099      0.0160702"));
        EXPECT_THAT(str, testing::HasSubstr("           3.07         0.1875        5.46875        1.07003        5.37765"));
		ifs.close();
		remove("doss1_pw.txt");
#ifdef __MPI
	}
#endif
}



TEST_F(DosTest,DosW1)
{
	//is,fa,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",-0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	EXPECT_LE(dosp.de_ev,0);
	GlobalV::ofs_warning.open("warning1.log");

    const int istep = 1;
	EXPECT_NO_THROW(ModuleIO::cal_dos(dosp.is,
			dosp.fa,
			dosp.de_ev,
			dosp.emax_ev,
			dosp.emin_ev,
			dosp.bcoeff,
			dosp.nks,
			dosp.nkstot,
			dosp.wk,
			dosp.isk,
			dosp.nbands,
			dosp.ekb,
			dosp.wg,
            istep));
	GlobalV::ofs_warning.close();
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("warning1.log");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("ModuleIO::cal_dos  warning : de <= 0"));
		ifs.close();
		remove("warning1.log");
		remove("doss1_pw.txt");
#ifdef __MPI
	}
#endif
}


TEST_F(DosTest,DosW2)
{
    //is,fa,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,-6,18,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	GlobalV::ofs_warning.open("warning2.log");

    const int istep = 1;
	EXPECT_NO_THROW(ModuleIO::cal_dos(dosp.is,
			dosp.fa,
			dosp.de_ev,
			dosp.emax_ev,
			dosp.emin_ev,
			dosp.bcoeff,
			dosp.nks,
			dosp.nkstot,
			dosp.wk,
			dosp.isk,
			dosp.nbands,
			dosp.ekb,
			dosp.wg,
			istep));

	GlobalV::ofs_warning.close();
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("warning2.log");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("ModuleIO::cal_dos  warning : emax_ev < emin_ev"));
		ifs.close();
		remove("warning2.log");
		remove("doss1_pw.txt");
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

    // only test a certain one
    //::testing::GTEST_FLAG(filter) = "DosTest.DosW1";

	int result = RUN_ALL_TESTS();

	MPI_Finalize();

	return result;
}
#endif

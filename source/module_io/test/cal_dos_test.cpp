#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/cal_dos.h"
#include "module_base/global_variable.h"
#include <string>
#ifdef __MPI
#include "mpi.h"
#endif
#include "dos_test.h"

/************************************************
 *  unit test of calculate_dos
 ***********************************************/

/**
 * - Tested Functions:
 *   - calculate_dos()
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
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	ModuleIO::calculate_dos(dosp.is,
			dosp.fa,
			dosp.fa1,
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
			dosp.wg);
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

TEST_F(DosTest,DosW1)
{
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",-0.005,18,-6,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	EXPECT_LE(dosp.de_ev,0);
	GlobalV::ofs_warning.open("warninglog1");
	EXPECT_NO_THROW(ModuleIO::calculate_dos(dosp.is,
			dosp.fa,
			dosp.fa1,
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
			dosp.wg));
	GlobalV::ofs_warning.close();
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("warninglog1");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("ModuleIO::calculate_dos  warning : de <= 0"));
		ifs.close();
		remove("warninglog1");
		remove("DOS1_smearing.dat");
		remove("DOS1");
#ifdef __MPI
	}
#endif
}

TEST_F(DosTest,DosW2)
{
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",0.005,-6,18,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	GlobalV::ofs_warning.open("warninglog2");
	EXPECT_NO_THROW(ModuleIO::calculate_dos(dosp.is,
			dosp.fa,
			dosp.fa1,
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
			dosp.wg));
	GlobalV::ofs_warning.close();
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("warninglog2");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("ModuleIO::calculate_dos  warning : emax_ev < emin_ev"));
		ifs.close();
		remove("warninglog2");
		remove("DOS1_smearing.dat");
		remove("DOS1");
#ifdef __MPI
	}
#endif
}

TEST_F(DosTest,DosW3)
{
				//is,fa,fa1,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"DOS1","DOS1_smearing.dat",0.005,18,18,0.07,36,36,8);
	dosp.set_isk();
	dosp.read_wk();
	dosp.read_istate_info();
	EXPECT_EQ(dosp.is,0);
	GlobalV::ofs_warning.open("warninglog3");
	EXPECT_NO_THROW(ModuleIO::calculate_dos(dosp.is,
			dosp.fa,
			dosp.fa1,
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
			dosp.wg));
	GlobalV::ofs_warning.close();
#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open("warninglog3");
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("ModuleIO::calculate_dos  warning : npoints <= 0"));
		ifs.close();
		remove("warninglog3");
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

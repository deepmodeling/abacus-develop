#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../cal_dos.h"
#include "source_base/global_variable.h"
#include <string>
#include <sstream>
#ifdef __MPI
#include "mpi.h"
#endif
#include "./dos_test.h"

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
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,18,-6,0.07,dos_test_data::nks,dos_test_data::nkstot,dos_test_data::nbands);
	dosp.set_isk();
	dosp.set_wk();
	dosp.set_istate_info();
	EXPECT_EQ(dosp.is,0);

    const int istep = 1;
    const bool out_app_flag = false;
    const int bndpar = 1;
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
			istep,
			out_app_flag,
			bndpar);

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
#endif
		std::ifstream ifs;
		ifs.open(dosp.fa.c_str());
		std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
		EXPECT_THAT(str, testing::HasSubstr("4801 # number of points"));
	EXPECT_THAT(str, testing::HasSubstr("   -5.000000     0.500000     0.500000"));
	EXPECT_THAT(str, testing::HasSubstr("   18.000000     0.000000     8.000000"));
	// node-sampled smeared DOS: Gaussian evaluated at the exact eigenvalue
	// position, e.g. node -5.000 eV gets 0.5*G(0) = 0.5/(sqrt(2*pi)*0.07)
	EXPECT_THAT(str, testing::HasSubstr("2.849588"));
    std::istringstream rows(str);
    std::string line;
    std::getline(rows, line); // ionic step
    std::getline(rows, line); // declared grid size
    int declared_points = 0;
    std::istringstream(line) >> declared_points;
    EXPECT_EQ(declared_points, 4801);
    std::getline(rows, line); // column names
    int row_count = 0;
    while (std::getline(rows, line))
    {
        std::istringstream fields(line);
        double energy = 0.0;
        double raw = 0.0;
        double integral = 0.0;
        double smeared = 0.0;
        double smeared_integral = 0.0;
        ASSERT_TRUE(fields >> energy >> raw >> integral >> smeared >> smeared_integral);
        EXPECT_NEAR(energy, -6.0 + row_count * 0.005, 5.0e-7);
        if (row_count == 200)
        {
            EXPECT_DOUBLE_EQ(energy, -5.0);
            EXPECT_DOUBLE_EQ(raw, 0.5);
            EXPECT_DOUBLE_EQ(integral, 0.5);
            EXPECT_NEAR(smeared, 2.849588, 5.0e-7);
        }
        if (row_count == 4800)
        {
            EXPECT_DOUBLE_EQ(energy, 18.0);
            EXPECT_DOUBLE_EQ(integral, 8.0);
        }
        ++row_count;
    }
    EXPECT_EQ(row_count, declared_points);
	ifs.close();
	remove("doss1_pw.txt");
#ifdef __MPI
	}
#endif
}



TEST_F(DosTest,DosW1)
{
	//is,fa,de_ev,emax_ev,emin_ev,bcoeff,nks,nkstot,nbands
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",-0.005,18,-6,0.07,dos_test_data::nks,dos_test_data::nkstot,dos_test_data::nbands);
	dosp.set_isk();
	dosp.set_wk();
	dosp.set_istate_info();
	EXPECT_EQ(dosp.is,0);
	EXPECT_LE(dosp.de_ev,0);
	GlobalV::ofs_warning.open("warning1.log");

    const int istep = 1;
    const bool out_app_flag = false;
    const int bndpar = 1;
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
			istep,
			out_app_flag,
			bndpar));
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
	DosPrepare dosp = DosPrepare(0,"doss1_pw.txt",0.005,-6,18,0.07,dos_test_data::nks,dos_test_data::nkstot,dos_test_data::nbands);
	dosp.set_isk();
	dosp.set_wk();
	dosp.set_istate_info();
	EXPECT_EQ(dosp.is,0);
	GlobalV::ofs_warning.open("warning2.log");

    const int istep = 1;
    const bool out_app_flag = false;
    const int bndpar = 1;
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
			istep,
			out_app_flag,
			bndpar));

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

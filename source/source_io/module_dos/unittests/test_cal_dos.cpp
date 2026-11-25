#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../cal_dos.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>
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
    const int out_app_flag = 0;
    const int precision = 6;
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
			precision,
			bndpar, nullptr);

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
    // Skip the per-frame provenance comments, then the ionic-step marker.
    while (std::getline(rows, line))
    {
        if (line.find("# ABACUS version:") != std::string::npos)
        {
            continue;
        }
        if (line.find("# Written at ") != std::string::npos)
        {
            continue;
        }
        break;
    }
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



TEST_F(DosTest, CapturesSpectrumWithoutTextOutput)
{
    const pid_t pid = getpid();
    const std::string text_path = testing::TempDir() + "abacus_cal_dos_capture_" + std::to_string(pid) + ".txt";
    const double step = 0.5;
    const double sigma = 0.5;
    ModuleBase::matrix eigenvalues(1, 1);
    eigenvalues(0, 0) = 0.0;
    ModuleBase::matrix occupations(1, 1);
    occupations(0, 0) = 1.0;
    ModuleIO::DosSpectrum spectrum;
    const bool calculated = ModuleIO::cal_dos(0, text_path, step, 0.5, 0.0, sigma,
                                              1, 1, {1.0}, {0}, 1, eigenvalues, occupations,
                                              /*istep=*/-1, /*out_app_flag=*/0, /*precision=*/6, /*bndpar=*/1, &spectrum);
    ASSERT_TRUE(calculated);
    const std::ifstream text_output(text_path);
    EXPECT_FALSE(text_output.is_open());

    // One state at 0 eV, sampled at the Gaussian peak and one sigma above it.
    const double peak = 1.0 / (sigma * std::sqrt(ModuleBase::TWO_PI));
    const double shoulder = peak * std::exp(-0.5);
    const ModuleIO::DosSpectrum expected = {
        {0.0, 1.0, 1.0, peak, step * peak},
        {0.5, 0.0, 1.0, shoulder, step * (peak + shoulder)}
    };
    ASSERT_EQ(spectrum.size(), expected.size());
    for (size_t row = 0; row < expected.size(); ++row)
    {
        for (std::size_t col = 0; col < expected[row].size(); ++col)
        {
            EXPECT_NEAR(spectrum[row][col], expected[row][col], 1e-12) << "row " << row << ", column " << col;
        }
    }
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
    const int out_app_flag = 0;
    const int precision = 6;
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
			precision,
			bndpar, nullptr));
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
    const int out_app_flag = 0;
    const int precision = 6;
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
			precision,
			bndpar, nullptr));

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

// ---------------------------------------------------------------------------
// out_app_flag = 1: cal_dos appends successive records to the same file.
// ---------------------------------------------------------------------------
TEST_F(DosTest, Mode1Appends)
{
    const std::string fname = "doss1_mode1.txt";
    DosPrepare dosp(0, fname, 0.005, 18, -6, 0.07,
                    dos_test_data::nks, dos_test_data::nkstot, dos_test_data::nbands);
    dosp.set_isk();
    dosp.set_wk();
    dosp.set_istate_info();

    const int out_app_flag = 1;
    const int precision = 6;
    const int bndpar = 1;

    // First call: istep=0 -> "# Ionic_Step 1".
    ModuleIO::cal_dos(dosp.is, dosp.fa, dosp.de_ev, dosp.emax_ev, dosp.emin_ev,
                      dosp.bcoeff, dosp.nks, dosp.nkstot, dosp.wk, dosp.isk,
                      dosp.nbands, dosp.ekb, dosp.wg,
                      /*istep=*/0, out_app_flag, precision, bndpar, nullptr);
    // Second call: istep=1 -> "# Ionic_Step 2", appended.
    ModuleIO::cal_dos(dosp.is, dosp.fa, dosp.de_ev, dosp.emax_ev, dosp.emin_ev,
                      dosp.bcoeff, dosp.nks, dosp.nkstot, dosp.wk, dosp.isk,
                      dosp.nbands, dosp.ekb, dosp.wg,
                      /*istep=*/1, out_app_flag, precision, bndpar, nullptr);

#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
#endif
    {
        std::ifstream ifs(fname);
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        ifs.close();
        // Two Ionic_Step headers present.
        EXPECT_THAT(str, testing::HasSubstr("# Ionic_Step 1"));
        EXPECT_THAT(str, testing::HasSubstr("# Ionic_Step 2"));
        remove(fname.c_str());
    }
}

// ---------------------------------------------------------------------------
// out_app_flag = 2: same overwrite behaviour as mode 0 (file is truncated each
// call; naming is handled by the caller write_dos_* which adds g{N}).
// ---------------------------------------------------------------------------
TEST_F(DosTest, Mode2Overwrites)
{
    const std::string fname = "doss1_mode2.txt";
    DosPrepare dosp(0, fname, 0.005, 18, -6, 0.07,
                    dos_test_data::nks, dos_test_data::nkstot, dos_test_data::nbands);
    dosp.set_isk();
    dosp.set_wk();
    dosp.set_istate_info();

    const int out_app_flag = 2;
    const int precision = 6;
    const int bndpar = 1;

    // First call writes istep 1, second overwrites with istep 2.
    ModuleIO::cal_dos(dosp.is, dosp.fa, dosp.de_ev, dosp.emax_ev, dosp.emin_ev,
                      dosp.bcoeff, dosp.nks, dosp.nkstot, dosp.wk, dosp.isk,
                      dosp.nbands, dosp.ekb, dosp.wg,
                      /*istep=*/0, out_app_flag, precision, bndpar, nullptr);
    ModuleIO::cal_dos(dosp.is, dosp.fa, dosp.de_ev, dosp.emax_ev, dosp.emin_ev,
                      dosp.bcoeff, dosp.nks, dosp.nkstot, dosp.wk, dosp.isk,
                      dosp.nbands, dosp.ekb, dosp.wg,
                      /*istep=*/1, out_app_flag, precision, bndpar, nullptr);

#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
#endif
    {
        std::ifstream ifs(fname);
        std::string str((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
        ifs.close();
        // Only the last Ionic_Step header must be present.
        EXPECT_THAT(str, testing::HasSubstr("# Ionic_Step 2"));
        // A simple count check: "# Ionic_Step 1" must appear only inside
        // the DOS grid lines as a value, not as a header. Check that the
        // header for step 1 is absent by verifying it is NOT found as a
        // standalone header line.
        EXPECT_THAT(str, testing::Not(testing::HasSubstr("\n# Ionic_Step 1")));
        remove(fname.c_str());
    }
}

// ---------------------------------------------------------------------------
// The precision argument controls the number of decimals in the DOS table.
// ---------------------------------------------------------------------------
TEST_F(DosTest, PrecisionControlsDecimals)
{
    const std::string fname = testing::TempDir() + "doss1_precision.txt";
    DosPrepare dosp(0, fname, 0.005, 18, -6, 0.07,
                    dos_test_data::nks, dos_test_data::nkstot, dos_test_data::nbands);
    dosp.set_isk();
    dosp.set_wk();
    dosp.set_istate_info();

    const int istep = 0;
    const int out_app_flag = 0;
    const int precision = 4;
    const int bndpar = 1;

    ModuleIO::cal_dos(dosp.is, dosp.fa, dosp.de_ev, dosp.emax_ev, dosp.emin_ev,
                      dosp.bcoeff, dosp.nks, dosp.nkstot, dosp.wk, dosp.isk,
                      dosp.nbands, dosp.ekb, dosp.wg,
                      istep, out_app_flag, precision, bndpar, nullptr);

    std::ifstream ifs(fname);
    ASSERT_TRUE(ifs.is_open());
    int rows = 0;
    std::string line;
    while (std::getline(ifs, line))
    {
        std::vector<std::string> fields;
        std::istringstream ls(line);
        for (std::string token; ls >> token;)
        {
            fields.push_back(token);
        }
        if (line.find('#') != std::string::npos || fields.size() != 5)
        {
            continue; // provenance, column headers and the npoints line
        }
        ++rows;
        // Every field must carry exactly four decimals; a substring check would
        // also accept five-decimal fields such as -5.00000 or 0.50000.
        for (const std::string& field : fields)
        {
            EXPECT_THAT(field, testing::MatchesRegex(R"(-?[0-9]+\.[0-9]{4})")) << line;
        }
    }
    ifs.close();
    remove(fname.c_str());
    EXPECT_GT(rows, 0);
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

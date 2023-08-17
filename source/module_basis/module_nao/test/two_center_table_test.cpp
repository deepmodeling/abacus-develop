#include "module_basis/module_nao/two_center_table.h"

#include "gtest/gtest.h"
#include "module_base/constants.h"
#include "module_base/spherical_bessel_transformer.h"

#include "module_base/math_integral.h"

#include <cstring>
#include <chrono>
using iclock = std::chrono::high_resolution_clock;

//////////////////////////////////////////
//! for comparison with module_ao
#include "module_basis/module_ao/ORB_table_phi.h"
//////////////////////////////////////////

#ifdef __MPI
#include <mpi.h>
#endif

/***********************************************************
 *      Unit test of class "TwoCenterTable"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - builds a two-center integral radial table from two RadialCollection objects
 *                                                                      */
class TwoCenterTableTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();

    TwoCenterTable S_tab;
    TwoCenterTable T_tab;

    TwoCenterTable betapsi_tab;

    RadialCollection orb;
    int nfile = 0;                                              //! number of orbital files
    std::string* file = nullptr;                                //!< orbital files to read from
    std::string log_file = "./test_files/two_center_table.log"; //!< file for logging

    double tol = 1e-6;       /// tolerance for comparison between new and legacy tables
    double tol_d_abs = 1e-5; /// absolute tolerance for derivative table (compared to finite difference)
    double tol_d_rel = 1e-2; /// relative tolerance for derivative table (compared to finite difference)
};

void TwoCenterTableTest::SetUp()
{
#ifdef __MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    std::string dir = "../../../../../tests/PP_ORB/";
    nfile = 8;
    file = new std::string[nfile];
    file[0] = dir + "C_gga_8au_100Ry_2s2p1d.orb";
    file[1] = dir + "Fe_gga_9au_100Ry_4s2p2d1f.orb";
    file[2] = dir + "O_gga_10au_100Ry_2s2p1d.orb";
    file[3] = dir + "H_gga_8au_60Ry_2s1p.orb";
    file[4] = dir + "Cs_gga_10au_100Ry_4s2p1d.orb";
    file[5] = dir + "Pb_gga_7au_100Ry_2s2p2d1f.orb";
    file[6] = dir + "F_gga_7au_100Ry_2s2p1d.orb";
    file[7] = dir + "I_gga_7au_100Ry_2s2p2d1f.orb";
}

void TwoCenterTableTest::TearDown()
{
    delete[] file;
}

TEST_F(TwoCenterTableTest, BuildOverlapAndKinetic)
{
    orb.build(nfile, file, 'o');
    double rmax = orb.rcut_max() * 2.0;
    double dr = 0.01;
    int nr = static_cast<int>(rmax / dr) + 1;
    
    //ModuleBase::SphericalBesselTransformer sbt;
    //sbt.set_fftw_plan_flag(FFTW_MEASURE); // not necessarily worth it!
    //orb.set_transformer(&sbt, 0);
    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    iclock::time_point start;
    std::chrono::duration<double> dur;

    start = iclock::now();
    S_tab.build(orb, orb, 'S', nr, rmax, true);
    T_tab.build(orb, orb, 'T', nr, rmax, true);
    dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " s" << std::endl;

    // check whether the derivative table agrees with the finite difference of table
    int ntype = nfile;
	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		for (int T2 = T1 ; T2 < ntype ; T2++)
		{
			for (int L1 = 0; L1 <= orb(T1).lmax(); L1++)
			{
				for (int N1 = 0; N1 < orb(T1).nzeta(L1); N1++)
				{
					for (int L2 = 0; L2 <= orb(T2).lmax(); L2++)
					{
						for (int N2 = 0; N2 < orb(T2).nzeta(L2); N2++)
						{
							for (int L = std::abs(L1-L2); L <= (L1+L2) ; L += 2)
                            {
                                const double* f = S_tab.table(T1, L1, N1, T2, L2, N2, L, false);
                                const double* df  = S_tab.table(T1, L1, N1, T2, L2, N2, L, true);

                                for (int ir = 4; ir != S_tab.nr() - 4; ++ir)
                                {
                                    double df_fd = ( - 1.0/280 *(f[ir+4] - f[ir-4])
                                                     + 4.0/105 *(f[ir+3] - f[ir-3])
                                                     - 0.2 * (f[ir+2] - f[ir-2])
                                                     + 0.8 * (f[ir+1] - f[ir-1]) ) / dr;

                                    // df is (d/dR)(S/R^l), it should be rescaled to have the
                                    // same unit as dS/dR in order to have meaningful error comparison
                                    double err_abs = std::abs(df_fd - df[ir]) * std::pow(ir * dr, L);
                                    double err_rel = std::abs((df_fd - df[ir]) / df[ir]);
                                    if ( err_abs > tol_d_abs && err_rel > tol_d_rel )
                                    {
                                        printf("T1 = %i   L1 = %i   N1 = %i   T2 = %i   L2 = %i   N2 = %i   L = %i   ir = %2i   df_fd = % 8.5e   df_tab = % 8.5e   err_abs = %8.5e   err_rel = %8.5e\n",
                                                T1, L1, N1, T2, L2, N2, L, ir, df_fd, df[ir], err_abs, err_rel);
                                    }

                                    EXPECT_TRUE( err_abs < tol_d_abs || err_rel < tol_d_rel );
                                }

                                f  = T_tab.table(T1, L1, N1, T2, L2, N2, L, false);
                                df = T_tab.table(T1, L1, N1, T2, L2, N2, L, true);

                                for (int ir = 4; ir != T_tab.nr() - 4; ++ir)
                                {
                                    double df_fd = ( - 1.0/280 *(f[ir+4] - f[ir-4])
                                                     + 4.0/105 *(f[ir+3] - f[ir-3])
                                                     - 0.2 * (f[ir+2] - f[ir-2])
                                                     + 0.8 * (f[ir+1] - f[ir-1]) ) / dr;

                                    // df is (d/dR)(S/R^l), it should be rescaled to have the
                                    // same unit as dS/dR in order to have meaningful error comparison
                                    double err_abs = std::abs(df_fd - df[ir]) * std::pow(ir * dr, L);
                                    double err_rel = std::abs((df_fd - df[ir]) / df[ir]);
                                    if ( err_abs > tol_d_abs && err_rel > tol_d_rel )
                                    {
                                        printf("T1 = %i   L1 = %i   N1 = %i   T2 = %i   L2 = %i   N2 = %i   L = %i   ir = %2i   df_fd = % 8.5e   df_tab = % 8.5e   err_abs = %8.5e   err_rel = %8.5e\n",
                                                T1, L1, N1, T2, L2, N2, L, ir, df_fd, df[ir], err_abs, err_rel);
                                    }

                                    EXPECT_TRUE( err_abs < tol_d_abs || err_rel < tol_d_rel );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    EXPECT_EQ(S_tab.with_deriv(), true);
    EXPECT_EQ(T_tab.with_deriv(), true);

    EXPECT_EQ(S_tab.nr(), nr);
    EXPECT_EQ(T_tab.nr(), nr);

    EXPECT_EQ(S_tab.rmax(), rmax);
    EXPECT_EQ(T_tab.rmax(), rmax);
}

TEST_F(TwoCenterTableTest, LegacyConsistency)
{
    // use less files so that the test wouldn't take too long
    nfile = 3;

    orb.build(nfile, file, 'o');
    double rmax = orb.rcut_max() * 2.0;
    double dr = 0.01;
    int nr = static_cast<int>(rmax / dr) + 1;
    
    //ModuleBase::SphericalBesselTransformer sbt;
    //orb.set_transformer(&sbt, 0);
    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    iclock::time_point start;
    std::chrono::duration<double> dur;

    start = iclock::now();
    S_tab.build(orb, orb, 'S', nr, rmax, true);
    T_tab.build(orb, orb, 'T', nr, rmax, true);
    dur = iclock::now() - start;
    std::cout << "radfft time elapsed = " << dur.count() << " s" << std::endl;

    //////////////////////////////////////////////
    //      module_ao ORB_table_phi
    //////////////////////////////////////////////

    ORB_table_phi otp;
    LCAO_Orbitals lcao;
    ModuleBase::Sph_Bessel_Recursive::D2 sbr;

    // prepare data required to initialized ORB_table_phi
    std::ofstream ofs_log;
    int ntype;
    int lmax;
    int out_mat_r; // unused variable
    bool force_flag;
    int my_rank;
    bool deepks_setorb;
    bool read_in_flag;
    std::string descriptor_file;
    std::vector<std::string> orbital_file;
    double ecutwfc;
    double dk;
    double dR;
    double Rmax;

    ofs_log.open("ORB_table_phi_test.log");
    ntype = nfile;
    lmax = 3;
    out_mat_r = 0; // unused variable
    force_flag = true;
    my_rank = 0;
    deepks_setorb = false;

    read_in_flag = true;
    for (int i = 0; i != nfile; ++i) {
        orbital_file.push_back(file[i]);
    }

    // below we mimic ORB_control::read_orb_first
    ecutwfc = 10000.0; // ecutwfc has to be very large in order to reach the agreement between new & old tables
    dk = 0.01;
    dR = 0.01;
    Rmax = 30;

    lcao.read_in_flag = read_in_flag;
    lcao.orbital_file = orbital_file;
#ifdef __MPI
    lcao.bcast_files(ntype, GlobalV::MY_RANK);
#endif

    // see ORB_control::read_orb_first
    lcao.ecutwfc = ecutwfc;
    lcao.dk = dk;
    lcao.dR = dR;
    lcao.Rmax = Rmax;

    lcao.Read_Orbitals(ofs_log, ntype, lmax, deepks_setorb, out_mat_r,
            force_flag, my_rank);

    otp.allocate(ntype, lmax, lcao.get_kmesh(), Rmax, dR, dk);

    auto calc_nr = [](double const& Rmax, double const& dR) {
        int nr = static_cast<int>(Rmax / dR) + 4;
        return (nr%2) ? nr : nr+1;
    };

    // initialize spherical bessel table
    sbr.set_dx(dR*dk);

    // max(l+l) = 4, but the derivative calculation of j_l relies on j_{l+1}
    sbr.cal_jlx(2*lmax+1, calc_nr(Rmax, dR), lcao.get_kmesh());

    otp.pSB = &sbr;
    otp.init_OV_Tpair(lcao);
    otp.init_OV_Opair(lcao);

    start = iclock::now();
    otp.init_Table(lcao);
    dur = iclock::now() - start;
    std::cout << "legacy time elapsed = " << dur.count() << " s" << std::endl;

	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		for (int T2 = T1 ; T2 < ntype ; T2++)
		{
			int Tpair=otp.OV_Tpair(T1,T2);

			for (int L1 = 0; L1 <= lcao.Phi[T1].getLmax(); L1++)
			{
				for (int N1 = 0; N1 < lcao.Phi[T1].getNchi(L1); N1++)
				{
					for (int L2 = 0; L2 <= lcao.Phi[T2].getLmax(); L2 ++)
					{
						for (int N2 = 0; N2 < lcao.Phi[T2].getNchi(L2); N2++)
						{
							int Opair = otp.OV_Opair(Tpair,L1,L2,N1,N2);

							for (int L = std::abs(L1-L2); L <= (L1+L2) ; L += 2)
							{
                                // table whose integrand has a higher order of k is less accurate
                                // as it requires a larger ecutwfc to converge.

                                //if (std::abs(S_tab.table(T1, L1, N1, T2, L2, N2, L)[0] -
                                //        otp.Table_SR[0][Tpair][Opair][L][0]) > 1e-4) {
                                //    std::cout << "T1 = " << T1 << ", L1 = " << L1 << ", N1 = " << N1 << std::endl;
                                //    std::cout << "T2 = " << T2 << ", L2 = " << L2 << ", N2 = " << N2 << std::endl;
                                //    std::cout << "L = " << L << std::endl;

                                //    for (int ir = 1; ir < 10; ++ir) {
                                //        double rl = std::pow(ir * dr, L);
                                //        printf("%20.15e   %20.15e\n", S_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl, otp.Table_SR[0][Tpair][Opair][L][ir]);
                                //    }
                                //}

                                // R == 0
                                //EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L)[0],
                                //        otp.Table_SR[0][Tpair][Opair][L][0], tol);
                                //EXPECT_NEAR(S_dtab.table(T1, L1, N1, T2, L2, N2, L)[0], 
                                //        otp.Table_SR[1][Tpair][Opair][L][0], tol);
                                //EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L)[0], 
                                //        otp.Table_TR[0][Tpair][Opair][L][0], tol * 10);
                                //EXPECT_NEAR(T_dtab.table(T1, L1, N1, T2, L2, N2, L)[0], 
                                //        otp.Table_TR[1][Tpair][Opair][L][0], tol * 100);

                                // R > 0
                                for (int ir = 1; ir != 1600; ++ir) {
                                    double rl = std::pow(ir * dr, L);
                                    EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl,
                                            otp.Table_SR[0][Tpair][Opair][L][ir], tol);
                                    //EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L, true)[ir], 
                                    //        otp.Table_SR[1][Tpair][Opair][L][ir], tol);
                                    EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl, 
                                            otp.Table_TR[0][Tpair][Opair][L][ir], tol * 10);
                                    //EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L, true)[ir], 
                                    //        otp.Table_TR[1][Tpair][Opair][L][ir], tol * 100);
                                }
							}
						}
					}
				}
			}
		}
	}

    otp.Destroy_Table(lcao);
}

int main(int argc, char** argv)
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

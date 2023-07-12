#include "module_basis/module_nao/two_center_table.h"

#include "gtest/gtest.h"
#include "module_base/constants.h"
#include "module_base/spherical_bessel_transformer.h"

#include "module_base/math_integral.h"

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
    TwoCenterTable S_dtab;
    TwoCenterTable T_tab;
    TwoCenterTable T_dtab;

    TwoCenterTable betapsi_tab;
    TwoCenterTable betapsi_dtab;

    RadialCollection orb;
    int nfile = 0;                                              //! number of orbital files
    std::string* file = nullptr;                                //!< orbital files to read from
    std::string log_file = "./test_files/two_center_table.log"; //!< file for logging

    double table_tol = 1e-6; //! tolerance for comparison between new and legacy tables
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

    const double* rgrid = orb(0,0,0).ptr_rgrid();

    iclock::time_point start;
    std::chrono::duration<double> dur;

    start = iclock::now();
    S_tab.build(orb, orb, 'S', nr, rgrid, false);
    T_tab.build(orb, orb, 'T', nr, rgrid, false);
    S_dtab.build(orb, orb, 'S', nr, rgrid, true);
    T_dtab.build(orb, orb, 'T', nr, rgrid, true);
    dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " s" << std::endl;

    // check whether analytical derivative and finite difference agree
    for (int i = 0; i != S_tab.ntab(); ++i)
    {
        const double* f = S_tab.ptr_table(0, 0, 0, 0, 0, 0, 0) + i * S_tab.nr();
        const double* df = S_dtab.ptr_table(0, 0, 0, 0, 0, 0, 0) + i * S_dtab.nr();
        for (int ir = 5; ir != S_tab.nr() - 4; ++ir)
        {
            double df_fd = ( +1.0/280*f[ir-4] - 4.0/105*f[ir-3] 
                      +1.0/5*f[ir-2] - 4.0/5*f[ir-1] 
                      -1.0/280*f[ir+4] + 4.0/105*f[ir+3] 
                      -1.0/5*f[ir+2] + 4.0/5*f[ir+1]
                    ) / 0.01;

            EXPECT_NEAR(df_fd, df[ir], 1e-5);
        }
    }
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

    const double* rgrid = orb(0,0,0).ptr_rgrid();

    iclock::time_point start;
    std::chrono::duration<double> dur;

    start = iclock::now();
    S_tab.build(orb, orb, 'S', nr, rgrid, false);
    T_tab.build(orb, orb, 'T', nr, rgrid, false);
    S_dtab.build(orb, orb, 'S', nr, rgrid, true);
    T_dtab.build(orb, orb, 'T', nr, rgrid, true);
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

    //double* tmp = new double[2001];
    //const NumericalRadial& f = orb(1,3,0);
    //f.radtab('T', f, 6, tmp);
    //for (int ir = 0; ir != 21; ++ir) {
    //    std::cout << tmp[ir] << std::endl;
    //}
    //std::cout << std::endl;

//    double* ktmp = new double[f.nk()];
//    for (int ik = 0; ik != f.nk(); ++ik) {
//        ktmp[ik] = std::pow(f.ptr_kvalue()[ik] * f.ptr_kgrid()[ik] * f.ptr_kgrid()[ik],2);
//    }
//    double r1 = f.ptr_rgrid()[1];


    std::cout << std::endl;

    // compare the table
    // length of table depends on the cutoff radius
    //for (int ir = 0; ir != 1600; ++ir) {
    //    double err = std::abs(T_tab.ptr_table(1,3,0,1,3,0,6)[ir] - otp.Table_TR[0][3][80][6][ir]);
    //    if (err > table_tol) {
    //        std::cout << "                 ir = " << ir << std::endl;
    //    }
    //    EXPECT_NEAR(T_tab.ptr_table(1,3,0,1,3,0,6)[ir], otp.Table_TR[0][3][80][6][ir], table_tol);
    //    //EXPECT_NEAR(S_tab.ptr_table(0,0,0,0,0,0,0)[ir], otp.Table_SR[0][0][0][0][ir], table_tol);
    //    //EXPECT_NEAR(T_dtab.ptr_table(0,0,0,0,0,0,0)[ir], otp.Table_TR[1][0][0][0][ir], table_tol);
    //    //EXPECT_NEAR(S_dtab.ptr_table(0,0,0,0,0,0,0)[ir], otp.Table_SR[1][0][0][0][ir], table_tol);
    //}

    //std::cout << T_tab.nr() << std::endl;
    //std::cout << otp.OV_Tpair(1,1) << std::endl;
    //std::cout << otp.OV_Opair(otp.OV_Tpair(1,1), 3, 3, 0, 0) << std::endl;

    ////////////////////////////////////////////
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
                                //if (T1 == 0 && L1 == 1 && N1 == 1 && T2 == 0 && L2 == 1 && N2 == 1 && L == 0) {
                                //    for (int ir = 1; ir != 16; ++ir) {
                                //        printf("%12.8f   %12.8f\n",
                                //                T_dtab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir], 
                                //                otp.Table_TR[1][Tpair][Opair][L][ir]);
                                //    }
                                //}
                                //break;
                                double err_max = -1.0;
                                for (int ir = 1; ir != 1600; ++ir) {
                                    //double err = std::abs(T_dtab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir]
                                    //            - otp.Table_TR[1][Tpair][Opair][L][ir]);
                                    //if ( err > table_tol) {
                                    //    if (err_max < 0) {
                                    //        printf("%i  %i  %i  %i  %i  %i  %i  %5i", T1, L1, N1, T2, L2, N2, L, ir);
                                    //        err_max = err;
                                    //    } else {
                                    //        err_max = std::max(err_max, err);
                                    //    }
                                    //}

                                    // table whose integrand has a higher order of k is less accurate
                                    // as it requires a larger ecutwfc to converge.
                                    EXPECT_NEAR(T_tab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir], 
                                            otp.Table_TR[0][Tpair][Opair][L][ir], table_tol * 10);
                                    EXPECT_NEAR(S_tab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir], 
                                            otp.Table_SR[0][Tpair][Opair][L][ir], table_tol);
                                    EXPECT_NEAR(T_dtab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir], 
                                            otp.Table_TR[1][Tpair][Opair][L][ir], table_tol * 100);
                                    EXPECT_NEAR(S_dtab.ptr_table(T1, L1, N1, T2, L2, N2, L)[ir], 
                                            otp.Table_SR[1][Tpair][Opair][L][ir], table_tol);
                                }
                                if (err_max > 0) {
                                    printf("   %8.5e\n", err_max);
                                }
							}
						}
					}
				}
			}
		}
	}


    //////////////////////////////////////////////
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

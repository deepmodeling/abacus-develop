#include "module_basis/module_nao/two_center_integrator.h"

#include "gtest/gtest.h"
#include "module_base/constants.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_base/vector3.h"
#include "module_base/ylm.h"

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
 *      Unit test of class "TwoCenterIntegrator"
 ***********************************************************/
/*!
 *  Tested functions:
 *
 *  - build
 *      - builds an object for doing a specific two-center integral
 *                                                                      */
class TwoCenterIntegratorTest : public ::testing::Test
{
  protected:
    void SetUp();
    void TearDown();

    TwoCenterIntegrator S_intor;
    TwoCenterIntegrator T_intor;

    RadialCollection orb;
    int nfile = 0;                                                   //! number of orbital files
    std::string* file = nullptr;                                     //!< orbital files to read from
    std::string log_file = "./test_files/two_center_integrator.log"; //!< file for logging

    double elem_tol = 1e-6; //! tolerance for comparison between new and legacy matrix elements
};

void TwoCenterIntegratorTest::SetUp()
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

    ModuleBase::Ylm::set_coefficients();
}

void TwoCenterIntegratorTest::TearDown()
{
    delete[] file;
}

TEST_F(TwoCenterIntegratorTest, FiniteDifference)
{
    nfile = 3;
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

    S_intor.tabulate(orb, orb, 'S', nr, rmax);
    T_intor.tabulate(orb, orb, 'T', nr, rmax);

    dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " s" << std::endl;

    // check whether analytical derivative and finite difference agree
    int ntype = nfile;
    double tol_d = 1e-4;

    ModuleBase::Vector3<double> vR0 = {1.0, 2.0, 3.0};
    ModuleBase::Vector3<double> vR;

	for (int t1 = 0;  t1 < ntype ; t1++)
	{
		for (int l1 = 0; l1 <= orb(t1).lmax(); l1++)
		{
			for (int izeta1 = 0; izeta1 < orb(t1).nzeta(l1); izeta1++)
			{
                for (int m1 = -l1; m1 <= l1; ++m1)
                {
		            for (int t2 = t1 ; t2 < ntype ; t2++)
		            {
					    for (int l2 = 0; l2 <= orb(t2).lmax(); l2++)
					    {
					    	for (int izeta2 = 0; izeta2 < orb(t2).nzeta(l2); izeta2++)
					    	{
                                for (int m2 = -l2; m2 <= l2; ++m2)
                                {
                                    double dx = 1e-4;
                                    double elem_p;
                                    double elem_m;
                                    double grad_elem[3];

                                    // S
                                    vR = vR0;
                                    vR[2] += dx;
                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, false, &elem_p);

                                    vR = vR0;
                                    vR[2] -= dx;
                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, false, &elem_m);

                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, true, grad_elem);

                                    EXPECT_NEAR( (elem_p - elem_m) / (2. * dx), grad_elem[2], tol_d);

                                    // T
                                    vR = vR0;
                                    vR[2] += dx;
                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, false, &elem_p);

                                    vR = vR0;
                                    vR[2] -= dx;
                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, false, &elem_m);

                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, true, grad_elem);

                                    EXPECT_NEAR( (elem_p - elem_m) / (2. * dx), grad_elem[2], tol_d);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//TEST_F(TwoCenterIntegratorTest, LegacyConsistency)
//{
//    // use less files so that the test wouldn't take too long
//    nfile = 3;
//
//    orb.build(nfile, file, 'o');
//    double rmax = orb.rcut_max() * 2.0;
//    double dr = 0.01;
//    int nr = static_cast<int>(rmax / dr) + 1;
//    
//    //ModuleBase::SphericalBesselTransformer sbt;
//    //orb.set_transformer(&sbt, 0);
//    orb.set_uniform_grid(true, nr, rmax, 'i', true);
//
//    iclock::time_point start;
//    std::chrono::duration<double> dur;
//
//    start = iclock::now();
//    S_tab.build(orb, orb, 'S', nr, rmax, true);
//    T_tab.build(orb, orb, 'T', nr, rmax, true);
//    dur = iclock::now() - start;
//    std::cout << "radfft time elapsed = " << dur.count() << " s" << std::endl;
//
//    //////////////////////////////////////////////
//    //      module_ao ORB_table_phi
//    //////////////////////////////////////////////
//
//    ORB_table_phi otp;
//    LCAO_Orbitals lcao;
//    ModuleBase::Sph_Bessel_Recursive::D2 sbr;
//
//    // prepare data required to initialized ORB_table_phi
//    std::ofstream ofs_log;
//    int ntype;
//    int lmax;
//    int out_mat_r; // unused variable
//    bool force_flag;
//    int my_rank;
//    bool deepks_setorb;
//    bool read_in_flag;
//    std::string descriptor_file;
//    std::vector<std::string> orbital_file;
//    double ecutwfc;
//    double dk;
//    double dR;
//    double Rmax;
//
//    ofs_log.open("ORB_table_phi_test.log");
//    ntype = nfile;
//    lmax = 3;
//    out_mat_r = 0; // unused variable
//    force_flag = true;
//    my_rank = 0;
//    deepks_setorb = false;
//
//    read_in_flag = true;
//    for (int i = 0; i != nfile; ++i) {
//        orbital_file.push_back(file[i]);
//    }
//
//    // below we mimic ORB_control::read_orb_first
//    ecutwfc = 10000.0; // ecutwfc has to be very large in order to reach the agreement between new & old tables
//    dk = 0.01;
//    dR = 0.01;
//    Rmax = 30;
//
//    lcao.read_in_flag = read_in_flag;
//    lcao.orbital_file = orbital_file;
//#ifdef __MPI
//    lcao.bcast_files(ntype, GlobalV::MY_RANK);
//#endif
//
//    // see ORB_control::read_orb_first
//    lcao.ecutwfc = ecutwfc;
//    lcao.dk = dk;
//    lcao.dR = dR;
//    lcao.Rmax = Rmax;
//
//    lcao.Read_Orbitals(ofs_log, ntype, lmax, deepks_setorb, out_mat_r,
//            force_flag, my_rank);
//
//    otp.allocate(ntype, lmax, lcao.get_kmesh(), Rmax, dR, dk);
//
//    auto calc_nr = [](double const& Rmax, double const& dR) {
//        int nr = static_cast<int>(Rmax / dR) + 4;
//        return (nr%2) ? nr : nr+1;
//    };
//
//    // initialize spherical bessel table
//    sbr.set_dx(dR*dk);
//
//    // max(l+l) = 4, but the derivative calculation of j_l relies on j_{l+1}
//    sbr.cal_jlx(2*lmax+1, calc_nr(Rmax, dR), lcao.get_kmesh());
//
//    otp.pSB = &sbr;
//    otp.init_OV_Tpair(lcao);
//    otp.init_OV_Opair(lcao);
//
//    start = iclock::now();
//    otp.init_Table(lcao);
//    dur = iclock::now() - start;
//    std::cout << "legacy time elapsed = " << dur.count() << " s" << std::endl;
//
//	for (int T1 = 0;  T1 < ntype ; T1++)
//	{
//		for (int T2 = T1 ; T2 < ntype ; T2++)
//		{
//			int Tpair=otp.OV_Tpair(T1,T2);
//
//			for (int L1 = 0; L1 <= lcao.Phi[T1].getLmax(); L1++)
//			{
//				for (int N1 = 0; N1 < lcao.Phi[T1].getNchi(L1); N1++)
//				{
//					for (int L2 = 0; L2 <= lcao.Phi[T2].getLmax(); L2 ++)
//					{
//						for (int N2 = 0; N2 < lcao.Phi[T2].getNchi(L2); N2++)
//						{
//							int Opair = otp.OV_Opair(Tpair,L1,L2,N1,N2);
//
//							for (int L = std::abs(L1-L2); L <= (L1+L2) ; L += 2)
//							{
//                                // table whose integrand has a higher order of k is less accurate
//                                // as it requires a larger ecutwfc to converge.
//
//                                //if (std::abs(S_tab.table(T1, L1, N1, T2, L2, N2, L)[0] -
//                                //        otp.Table_SR[0][Tpair][Opair][L][0]) > 1e-4) {
//                                //    std::cout << "T1 = " << T1 << ", L1 = " << L1 << ", N1 = " << N1 << std::endl;
//                                //    std::cout << "T2 = " << T2 << ", L2 = " << L2 << ", N2 = " << N2 << std::endl;
//                                //    std::cout << "L = " << L << std::endl;
//
//                                //    for (int ir = 1; ir < 10; ++ir) {
//                                //        double rl = std::pow(ir * dr, L);
//                                //        printf("%20.15e   %20.15e\n", S_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl, otp.Table_SR[0][Tpair][Opair][L][ir]);
//                                //    }
//                                //}
//
//                                // R == 0
//                                //EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L)[0],
//                                //        otp.Table_SR[0][Tpair][Opair][L][0], table_tol);
//                                //EXPECT_NEAR(S_dtab.table(T1, L1, N1, T2, L2, N2, L)[0], 
//                                //        otp.Table_SR[1][Tpair][Opair][L][0], table_tol);
//                                //EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L)[0], 
//                                //        otp.Table_TR[0][Tpair][Opair][L][0], table_tol * 10);
//                                //EXPECT_NEAR(T_dtab.table(T1, L1, N1, T2, L2, N2, L)[0], 
//                                //        otp.Table_TR[1][Tpair][Opair][L][0], table_tol * 100);
//
//                                // R > 0
//                                for (int ir = 1; ir != 1600; ++ir) {
//                                    double rl = std::pow(ir * dr, L);
//                                    EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl,
//                                            otp.Table_SR[0][Tpair][Opair][L][ir], table_tol);
//                                    EXPECT_NEAR(S_tab.table(T1, L1, N1, T2, L2, N2, L, true)[ir], 
//                                            otp.Table_SR[1][Tpair][Opair][L][ir], table_tol);
//                                    EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L)[ir] * rl, 
//                                            otp.Table_TR[0][Tpair][Opair][L][ir], table_tol * 10);
//                                    EXPECT_NEAR(T_tab.table(T1, L1, N1, T2, L2, N2, L, true)[ir], 
//                                            otp.Table_TR[1][Tpair][Opair][L][ir], table_tol * 100);
//                                }
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//    otp.Destroy_Table(lcao);
//}

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

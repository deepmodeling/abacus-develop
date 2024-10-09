#include "module_basis/module_nao/two_center_table.h"

#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/spherical_bessel_transformer.h"

#include "gtest/gtest.h"
#include <chrono>
#include <cstring>
using iclock = std::chrono::high_resolution_clock;

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

    ModuleBase::SphericalBesselTransformer sbt;
    orb.set_transformer(sbt);

    double rmax = orb.rcut_max() * 2.0;
    double dr = 0.01;
    int nr = static_cast<int>(rmax / dr) + 1;

    // ModuleBase::SphericalBesselTransformer sbt;
    // sbt.set_fftw_plan_flag(FFTW_MEASURE); // not necessarily worth it!
    // orb.set_transformer(&sbt, 0);
    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    iclock::time_point start;
    std::chrono::duration<double> dur;

    start = iclock::now();
    S_tab.build(orb, orb, 'S', nr, rmax);
    T_tab.build(orb, orb, 'T', nr, rmax);
    dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " s" << std::endl;

    // check whether the derivative table agrees with the finite difference of table
    int ntype = nfile;
    for (int T1 = 0; T1 < ntype; T1++)
    {
        for (int T2 = T1; T2 < ntype; T2++)
        {
            for (int L1 = 0; L1 <= orb(T1).lmax(); L1++)
            {
                for (int N1 = 0; N1 < orb(T1).nzeta(L1); N1++)
                {
                    for (int L2 = 0; L2 <= orb(T2).lmax(); L2++)
                    {
                        for (int N2 = 0; N2 < orb(T2).nzeta(L2); N2++)
                        {
                            for (int L = std::abs(L1 - L2); L <= (L1 + L2); L += 2)
                            {
                                const double* f = S_tab.table(T1, L1, N1, T2, L2, N2, L, false);
                                const double* df = S_tab.table(T1, L1, N1, T2, L2, N2, L, true);

                                for (int ir = 4; ir != S_tab.nr() - 4; ++ir)
                                {
                                    double df_fd
                                        = (-1.0 / 280 * (f[ir + 4] - f[ir - 4]) + 4.0 / 105 * (f[ir + 3] - f[ir - 3])
                                           - 0.2 * (f[ir + 2] - f[ir - 2]) + 0.8 * (f[ir + 1] - f[ir - 1]))
                                          / dr;

                                    // df is (d/dR)(S/R^l), it should be rescaled to have the
                                    // same unit as dS/dR in order to have meaningful error comparison
                                    double err_abs = std::abs(df_fd - df[ir]) * std::pow(ir * dr, L);
                                    double err_rel = std::abs((df_fd - df[ir]) / df[ir]);
                                    if (err_abs > tol_d_abs && err_rel > tol_d_rel)
                                    {
                                        printf("T1 = %i   L1 = %i   N1 = %i   T2 = %i   L2 = %i   N2 = %i   L = %i   "
                                               "ir = %2i   df_fd = % 8.5e   df_tab = % 8.5e   err_abs = %8.5e   "
                                               "err_rel = %8.5e\n",
                                               T1,
                                               L1,
                                               N1,
                                               T2,
                                               L2,
                                               N2,
                                               L,
                                               ir,
                                               df_fd,
                                               df[ir],
                                               err_abs,
                                               err_rel);
                                    }

                                    EXPECT_TRUE(err_abs < tol_d_abs || err_rel < tol_d_rel);
                                }

                                f = T_tab.table(T1, L1, N1, T2, L2, N2, L, false);
                                df = T_tab.table(T1, L1, N1, T2, L2, N2, L, true);

                                for (int ir = 4; ir != T_tab.nr() - 4; ++ir)
                                {
                                    double df_fd
                                        = (-1.0 / 280 * (f[ir + 4] - f[ir - 4]) + 4.0 / 105 * (f[ir + 3] - f[ir - 3])
                                           - 0.2 * (f[ir + 2] - f[ir - 2]) + 0.8 * (f[ir + 1] - f[ir - 1]))
                                          / dr;

                                    // df is (d/dR)(S/R^l), it should be rescaled to have the
                                    // same unit as dS/dR in order to have meaningful error comparison
                                    double err_abs = std::abs(df_fd - df[ir]) * std::pow(ir * dr, L);
                                    double err_rel = std::abs((df_fd - df[ir]) / df[ir]);
                                    if (err_abs > tol_d_abs && err_rel > tol_d_rel)
                                    {
                                        printf("T1 = %i   L1 = %i   N1 = %i   T2 = %i   L2 = %i   N2 = %i   L = %i   "
                                               "ir = %2i   df_fd = % 8.5e   df_tab = % 8.5e   err_abs = %8.5e   "
                                               "err_rel = %8.5e\n",
                                               T1,
                                               L1,
                                               N1,
                                               T2,
                                               L2,
                                               N2,
                                               L,
                                               ir,
                                               df_fd,
                                               df[ir],
                                               err_abs,
                                               err_rel);
                                    }

                                    EXPECT_TRUE(err_abs < tol_d_abs || err_rel < tol_d_rel);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    EXPECT_EQ(S_tab.nr(), nr);
    EXPECT_EQ(T_tab.nr(), nr);

    EXPECT_EQ(S_tab.rmax(), rmax);
    EXPECT_EQ(T_tab.rmax(), rmax);
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

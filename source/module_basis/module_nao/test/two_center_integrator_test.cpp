#include "module_basis/module_nao/two_center_integrator.h"

#include "module_base/constants.h"
#include "module_base/math_sphbes.h"
#include "module_base/spherical_bessel_transformer.h"
#include "module_base/vector3.h"
#include "module_base/ylm.h"

#include "gtest/gtest.h"
#include <chrono>
#include <cstring>
using iclock = std::chrono::high_resolution_clock;
using ModuleBase::Sphbes;

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

    S_intor.tabulate(orb, orb, 'S', nr, rmax);
    T_intor.tabulate(orb, orb, 'T', nr, rmax);

    dur = iclock::now() - start;
    std::cout << "time elapsed = " << dur.count() << " s" << std::endl;

    // check whether analytical derivative and finite difference agree
    int ntype = nfile;
    double tol_d = 1e-4;

    ModuleBase::Vector3<double> vR0 = {1.0, 2.0, 3.0};
    ModuleBase::Vector3<double> vR;

    for (int t1 = 0; t1 < ntype; t1++)
    {
        for (int l1 = 0; l1 <= orb(t1).lmax(); l1++)
        {
            for (int izeta1 = 0; izeta1 < orb(t1).nzeta(l1); izeta1++)
            {
                for (int m1 = -l1; m1 <= l1; ++m1)
                {
                    for (int t2 = t1; t2 < ntype; t2++)
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
                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, &elem_p);

                                    vR = vR0;
                                    vR[2] -= dx;
                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, &elem_m);

                                    S_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, nullptr, grad_elem);

                                    EXPECT_NEAR((elem_p - elem_m) / (2. * dx), grad_elem[2], tol_d);

                                    // T
                                    vR = vR0;
                                    vR[2] += dx;
                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, &elem_p);

                                    vR = vR0;
                                    vR[2] -= dx;
                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, &elem_m);

                                    T_intor.calculate(t1, l1, izeta1, m1, t2, l2, izeta2, m2, vR, nullptr, grad_elem);

                                    EXPECT_NEAR((elem_p - elem_m) / (2. * dx), grad_elem[2], tol_d);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_F(TwoCenterIntegratorTest, SphericalBessel)
{
    int lmax = 3;
    int nbes = 5;
    int rcut = 7.0;
    double sigma = 0.0;
    double dr = 0.005;
    // The truncated spherical Bessel function has discontinuous first and
    // second derivative at the cutoff, so a small "dr" is required in order
    // to achieve sufficient accuracy.
    //
    // for dr = 0.01, the error of kinetic matrix element is about 1.5e-3
    // for dr = 0.001, the error of kinetic matrix element is about 1.5e-4

    orb.build(lmax, nbes, rcut, sigma, dr);

    ModuleBase::SphericalBesselTransformer sbt;
    orb.set_transformer(sbt);

    double rmax = orb.rcut_max() * 2.0;
    int nr = static_cast<int>(rmax / dr) + 1;

    orb.set_uniform_grid(true, nr, rmax, 'i', true);

    S_intor.tabulate(orb, orb, 'S', nr, rmax);
    T_intor.tabulate(orb, orb, 'T', nr, rmax);

    ModuleBase::Vector3<double> R0 = {0.0, 0.0, 0.0};

    // zeros of spherical bessel functions
    double* zeros = new double[nbes * (lmax + 1)];
    Sphbes::sphbes_zeros(lmax, nbes, zeros, true);

    // checks the diagonal elements with analytical expression
    double elem, ref;
    for (int l = 0; l <= lmax; ++l)
    {
        for (int zeta = 0; zeta < nbes; ++zeta)
        {
            S_intor.calculate(0, l, zeta, 0, 0, l, zeta, 0, R0, &elem);
            ref = 0.5 * std::pow(rcut, 3) * std::pow(Sphbes::sphbesj(l + 1, zeros[l * nbes + zeta]), 2);
            EXPECT_NEAR(elem, ref, 1e-5);

            T_intor.calculate(0, l, zeta, 0, 0, l, zeta, 0, R0, &elem);
            ref = 0.5 * rcut * std::pow(zeros[l * nbes + zeta] * Sphbes::sphbesj(l + 1, zeros[l * nbes + zeta]), 2);
            EXPECT_NEAR(elem, ref, 1e-3);

            // orthogonality
            for (int zeta2 = 0; zeta2 < zeta; ++zeta2)
            {
                S_intor.calculate(0, l, zeta, 0, 0, l, zeta2, 0, R0, &elem);
                ref = 0.0;
                EXPECT_NEAR(elem, ref, 1e-5);
            }
        }
    }
    delete[] zeros;
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

#include "../math_sphbes.h"

#include <cmath>
#include <fstream>
#include <iostream>

#ifdef __MPI
#include "mpi.h"
#endif

#include "gtest/gtest.h"

#define doublethreshold 1e-7
#define MAX(x, y) ((x) > (y) ? (x) : (y))

/************************************************
 *  unit test of class Integral
 ***********************************************/

/**
 * Note: this unit test try to ensure the invariance
 * of the spherical Bessel produced by class Sphbes,
 * and the reference results are produced by ABACUS
 * at 2022-1-27.
 *
 * Tested function:
 *      - Spherical_Bessel.
 *      - Spherical_Bessel_Roots
 *      - overloading of Spherical_Bessel. This funnction sets sjp[i] to 1.0 when i < msh.
 *      - sphbesj
 *      - sphbes_zeros
 */

double mean(const double* vect, const int totN)
{
    double meanv = 0.0;
    for (int i = 0; i < totN; ++i)
    {
        meanv += vect[i] / totN;
    }
    return meanv;
}

class Sphbes : public testing::Test
{
  protected:
    int msh = 700;
    int l0 = 0;
    int l1 = 1;
    int l2 = 2;
    int l3 = 3;
    int l4 = 4;
    int l5 = 5;
    int l6 = 6;
    int l7 = 7;
    double q = 1.0;
    double* r = new double[msh];
    double* jl = new double[msh];
    double* djl = new double[msh];

    void SetUp()
    {
        for (int i = 0; i < msh; ++i)
        {
            r[i] = 0.01 * (i);
        }
    }

    void TearDown()
    {
        delete[] r;
        delete[] jl;
        delete[] djl;
    }
};

TEST_F(Sphbes, Constructor)
{
    EXPECT_NO_THROW(ModuleBase::Sphbes sb);
}

TEST_F(Sphbes, SphericalBessel)
{
    // int l = 0;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l0, jl);
    // reference result is from bessel_test.cpp which is calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1
    EXPECT_NEAR(mean(jl, msh) / 0.2084468748396, 1.0, doublethreshold);

    // int l = 1;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l1, jl);
    // reference result is from bessel_test.cpp which is calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1
    EXPECT_NEAR(mean(jl, msh) / 0.12951635180384, 1.0, doublethreshold);

    // int l = 2;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l2, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.124201140093879
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.12420114009271901456, 1.0, doublethreshold);

    // int l = 3;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l3, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.118268654505568
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.11826865448477598408, 1.0, doublethreshold);

    // int l = 4;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l4, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.0933871035384385
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.093387100084701621383, 1.0, doublethreshold);

    // int l = 5;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l5, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.0603800487910689
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.060380048719821471925, 1.0, doublethreshold);

    // int l = 6;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l6, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.0327117051555907
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.032711705053977857549, 1.0, doublethreshold);

    // int l = 7;
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l7, jl);
    // the result from bessel_test.cpp calculated by
    // ModuleBase::Sph_Bessel_Recursive::D1 is 0.0152155566653926
    // reference result is calculated by Sphbes::Spherical_Bessel(msh,r,q,l,jl)
    EXPECT_NEAR(mean(jl, msh) / 0.015215556095798710851, 1.0, doublethreshold);
}

TEST_F(Sphbes, dSpherical_Bessel_dx)
{
    double djl0;
    for (int il = 0; il <= l7; ++il)
    {
        if (il == 1)
            djl0 = 1.0 / 3.0;
        else
            djl0 = 0.0;
        ModuleBase::Sphbes::dSpherical_Bessel_dx(msh, r, q, il, djl);
        ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, il, jl);
        EXPECT_NEAR(djl[0], djl0, 1e-8);
        for (int i = 1; i < msh - 1; ++i)
        {
            if (jl[i - 1] < 1e-8)
                continue;
            double djl_diff = (jl[i + 1] - jl[i - 1]) / (q * (r[i + 1] - r[i - 1]));
            EXPECT_NEAR(djl[i], djl_diff, 1e-4);
        }
        ModuleBase::Sphbes::dSpherical_Bessel_dx(msh, r, 0, il, djl);
        for (int i = 0; i < msh; ++i)
        {
            EXPECT_NEAR(djl[i], djl0, 1e-8);
        }
    }
}

TEST_F(Sphbes, SphericalBesselRoots)
{
    int neign = 100;
    double** eign = new double*[8];
    for (int i = 0; i < 8; ++i)
    {
        eign[i] = new double[neign];
        ModuleBase::Sphbes::Spherical_Bessel_Roots(neign, i, 1.0e-12, eign[i], 10.0);
    }

    EXPECT_NEAR(eign[0][0] / 0.31415926535899563188, 1.0, doublethreshold);
    EXPECT_NEAR(eign[0][99] / 31.415926535896932847, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[0], 100) / 15.865042900628463229, 1.0, doublethreshold);
    EXPECT_NEAR(eign[1][0] / 0.44934094579091843347, 1.0, doublethreshold);
    EXPECT_NEAR(eign[1][99] / 31.572689440204385392, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[1], 100) / 16.020655759558295017, 1.0, doublethreshold);
    EXPECT_NEAR(eign[2][0] / 0.57634591968946913276, 1.0, doublethreshold);
    EXPECT_NEAR(eign[2][99] / 31.729140298172534784, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[2], 100) / 16.175128483074864505, 1.0, doublethreshold);
    EXPECT_NEAR(eign[3][0] / 0.69879320005004752492, 1.0, doublethreshold);
    EXPECT_NEAR(eign[3][99] / 31.885283678838447941, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[3], 100) / 16.328616567969248763, 1.0, doublethreshold);
    EXPECT_NEAR(eign[4][0] / 0.81825614525711076741, 1.0, doublethreshold);
    EXPECT_NEAR(eign[4][99] / 32.041124042016576823, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[4], 100) / 16.481221742387987206, 1.0, doublethreshold);
    EXPECT_NEAR(eign[5][0] / 0.93558121110426506473, 1.0, doublethreshold);
    EXPECT_NEAR(eign[5][99] / 32.196665741899131774, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[5], 100) / 16.633019118735202113, 1.0, doublethreshold);
    EXPECT_NEAR(eign[6][0] / 1.051283540809391015, 1.0, doublethreshold);
    EXPECT_NEAR(eign[6][99] / 32.351913030537232885, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[6], 100) / 16.784067905062840964, 1.0, doublethreshold);
    EXPECT_NEAR(eign[7][0] / 1.1657032192516516567, 1.0, doublethreshold);
    EXPECT_NEAR(eign[7][99] / 32.506870061157627561, 1.0, doublethreshold);
    EXPECT_NEAR(mean(eign[7], 100) / 16.934416735327332049, 1.0, doublethreshold);

    for (int i = 0; i < 8; ++i)
        delete[] eign[i];
    delete[] eign;
}

// TEST_F(Sphbes, SeriesAndRecurrence)
// {
//     // This test checks whether Spherical_Bessel and sphbesj agree with each other
//     // on a coarse grid for a range of l and q values.
//     //
//     // NOTE: this test should be removed once Spherical_Bessel is removed from the code.
//     int lmax = 8;
//     int nr = 5000;
//     double rcut = 50;
//     double dr = rcut / (nr - 1);
//     double* r = new double[nr];
//     for (int i = 0; i < nr; ++i)
//     {
//         r[i] = i * dr;
//     }

//     double* jl_old = new double[nr];
//     double* jl_new = new double[nr];
//     for (int l = 0; l <= lmax; ++l)
//     {
//         for (double q = 0.1; q < 10; q += 0.1)
//         {
//             ModuleBase::Sphbes::Spherical_Bessel(nr, r, q, l, jl_old);
//             ModuleBase::Sphbes::sphbesj(nr, r, q, l, jl_new);
//             for (int i = 0; i < nr; ++i)
//             {
//                 EXPECT_NEAR(jl_old[i], jl_new[i], 1e-12);
//             }

//             // derivative
//             ModuleBase::Sphbes::dSpherical_Bessel_dx(nr, r, q, l, jl_old);
//             ModuleBase::Sphbes::dsphbesj(nr, r, q, l, jl_new);
//             for (int i = 0; i < nr; ++i)
//             {
//                 EXPECT_NEAR(jl_old[i], jl_new[i], 1e-12);
//             }
//         }
//     }

//     delete[] r;
//     delete[] jl_old;
//     delete[] jl_new;
// }

TEST_F(Sphbes, SphericalBesselPrecisionGrid)
{
    // This test checks whether sphbesj agrees with the Octave implementation
    // on a coarse grid for a range of l and q values.
    const double q = 0.1;
    const double rcut = 50;
    const int nr = 5000;
    const int l_lo = 0, l_hi = 11;
    const double dr = rcut / nr;
    double* r = new double[nr + 10];
    double* Y = new double[nr * (l_hi - l_lo + 1) + 10];
    // save errs to binary file
    std::ofstream file_o("data/sphj_old_out.bin", std::ios::binary);
    std::ofstream file_n("data/sphj_new_out.bin", std::ios::binary);

    // case 0: x = 0, l = 0
    EXPECT_NEAR(ModuleBase::Sphbes::sphbesj(0, 0), 1.0, 1e-12);
    // case 1: x = 0, l = 1, ... , 11
    for (int l = 1; l <= l_hi; ++l)
    {
        EXPECT_NEAR(ModuleBase::Sphbes::sphbesj(l, 0), 0.0, 1e-12);
    }
    // case 2: wide range of x and l
    double y;
    // read reference data
    std::ifstream fin("data/bjo.bin", std::ios::binary);
    int i = 0;
    while (fin.read(reinterpret_cast<char*>(&y), sizeof(double)))
    {
        Y[i] = y;
        ++i;
    }
    fin.close();

    for (int i = 0; i < nr; ++i)
    {
        r[i] = (i + 1) * dr;
    }

    // test for new sphbesj
    for (int l = l_lo; l <= l_hi; ++l)
    {
        for (int i = 0; i < nr; ++i)
        {
            EXPECT_NEAR(ModuleBase::Sphbes::sphbesj(l, r[i] * q), Y[l * nr + i], 1e-12);
            double tmp = std::abs(Y[l * nr + i] - ModuleBase::Sphbes::sphbesj(l, r[i] * q));
            file_n.write(reinterpret_cast<char*>(&tmp), sizeof(double));
        }
    }
    // test for old Bessel
    // most of l cases precision failed to achieve 1e-12
    double* jl_old = new double[nr + 10];
    for (int l = l_lo; l <= l_hi; ++l)
    {
        ModuleBase::Sphbes::Spherical_Bessel(nr, r, q, l, jl_old);
        for (int i = 0; i < nr; ++i)
        {
            double tmp = std::abs(jl_old[i] - Y[l * nr + i]);
            file_o.write(reinterpret_cast<char*>(&tmp), sizeof(double));
        }
    }

    delete[] r;
    delete[] Y;
    delete[] jl_old;
    file_o.close();
    file_n.close();
}

TEST_F(Sphbes, SphericalBesselPrecisionNearZero)
{
    // This test checks whether sphbesj agrees with the Octave implementation
    // when x is near zero point for a range of l.
    const int n = 16;
    const int l_lo = 0, l_hi = 11;
    double* x = new double[n + 10];
    double* Y = new double[n * (l_hi - l_lo + 1) + 10];

    // read reference data
    double y;
    std::ifstream fin("data/bjxo.bin", std::ios::binary);
    int i = 0;
    while (fin.read(reinterpret_cast<char*>(&y), sizeof(double)))
    {
        Y[i] = y;
        ++i;
    }
    fin.close();
    // generate x
    x[0] = 1.0 / (1 << 5);
    for (int i = 1; i < n; i++)
    {
        x[i] = x[i - 1] / 2;
    }

    // test for sphbesj near zero
    for (int l = l_lo; l <= l_hi; ++l)
    {
        for (int i = 0; i < n; ++i)
        {
            EXPECT_NEAR(ModuleBase::Sphbes::sphbesj(l, x[i]), Y[l * n + i], 1e-12);
        }
    }
    delete[] x;
    delete[] Y;
}

TEST_F(Sphbes, Zeros)
{
    // This test checks whether sphbes_zeros properly computes the zeros of sphbesj.

    int lmax = 20;
    int nzeros = 500;
    double* zeros = new double[nzeros*(lmax+1)];
    for (int l = 0; l <= lmax; ++l)
    {
        ModuleBase::Sphbes::sphbes_zeros(l, nzeros, zeros, false);
        for (int i = 0; i < nzeros; ++i)
        {
            EXPECT_LT(std::abs(ModuleBase::Sphbes::sphbesj(l, zeros[i])), 1e-14);
        }
    }


    ModuleBase::Sphbes::sphbes_zeros(lmax, nzeros, zeros, true);
    for (int l = 0; l <= lmax; ++l)
    {
        for (int i = 0; i < nzeros; ++i)
        {
            EXPECT_LT(std::abs(ModuleBase::Sphbes::sphbesj(l, zeros[l*nzeros+i])), 1e-14);
        }
    }

    delete[] zeros;
}

TEST_F(Sphbes, ZerosOld)
{
    // This test checks whether Spherical_Bessel_Roots properly computes the zeros of sphbesj.

    int lmax = 7;
    int nzeros = 50;
    double* zeros = new double[nzeros];
    for (int l = 0; l <= lmax; ++l)
    {
        ModuleBase::Sphbes::Spherical_Bessel_Roots(nzeros, l, 1e-7, zeros, 1.0);
        for (int i = 0; i < nzeros; ++i)
        {
            EXPECT_LT(std::abs(ModuleBase::Sphbes::sphbesj(l, zeros[i])), 1e-7);
        }
    }
}

TEST_F(Sphbes, Derivatives)
{
    int lmax = 20;
    int numr = 20;
    double* r = new double[numr];
    double* djl = new double[numr];
    double q = 0.001;
    r[0] = 1.0;
    for (int i = 0; i < numr; ++i)
    {
        r[i + 1] = r[i] * 2.0;
    }

    for (int l = 0; l <= lmax; ++l)
    {
        ModuleBase::Sphbes::dsphbesj(numr, r, q, l, djl);
        for (int i = 0; i < numr; ++i)
        {
            double h = 1e-8;
            EXPECT_LT(
                abs(djl[i] * 2 * h
                    - (ModuleBase::Sphbes::sphbesj(l, q * r[i] + h) - ModuleBase::Sphbes::sphbesj(l, q * r[i] - h))),
                1e-14);
        }
    }
}

TEST_F(Sphbes, DerivativesOld)
{
    int lmax = 20;
    int numr = 20;
    double* r = new double[numr];
    double* djl = new double[numr];
    double q = 0.001;
    r[0] = 1.0;
    for (int i = 0; i < numr; ++i)
    {
        r[i + 1] = r[i] * 2.0;
    }

    for (int l = 0; l < lmax; l++)
    {
        ModuleBase::Sphbes::dSpherical_Bessel_dx(numr, r, q, l, djl);
        for (int i = 0; i < numr; i++)
        {
            double h = 1e-8;
            double errs
                = abs(djl[i] * 2 * h
                      - (ModuleBase::Sphbes::sphbesj(l, q * r[i] + h) - ModuleBase::Sphbes::sphbesj(l, q * r[i] - h)));
            if (errs > 1e-14)
            {
                std::cout << "l = " << l << ", r = " << r[i] << ", errs = " << errs << std::endl;
            }
        }
    }
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

TEST_F(Sphbes, SphericalBesselsjp)
{
    int iii;
    double* sjp = new double[msh];
    ModuleBase::Sphbes::Spherical_Bessel(msh, r, q, l0, jl, sjp);
    EXPECT_NEAR(mean(jl, msh) / 0.2084468748396, 1.0, doublethreshold);
    for (int iii = 0; iii < msh; ++iii)
    {
        EXPECT_EQ(sjp[iii], 1.0);
    }
    delete[] sjp;
}

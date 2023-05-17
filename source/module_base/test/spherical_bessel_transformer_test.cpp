#include "module_base/spherical_bessel_transformer.h"

#include "gtest/gtest.h"
#include "module_base/constants.h"

#ifdef __MPI
#include <mpi.h>
#endif

using ModuleBase::PI;
using ModuleBase::SphericalBesselTransformer;

/***********************************************************
 *      Unit test of class SphericalBesselTransform
 ***********************************************************/

/*! Tested functions:
 *
 *  - radrfft
 *      - Performs a spherical Bessel transform via fast
 *        Fourier transforms.
 *
 *  - set_fftw_plan_flag
 *      - Sets the planner flag for FFTW plan creation.
 *
 *                                                          */

class SphericalBesselTransformTest : public ::testing::Test
{

  protected:
    /// Allocates buffers
    void SetUp();

    /// Deallocates buffers
    void TearDown();

    /// Gets the maximum absolute element-wise difference between two arrays
    double max_diff(const int sz, const double* const arr1, const double* const arr2);

    int sz_max = 10000;      ///< size of each buffer
    double* f = nullptr;     ///< buffer for input array
    double* g = nullptr;     ///< buffer for output array
    double* g_ref = nullptr; ///< buffer for reference array

    double tol = 1e-9; ///< tolerance for element-wise numerical error
};

void SphericalBesselTransformTest::SetUp()
{
    f = new double[sz_max];
    g = new double[sz_max];
    g_ref = new double[sz_max];
}

void SphericalBesselTransformTest::TearDown()
{
    delete[] f;
    delete[] g;
    delete[] g_ref;
}

double SphericalBesselTransformTest::max_diff(const int sz, const double* const arr1, const double* const arr2)
{
    double diff = 0.0;
    double tmp = 0.0;
    for (int i = 0; i < sz; ++i)
    {
        tmp = std::abs(arr1[i] - arr2[i]);
        if (tmp > diff)
        {
            diff = tmp;
        }
    }
    return diff;
}

TEST_F(SphericalBesselTransformTest, BasicFunctionality)
{
    /*
     * Computes the zeroth, first and second order spherical Bessel
     * transforms of r*exp(-r) and compares the results with analytic
     * expressions:
     *
     * zeroth:  2*sqrt(2/pi) * (3-k^2) / (k^2+1)^3.
     * first :  2*sqrt(2/pi) *   4k    / (k^2+1)^3.
     * second:  2*sqrt(2/pi) *   4k^2  / (k^2+1)^3.
     *                                                              */
    int sz = 10000;
    assert(sz <= sz_max);

    double dr = 0.01;
    double rcut = dr * (sz - 1);
    double dk = PI / rcut;
    double pref = std::sqrt(2. / PI) * 2.;

    SphericalBesselTransformer sbt;

    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = r * std::exp(-r);
    }

    // zeroth-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * (3.0 - k * k) / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(0, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // first-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(1, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // second-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k * k / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(2, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);
}

TEST_F(SphericalBesselTransformTest, ImplicitExponent)
{
    /*
     * Computes the second order spherical Bessel transform of
     * r^2*exp(-r) with input given as r^(p+2)*exp(-r) instead of
     * bare r^2*exp(-r). Compares the results with the analytic
     * expression:
     *
     *      48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *                                                          */
    int sz = 5000;
    assert(sz <= sz_max);

    double dr = 0.02;
    double rcut = dr * (sz - 1);
    double dk = PI / rcut;
    double pref = std::sqrt(2. / PI) * 48.;

    SphericalBesselTransformer sbt;
    sbt.set_fftw_plan_flag(FFTW_MEASURE);

    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * k * k / std::pow(k * k + 1, 4);
    }

    for (int p = -2; p <= 5; ++p)
    {
        for (int i = 0; i != sz; ++i)
        {
            double r = i * dr;
            f[i] = std::pow(r, 2 + p) * std::exp(-r);
        }
        sbt.radrfft(2, sz, rcut, f, g, p);
        EXPECT_LT(max_diff(sz, g_ref, g), tol);
    }
}

TEST_F(SphericalBesselTransformTest, VariableSize)
{
    /*
     * Computes the second order spherical Bessel transform of
     * r^2*exp(-r) with various input sizes. Compares the results
     * with the analytic expression:
     *
     *      48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *                                                          */
    double dr = 0.02;
    double pref = std::sqrt(2. / PI) * 48.;

    SphericalBesselTransformer sbt;
    sbt.set_fftw_plan_flag(FFTW_ESTIMATE);

    for (int sz = 5000; sz <= sz_max; sz += 1000)
    {

        double rcut = dr * (sz - 1);
        for (int i = 0; i != sz; ++i)
        {
            double r = i * dr;
            f[i] = r * r * std::exp(-r);
        }

        double dk = PI / rcut;
        for (int i = 0; i != sz; ++i)
        {
            double k = dk * i;
            g_ref[i] = pref * k * k / std::pow(k * k + 1, 4);
        }
        sbt.radrfft(2, sz, rcut, f, g, 0);
        EXPECT_LT(max_diff(sz, g_ref, g), tol);
    }
}

TEST_F(SphericalBesselTransformTest, InPlace)
{
    /*
     * Performs an in-place second order spherical Bessel transform
     * on r^2*exp(-r^2). Compares the results with the analytic
     * expression:
     *
     *      sqrt(2)/16 * k^2 * exp(-k^2/4)
     *                                                          */
    double dr = 0.02;
    double pref = std::sqrt(2.) / 16.;

    SphericalBesselTransformer sbt;
    sbt.set_fftw_plan_flag(FFTW_ESTIMATE);

    double sz = 10000;
    double rcut = dr * (sz - 1);
    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = r * r * std::exp(-r * r);
    }

    double dk = PI / rcut;
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * k * k * std::exp(-k * k / 4);
    }
    sbt.radrfft(2, sz, rcut, f, f, 0);
    EXPECT_LT(max_diff(sz, g_ref, f), tol);
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    fftw_cleanup();

    return result;
}

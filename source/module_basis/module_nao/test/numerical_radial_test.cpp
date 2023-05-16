#include "gtest/gtest.h"
#include <cmath>

#ifdef __MPI
#include <mpi.h>
#endif

#include "module_base/constants.h"
#include "module_basis/module_nao/numerical_radial.h"

using ModuleBase::PI;

/***********************************************************
 *      unit test of class "Numerical_Orbital_Lm"
 ***********************************************************/

/*! Tested functions:
 *
 *  - build
 *      - Initializes the object by setting the grid & values in one space.
 *
 *  - set_transformer
 *      - Sets a SphericalBesselTransformer for the object.
 *
 *  - set_grid
 *      - Sets up a new grid.
 *
 *  - set_uniform_grid
 *      - Sets up a new uniform grid.
 *
 *  - set_value
 *      - Updates values on an existing grid.
 *
 *  - wipe
 *      - Removes the grid & values from one space.
 *  
 *  - radtab
 *      - Computes the radial table for two-center integrals 
 *        between "this" and another object.
 *
 *  - all "getters"
 *      Gets access to private members.
 *                                                          */


class NumericalRadialTest : public ::testing::Test
{
  protected:

    void SetUp();
    void TearDown();

    int sz_max = 10000;     //!< size of each buffer
    double* grid = nullptr; //!< buffer for input grid
    double* f = nullptr;    //!< buffer for input values
    double* g = nullptr;    //!< buffer for reference values

    NumericalRadial chi; //!< object under test

    double tol = 1e-8; //!< tolerance for element-wise numerical error
};

void NumericalRadialTest::SetUp()
{
    grid = new double[sz_max];
    f = new double[sz_max];
    g = new double[sz_max];
}

void NumericalRadialTest::TearDown()
{
    delete[] f;
    delete[] g;
    delete[] grid;
}


TEST_F(NumericalRadialTest, ConstructorAndAssignment) {
    /* 
     * Tests the copy constructor and copy assignment operator.
     *                                                                      */
    double dk = PI/50;
    int sz = 10000;
    int pk = -2;
    double pref = 48 * std::sqrt(2./PI);
    for (int ik = 0; ik != sz; ++ik) {
        double k = ik * dk;
        grid[ik] = k;
        f[ik] = pref / std::pow(k*k+1, 4);
    }

    chi.build(2, false, sz, grid, f, pk);
    chi.set_uniform_grid(true, sz, PI/dk, 't');
    
    NumericalRadial chi2(chi);
    EXPECT_EQ(chi.symbol(), chi2.symbol());
    EXPECT_EQ(chi.ichi(), chi2.ichi());
    EXPECT_EQ(chi.itype(), chi2.itype());
    EXPECT_EQ(chi.l(), chi2.l());

    EXPECT_EQ(chi.nr(), chi2.nr());
    EXPECT_EQ(chi.nk(), chi2.nk());
    EXPECT_EQ(chi.rcut(), chi2.rcut());
    EXPECT_EQ(chi.kcut(), chi2.kcut());

    ASSERT_NE(chi2.ptr_rgrid(), nullptr);
    ASSERT_NE(chi2.ptr_rvalue(), nullptr);
    for (int ir = 0; ir != sz; ++ir) {
        EXPECT_EQ(chi.ptr_rgrid()[ir], chi2.ptr_rgrid()[ir]);
        EXPECT_EQ(chi.ptr_rvalue()[ir], chi2.ptr_rvalue()[ir]);
    }

    ASSERT_NE(chi2.ptr_kgrid(), nullptr);
    ASSERT_NE(chi2.ptr_kvalue(), nullptr);
    for (int ik = 0; ik != sz; ++ik) {
        EXPECT_EQ(chi.ptr_kgrid()[ik], chi2.ptr_kgrid()[ik]);
        EXPECT_EQ(chi.ptr_kvalue()[ik], chi2.ptr_kvalue()[ik]);
    }

    EXPECT_EQ(chi.pr(), chi2.pr());
    EXPECT_EQ(chi.pk(), chi2.pk());
    EXPECT_EQ(chi.is_fft_compliant(), chi2.is_fft_compliant());

    NumericalRadial chi3;
    chi3 = chi;
    EXPECT_EQ(chi.symbol(), chi3.symbol());
    EXPECT_EQ(chi.ichi(), chi3.ichi());
    EXPECT_EQ(chi.itype(), chi3.itype());
    EXPECT_EQ(chi.l(), chi3.l());

    EXPECT_EQ(chi.nr(), chi3.nr());
    EXPECT_EQ(chi.nk(), chi3.nk());
    EXPECT_EQ(chi.rcut(), chi3.rcut());
    EXPECT_EQ(chi.kcut(), chi3.kcut());

    ASSERT_NE(chi3.ptr_rgrid(), nullptr);
    ASSERT_NE(chi3.ptr_rvalue(), nullptr);
    for (int ir = 0; ir != sz; ++ir) {
        EXPECT_EQ(chi.ptr_rgrid()[ir], chi3.ptr_rgrid()[ir]);
        EXPECT_EQ(chi.ptr_rvalue()[ir], chi3.ptr_rvalue()[ir]);
    }

    ASSERT_NE(chi3.ptr_kgrid(), nullptr);
    ASSERT_NE(chi3.ptr_kvalue(), nullptr);
    for (int ik = 0; ik != sz; ++ik) {
        EXPECT_EQ(chi.ptr_kgrid()[ik], chi3.ptr_kgrid()[ik]);
        EXPECT_EQ(chi.ptr_kvalue()[ik], chi3.ptr_kvalue()[ik]);
    }

    EXPECT_EQ(chi.pr(), chi3.pr());
    EXPECT_EQ(chi.pk(), chi3.pk());
    EXPECT_EQ(chi.is_fft_compliant(), chi3.is_fft_compliant());
}


TEST_F(NumericalRadialTest, BuildAndGet) {
    /* 
     * Builds a NumericalRadial object and gets access to its members. 
     *                                                                      */
    int l = 1;
    double dr = 0.01;
    int sz = 5000;
    int pr = -1;
    int itype = 3;
    int ichi = 5;
    std::string symbol = "Au";
    for (int ir = 0; ir != sz; ++ir) {
        double r = ir * dr;
        grid[ir] = r;
        f[ir] = std::exp(-r);
    }

    chi.build(l, true, sz, grid, f, pr, itype, ichi, symbol);

    EXPECT_EQ(chi.symbol(), symbol);
    EXPECT_EQ(chi.ichi(), ichi);
    EXPECT_EQ(chi.itype(), itype);
    EXPECT_EQ(chi.l(), l);

    EXPECT_EQ(chi.nr(), sz);
    EXPECT_EQ(chi.nk(), 0);
    EXPECT_EQ(chi.rcut(), grid[sz-1]);
    EXPECT_EQ(chi.kcut(), 0);

    ASSERT_NE(chi.ptr_rgrid(), nullptr);
    ASSERT_NE(chi.ptr_rvalue(), nullptr);
    for (int ir = 0; ir != sz; ++ir) {
        EXPECT_EQ(chi.ptr_rgrid()[ir], grid[ir]);
        EXPECT_EQ(chi.ptr_rvalue()[ir], f[ir]);
    }

    EXPECT_EQ(chi.ptr_kgrid(), nullptr);
    EXPECT_EQ(chi.ptr_kvalue(), nullptr);

    EXPECT_EQ(chi.pr(), pr);
    EXPECT_EQ(chi.pk(), 0);
    EXPECT_EQ(chi.is_fft_compliant(), false);
}


TEST_F(NumericalRadialTest, SetGridAndWipe) {
    /* 
     * Builds a NumericalRadial object with r-space values of
     *
     *                  r*exp(-r)
     *
     * Sets up a k-space grid and checks whether the k-space values agree
     * with the analytic expression
     *
     *          sqrt(2/pi) * 8 * k / (k^2+1)^3.
     *
     * NOTE Currently only FFT-compliant grid is supported, so the test grids
     * are FFT-compliant.
     *
     * Finally wipe off the grid & values in r & k space.
     *                                                                      */
    double dr = 0.01;
    int sz = 5000;
    int pr = -1;
    for (int ir = 0; ir != sz; ++ir) {
        double r = ir * dr;
        grid[ir] = r;
        f[ir] = std::exp(-r);
    }

    chi.build(1, true, sz, grid, f, pr);

    double* kgrid = new double[sz];
    double dk = PI / chi.rcut();

    for (int ik = 0; ik != sz; ++ik) {
        kgrid[ik] = ik * dk;
    }

    chi.set_grid(false, sz, kgrid, 't');
    
    double pref = 8 * std::sqrt(2./PI);
    for (int ik = 0; ik != sz; ++ik)
    {
        double k = ik * dk;
        EXPECT_NEAR(pref * k / std::pow(k * k + 1, 3), chi.ptr_kvalue()[ik], tol);
    }

    EXPECT_EQ(chi.is_fft_compliant(), true);

    chi.wipe(true);
    EXPECT_EQ(chi.ptr_rgrid(), nullptr);
    EXPECT_EQ(chi.ptr_rvalue(), nullptr);
    EXPECT_EQ(chi.nr(), 0);
    EXPECT_EQ(chi.is_fft_compliant(), false);

    chi.wipe(false);
    EXPECT_EQ(chi.ptr_kgrid(), nullptr);
    EXPECT_EQ(chi.ptr_kvalue(), nullptr);
    EXPECT_EQ(chi.nk(), 0);
}


TEST_F(NumericalRadialTest, SetUniformGrid) {
    /* 
     * Builds a NumericalRadial object with k-space values of
     *
     *          48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *
     * Sets up a uniform r-space grid and checks whether the r-space values 
     * agree with the analytic expression
     *
     *          r^2 * exp(-r)
     *                                                                      */
    double dk = PI/50;
    int sz = 10000;
    int pk = -2;
    double pref = 48 * std::sqrt(2./PI);
    for (int ik = 0; ik != sz; ++ik) {
        double k = ik * dk;
        grid[ik] = k;
        f[ik] = pref / std::pow(k*k+1, 4);
    }

    chi.build(2, false, sz, grid, f, pk);
    chi.set_uniform_grid(true, sz, PI/dk, 't');
    
    double dr = PI / chi.kcut();
    for (int ir = 0; ir != sz; ++ir)
    {
        double r = ir * dr;
        EXPECT_NEAR(r*r*std::exp(-r), chi.ptr_rvalue()[ir], tol);
    }
}


//TEST_F(NumericalRadialTest, SetFFTGrid) {
//    /* 
//     * Builds a NumericalRadial object with k-space values of
//     *
//     *          48*sqrt(2/pi) * k^2  / (k^2+1)^4.
//     *
//     * Sets up a uniform r-space grid and checks whether the r-space values 
//     * agree with the analytic expression
//     *
//     *          r^2 * exp(-r)
//     *                                                                      */
//    double dk = PI/50;
//    int sz = 10000;
//    int pk = -2;
//    double pref = 48 * std::sqrt(2./PI);
//    for (int ik = 0; ik != sz; ++ik) {
//        double k = ik * dk;
//        grid[ik] = k;
//        f[ik] = pref / std::pow(k*k+1, 4);
//    }
//
//    chi.build(2, false, sz, grid, f, pk);
//    chi.set_uniform_grid(false, sz, PI/dk, 'i', true);
//    
//    double dr = PI / chi.kcut();
//    for (int ir = 100; ir != sz; ++ir)
//    {
//        double r = ir * dr;
//        EXPECT_NEAR(r*r*std::exp(-r), chi.ptr_rvalue()[ir], tol);
//    }
//}

TEST_F(NumericalRadialTest, SetValue) {
    /* 
     * Updates values in a NumericalRadial object.
     *
     *                                                                      */
    double dx = 0.01;
    int sz = 5000;
    int p = -1;
    for (int i = 0; i != sz; ++i) {
        double r = i * dx;
        grid[i] = r;
        f[i] = std::exp(-r);
    }

    chi.build(1, true, sz, grid, f, p);
    for (int ir = 0; ir != sz; ++ir) {
        f[ir] *= 2;
    }
    chi.set_value(true, f, p);

    for (int i = 0; i != sz; ++i) {
        EXPECT_EQ(chi.ptr_rvalue()[i], f[i]);
    }

    chi.build(1, false, sz, grid, f, p);
    for (int i = 0; i != sz; ++i) {
        f[i] *= 2;
    }
    chi.set_value(false, f, p);

    for (int i = 0; i != sz; ++i) {
        EXPECT_EQ(chi.ptr_kvalue()[i], f[i]);
    }

}




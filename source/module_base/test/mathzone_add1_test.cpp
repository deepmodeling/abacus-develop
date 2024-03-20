#include "../mathzone_add1.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of class Mathzone_Add1
 ***********************************************/

/**
 * - Tested Functions:
 *   - CubicSplineBoundary1
 *     - call SplineD2 and Cubic_Spline_Interpolation
 *     - to interpolate a function using the first
 *     - kind of boundary condition:f'(0) = f'(n) = 0.0
 *   - CubicSplineBoundary2
 *     - call SplineD2 and Cubic_Spline_Interpolation
 *     - to interpolate a function using the second
 *     - kind of boundary condition:f''(0) = f''(n) = 0.0
 *   - UniRadialF
 *     - call Uni_RadialF to interpolate the radial part
 *     - of an atomic orbital function, whose discrete r
 *     - points are uniform.
 *   - Factorial
 *     - calculate the factorial of an integer
 *   - DualFac
 *     - calculate the double factorial or
 *     - semifactorial of an integer
 */

class MathzoneAdd1Test : public testing::Test
{
  protected:
    const int MaxInt = 100;
    int nr_in = 17;
    int nr_out = 161;
    double* r_in;
    double* r_out;
    double* y2;
    double* psi_in;
    double* psi_out;
    double* dpsi;
    void SetUp()
    {
        r_in = new double[nr_in];
        r_out = new double[nr_out];
        y2 = new double[nr_in];
        psi_in = new double[nr_in];
        psi_out = new double[nr_out];
        dpsi = new double[nr_out];
    }
    void TearDown()
    {
        delete[] r_in;
        delete[] r_out;
        delete[] y2;
        delete[] psi_in;
        delete[] psi_out;
        delete[] dpsi;
    }
};

TEST_F(MathzoneAdd1Test, Constructor)
{
    EXPECT_NO_THROW(ModuleBase::Mathzone_Add1 MA1);
}

#define PI 3.1415926535897932
/// first kind boundary condition: f'(0) = f'(n) = 0.0
TEST_F(MathzoneAdd1Test, CubicSplineBoundary1)
{
    // data from abacus/tests/integrate/tools/PP_ORB/Si_gga_8au_60Ry_2s2p1d.orb
    // data for d orbital of Si : L = 2, N = 0
    psi_in[0] = 0;
    psi_in[1] = -2.583946346740e-01;
    psi_in[2] = -4.570087269049e-01;
    psi_in[3] = -4.374680500187e-01;
    psi_in[4] = -3.587829079989e-01;
    psi_in[5] = -2.581772323753e-01;
    psi_in[6] = -1.616203660437e-01;
    psi_in[7] = -9.108838081645e-02;
    psi_in[8] = -5.202206559586e-02;
    psi_in[9] = -3.126875315134e-02;
    psi_in[10] = -1.860199873973e-02;
    psi_in[11] = -8.049945178799e-03;
    psi_in[12] = -1.652010824028e-03;
    psi_in[13] = 1.495515249035e-03;
    psi_in[14] = 3.221037475903e-03;
    psi_in[15] = 3.802139894646e-03;
    psi_in[16] = 0;
    for (int i = 0; i < nr_in; i++)
    {
        r_in[i] = i * 0.5;
        // std::cout<< r_in[i] << " " << psi_in[i] << std::endl; // for plotting
    }
    for (int i = 0; i < nr_out; i++)
    {
        r_out[i] = i * 0.05;
    }
    ModuleBase::Mathzone_Add1::SplineD2(r_in, psi_in, nr_in, 0.0, 0.0, y2);
    // std::cout << "y2[0] "<< y2[0] << " y2[nr_in] "<< y2[nr_in-1] << std::endl; // for checking
    ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in, psi_in, y2, nr_in, r_out, nr_out, psi_out, dpsi);
    for (int i = 0; i < nr_out; i++)
    {
        int j = i / 10;
        if (i % 10 == 0)
        {
            EXPECT_EQ(psi_in[j], psi_out[i]);
        }
        // std::cout<< r_out[i] << " " << psi_out[i] << std::endl; // for plotting
    }
    EXPECT_NEAR(dpsi[0], 0.0, 1e-15);
    EXPECT_NEAR(dpsi[nr_out - 1], 0.0, 1e-15);
    // std::cout<<dpsi[0] << " " << dpsi[nr_out-1] << std::endl; // for checking
}
#define MAX(x, y) ((x) > (y) ? (x) : (y))

double count_err(double* x, std::function<double(double)> func, int n, int derivative)
{
    double maxf4x = -10000000000;
    double maxh = 0;

    for (int i = 0; i != n; ++i)
    {
        maxf4x = MAX(abs(func(x[i])), maxf4x);
    }
    // printf("maxf4x = %.8lf\n", maxf4x);
    for (int i = 1; i != n; ++i)
    {
        maxh = MAX(x[i + 1] - x[i], maxh);
    }

    // if(func==2) printf("maxh = %.8lf,maxf4x = %.8lf\n", maxh,maxf4x);
    double err = (5.0 / 384.0) * maxf4x * pow(maxh, 4 - derivative);
    return MAX(err, 1e-15);
    ;
}

/// first kind boundary condition: f'(0) = f'(n) = 0.0
TEST_F(MathzoneAdd1Test, sinx_Boundary1)
{
    for (int i = 0; i <= 10; i++)
    {
        r_in[i] = i * 0.1 * PI + 0.5 * PI;
        psi_in[i] = sin(r_in[i]);
    }
    auto f = [](double x) -> double { return sin(x); };
    double err = count_err(r_in, f, 11, 0);
    for (int i = 0; i < 100; i++)
    {
        r_out[i] = i * 0.01 * PI + 0.5 * PI;
    }
    ModuleBase::Mathzone_Add1::SplineD2(r_in, psi_in, 11, 0.0, 0.0, y2);
    ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in, psi_in, y2, 11, r_out, 100, psi_out, dpsi);
    for (int i = 0; i < 100; i++)
    {
        double ans = sin(r_out[i]);
        EXPECT_NEAR(psi_out[i], ans, err);
    }
}

/// second kind boundary condition: f''(0) = f''(n) = 0.0
TEST_F(MathzoneAdd1Test, sinx_Boundary2)
{
    for (int i = 0; i <= 10; i++)
    {
        r_in[i] = i * 0.1 * PI;
        psi_in[i] = sin(r_in[i]);
    }
    auto f = [](double x) -> double { return sin(x); };
    double err = count_err(r_in, f, 11, 0);
    for (int i = 0; i < 100; i++)
    {
        r_out[i] = i * 0.01 * PI;
    }
    ModuleBase::Mathzone_Add1::SplineD2(r_in, psi_in, 11, 0.0, 0.0, y2);
    ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in, psi_in, y2, 11, r_out, 100, psi_out, dpsi);
    for (int i = 0; i < 100; i++)
    {
        double ans = sin(r_out[i]);
        EXPECT_NEAR(psi_out[i], ans, err);
    }
}

TEST_F(MathzoneAdd1Test, expx)
{
    for (int i = 0; i <= 10; i++)
    {
        r_in[i] = i;
        psi_in[i] = exp(r_in[i]);
    }
    auto f = [](double x) -> double { return exp(x); };
    double err = count_err(r_in, f, 11, 0);
    for (int i = 0; i < 100; i++)
    {
        r_out[i] = i * 0.1;
    }
    double d2 = exp(10);
    ModuleBase::Mathzone_Add1::SplineD2(r_in, psi_in, 11, 1, d2, y2);
    ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in, psi_in, y2, 11, r_out, 100, psi_out, dpsi);
    for (int i = 0; i < 100; i++)
    {
        double ans = exp(r_out[i]);
        EXPECT_NEAR(psi_out[i], ans, err);
    }
}

// /// second kind boundary condition: f''(0) = f''(n) = 0.0
// TEST_F(MathzoneAdd1Test, expx_Boundary2)
// {
// 	for (int i=0; i<=10; i++)
// 	{
// 		r_in[i] = i;
// 		psi_in[i] = exp(r_in[i]);

// 	}
// 	double err = count_err(r_in, 2, 11, 0);
// 	for (int i=0; i< 100; i++)
// 	{
// 		r_out[i] = i*0.1;
// 	}
// 	double d2 = exp(10);
// 	ModuleBase::Mathzone_Add1::SplineD2(r_in,psi_in,11,1,d2,y2);
// 	ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in,psi_in,y2,11,r_out,100,psi_out,dpsi);
// 	for (int i=0; i< 100; i++)
// 	{
// 		double ans = exp(r_out[i]);
// 		EXPECT_NEAR(psi_out[i],ans,err);
// 	}
// }

/// second kind boundary condition: f''(0) = f''(n) = 0.0
TEST_F(MathzoneAdd1Test, CubicSplineBoundary2)
{
    // data from abacus/tests/integrate/tools/PP_ORB/Si_gga_8au_60Ry_2s2p1d.orb
    // data for 1st p orbital of Si: L = 1, N= 0
    psi_in[0] = 0;
    psi_in[1] = 2.023466616834e-01;
    psi_in[2] = 3.318755771343e-01;
    psi_in[3] = 3.648245752646e-01;
    psi_in[4] = 3.224822944963e-01;
    psi_in[5] = 2.491332695240e-01;
    psi_in[6] = 1.807635173291e-01;
    psi_in[7] = 1.266983610308e-01;
    psi_in[8] = 8.649140297968e-02;
    psi_in[9] = 5.949030687701e-02;
    psi_in[10] = 4.039513774853e-02;
    psi_in[11] = 2.778453347548e-02;
    psi_in[12] = 1.985533549037e-02;
    psi_in[13] = 1.345471632235e-02;
    psi_in[14] = 9.880871041599e-03;
    psi_in[15] = 7.795456942712e-03;
    psi_in[16] = 0;
    for (int i = 0; i < nr_in; i++)
    {
        r_in[i] = i * 0.5;
        // std::cout<< r_in[i] << " " << psi_in[i] << std::endl; // for plotting
    }
    for (int i = 0; i < nr_out; i++)
    {
        r_out[i] = i * 0.05;
    }
    ModuleBase::Mathzone_Add1::SplineD2(r_in, psi_in, nr_in, 100000.0, 100000.0, y2);
    EXPECT_EQ(y2[0], 0.0);
    EXPECT_EQ(y2[nr_in - 1], 0.0);
    // std::cout << "y2[0] "<< y2[0] << " y2[nr_in] "<< y2[nr_in-1] << std::endl; // for checking
    ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(r_in, psi_in, y2, nr_in, r_out, nr_out, psi_out, dpsi);
    for (int i = 0; i < nr_out; i++)
    {
        int j = i / 10;
        if (i % 10 == 0)
        {
            EXPECT_EQ(psi_in[j], psi_out[i]);
        }
        // std::cout<< r_out[i] << " " << psi_out[i] << std::endl; // for plotting
    }
    // std::cout<<dpsi[0] << " " << dpsi[nr_out-1] << std::endl; // for checking
}

TEST_F(MathzoneAdd1Test, UniRadialF)
{
    // data from abacus/tests/integrate/tools/PP_ORB/Si_gga_8au_60Ry_2s2p1d.orb
    // data for 1st p orbital of Si: L = 1, N= 0
    psi_in[0] = 0;
    psi_in[1] = 2.023466616834e-01;
    psi_in[2] = 3.318755771343e-01;
    psi_in[3] = 3.648245752646e-01;
    psi_in[4] = 3.224822944963e-01;
    psi_in[5] = 2.491332695240e-01;
    psi_in[6] = 1.807635173291e-01;
    psi_in[7] = 1.266983610308e-01;
    psi_in[8] = 8.649140297968e-02;
    psi_in[9] = 5.949030687701e-02;
    psi_in[10] = 4.039513774853e-02;
    psi_in[11] = 2.778453347548e-02;
    psi_in[12] = 1.985533549037e-02;
    psi_in[13] = 1.345471632235e-02;
    psi_in[14] = 9.880871041599e-03;
    psi_in[15] = 7.795456942712e-03;
    psi_in[16] = 0;
    for (int i = 0; i < nr_in; i++)
    {
        r_in[i] = i * 0.5;
        // std::cout<< r_in[i] << " " << psi_in[i] << std::endl; // for plotting
    }
    double dr = 0.5;
    for (int i = 0; i < nr_out; i++)
    {
        int j = i / 10;
        r_out[i] = i * 0.05;
        psi_out[i] = ModuleBase::Mathzone_Add1::Uni_RadialF(psi_in, nr_in, dr, r_out[i]);
        if (i % 10 == 0)
        {
            EXPECT_NEAR(psi_in[j], psi_out[i], 1e-15);
        }
        // std::cout<< r_out[i] << " " << psi_out[i] << std::endl; // for plotting
    }
}

TEST_F(MathzoneAdd1Test, Factorial)
{
    double fac[MaxInt];
    double fac1;
    fac[1] = fac[0] = 1.0;
    fac1 = ModuleBase::Mathzone_Add1::factorial(0);
    EXPECT_EQ(fac[0], fac1);
    fac1 = ModuleBase::Mathzone_Add1::factorial(1);
    EXPECT_EQ(fac[1], fac1);
    for (int i = 2; i < MaxInt; i++)
    {
        fac[i] = i * fac[i - 1];
        fac1 = ModuleBase::Mathzone_Add1::factorial(i);
        EXPECT_EQ(fac[i], fac1);
    }
}

TEST_F(MathzoneAdd1Test, DualFac)
{
    double dualfac[MaxInt];
    double dualfacm1;
    double dualfac1;
    dualfacm1 = ModuleBase::Mathzone_Add1::dualfac(-1);
    EXPECT_EQ(dualfacm1, 1.0);
    dualfac[0] = 1.0;
    dualfac1 = ModuleBase::Mathzone_Add1::dualfac(0);
    EXPECT_EQ(dualfac[0], dualfac1);
    dualfac[1] = 1.0;
    dualfac1 = ModuleBase::Mathzone_Add1::dualfac(1);
    EXPECT_EQ(dualfac[1], dualfac1);
    for (int i = 2; i < MaxInt; i++)
    {
        dualfac[i] = i * dualfac[i - 2];
        dualfac1 = ModuleBase::Mathzone_Add1::dualfac(i);
        EXPECT_EQ(dualfac[i], dualfac1);
    }
}

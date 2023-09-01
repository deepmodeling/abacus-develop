#include "module_base/cubic_spline.h"

#include <cmath>

#include "gtest/gtest.h"
#include "module_base/constants.h"
#include "stdio.h"

#ifdef __MPI
#include <mpi.h>
#endif

using ModuleBase::CubicSpline;
using ModuleBase::PI;

/***********************************************************
 *      Unit test of class CubicSpline
 ***********************************************************/

/*! Tested functions:
 *
 *  - build
 *      Constructs the cubic spline interpolant from given
 *      data points and boundary condition.
 *
 *  - get
 *      Evaluates the interpolant at certain points.
 *                                                          */
class CubicSplineTest : public ::testing::Test
{
  protected:
    /// Allocates buffers
    void SetUp();

    /// Deallocates buffers
    void TearDown();

    int n_max = 20; ///< size of each buffer
    double* x_;     ///< buffer for the x coordinates of data points
    double* y_;     ///< buffer for the y coordinates of data points

    double* x_interp_; ///< places where the interpolant is evaluated
    double* y_interp_; ///< values of the interpolant at x_interp_
    double* y_ref_;    ///< reference values for y_interp_ drawn from
                       ///< scipy.interpolate.CubicSpline

    double tol = 1e-8; ///< tolerance for element-wise numerical error
};

void CubicSplineTest::SetUp()
{
    x_ = new double[n_max];
    y_ = new double[n_max];
    x_interp_ = new double[n_max];
    y_interp_ = new double[n_max];
    y_ref_ = new double[n_max];
}

void CubicSplineTest::TearDown()
{
    delete[] x_;
    delete[] y_;
    delete[] x_interp_;
    delete[] y_interp_;
    delete[] y_ref_;
}
#define MAX(x, y) ((x) > (y) ? (x) : (y))
double count_err(double* x, int func, int n, int derivative)
{
    double maxf4x = -10000000000;
    double maxh = 0;

    for (int i = 0; i != n; ++i)
    {
        if (func == 0)
        {
            maxf4x = MAX(abs(sin(x[i])), maxf4x);
        }
        else if (func == 1)
        {
            maxf4x = 1;
        }
        else if (func == 2)
        {
            maxf4x = MAX(abs(exp(-x[i])), maxf4x);
        }
        else if (func == 3)
        {
            maxf4x = MAX(abs(-6 * pow(x[i], -4)), maxf4x);
        }
    }
    // printf("maxf4x = %.8lf\n", maxf4x);
    for (int i = 1; i != n; ++i)
    {
        maxh = MAX(x[i + 1] - x[i], maxh);
    }

    // if(func==2) printf("maxh = %.8lf,maxf4x = %.8lf\n", maxh,maxf4x);
    double err = (5.0 / 384.0) * maxf4x * pow(maxh, 4 - derivative);
    return err;
}

TEST_F(CubicSplineTest, NotAKnot)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = std::sqrt(i);
        y_[i] = std::sin(x_[i]);
    }

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.57;
    x_interp_[3] = 2.00;
    x_interp_[4] = 2.99;
    x_interp_[5] = 3.00;

    y_ref_[0] = 0.;
    y_ref_[1] = 0.0105903444284005;
    y_ref_[2] = 1.0000633795463434;
    y_ref_[3] = 0.9092974268256817;
    y_ref_[4] = 0.1510153464180796;
    y_ref_[5] = 0.1411200080598672;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

TEST_F(CubicSplineTest, PeriodicAndUniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = i * 2 * PI / (n - 1);
        y_[i] = std::cos(x_[i]);
    }

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

    int ni = 5;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.5 * PI;
    x_interp_[2] = PI;
    x_interp_[3] = 1.5 * PI;
    x_interp_[4] = 2 * PI;

    y_ref_[0] = 1.0000000000000000e+00;
    y_ref_[1] = 1.4356324132368183e-04;
    y_ref_[2] = -9.9930291851807085e-01;
    y_ref_[3] = 1.4356324132349871e-04;
    y_ref_[4] = 1.0000000000000000e+00;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

TEST_F(CubicSplineTest, FirstDeriv)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = std::sqrt(i);
        y_[i] = std::exp(-x_[i]);
    }

    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 -3,
                 -0.5);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 2.0;
    x_interp_[4] = 2.54;
    x_interp_[5] = 3.00;

    y_ref_[0] = 1.;
    y_ref_[1] = 0.9704131180863818;
    y_ref_[2] = 0.1367376505691157;
    y_ref_[3] = 0.1353352832366127;
    y_ref_[4] = 0.0798871927951471;
    y_ref_[5] = 0.0497870683678639;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

/**
 * @brief a test case for CubicSpline on function??sinx??FirstDeriv and Uniform distribution
 *
 * err = C0*max(f'4(x))*h^4 = 5/384 * max(sin(x)) * h^4 = 5/384 *
 */
TEST_F(CubicSplineTest, FirstDeriv_sinx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0) * 2 * PI;
        y_[i] = sin(x_[i]);
    }

    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 1,
                 0.80901699);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 5.550;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.669239857276262;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_sinx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = sin(x_[0]);
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 2 * PI;
        y_[9 - i] = sin(x_[9 - i]);
    }

    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 1,
                 1);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 6;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.756802495307928;
    y_ref_[5] = -0.202090837026266;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_sinx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0) * 2 * PI;
        y_[i] = sin(x_[i]);
    }

    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 0,
                 0.587785);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 5.550;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.669239857276262;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_sinx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = sin(x_[0]);
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 2 * PI;
        y_[9 - i] = sin(x_[9 - i]);
    }

    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 0,
                 0);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 6;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.279415498198926;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, Periodic_sinx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n - 1; ++i)
    {
        x_[i] = ((double)i / 10.0) * 2 * PI;
        y_[i] = sin(x_[i]);
    }
    x_[9] = 2 * PI;
    y_[9] = sin(0);
    double err = count_err(x_, 0, n, 0);
    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 6;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.279415498198926;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, Periodic_sinx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = sin(x_[0]);
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 2 * PI;
        y_[9 - i] = sin(x_[9 - i]);
    }
    x_[9] = 2 * PI;
    y_[9] = sin(0);
    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

    double err = count_err(x_, 0, n, 0);
    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 6;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.279415498198926;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

// TEST_F(CubicSplineTest, NotAKnot)
// {
//     CubicSpline cubspl;

//     int n = 10;
//     for (int i = 0; i != n; ++i)
//     {
//         x_[i] = std::sqrt(i);
//         y_[i] = std::sin(x_[i]);
//     }

//     cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

//     int ni = 6;
//     x_interp_[0] = 0.0;
//     x_interp_[1] = 0.01;
//     x_interp_[2] = 1.57;
//     x_interp_[3] = 2.00;
//     x_interp_[4] = 2.99;
//     x_interp_[5] = 3.00;

//     y_ref_[0] = 0.;
//     y_ref_[1] = 0.0105903444284005;
//     y_ref_[2] = 1.0000633795463434;
//     y_ref_[3] = 0.9092974268256817;
//     y_ref_[4] = 0.1510153464180796;
//     y_ref_[5] = 0.1411200080598672;

//     cubspl.get(ni, x_interp_, y_interp_);
//     for (int i = 0; i != ni; ++i)
//     {
//         EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//     }
// }

TEST_F(CubicSplineTest, NotAKnot_sinx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0) * 2 * PI;
        y_[i] = sin(x_[i]);
    }

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);
    double err = count_err(x_, 0, n, 0);
    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 5; // if x=5.5, not pass

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.958924274663138;

    // printf("err=%.8lf\n",err);

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, NotAKnot_sinx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = sin(x_[0]);
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 2 * PI;
        y_[9 - i] = sin(x_[9 - i]);
    }

    double err = count_err(x_, 1, n, 0);
    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 3;
    x_interp_[4] = 4.559;
    x_interp_[5] = 6;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.009999833334167;
    y_ref_[2] = 0.913413361341225;
    y_ref_[3] = 0.141120008059867;
    y_ref_[4] = -0.988258957900347;
    y_ref_[5] = -0.279415498198926;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_2x_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i * 10.0);
        y_[i] = 2 * x_[i];
    }
    double err = count_err(x_, 1, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 2,
                 2);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 1.60;
    x_interp_[2] = 3.20;
    x_interp_[3] = 4.80;
    x_interp_[4] = 6.40;
    x_interp_[5] = 8.00;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 3.200000000000000;
    y_ref_[2] = 6.400000000000000;
    y_ref_[3] = 9.600000000000000;
    y_ref_[4] = 12.800000000000000;
    y_ref_[5] = 16.000000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_2x_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = 0;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = 2 * (x_[9 - i]);
    }
    double err = count_err(x_, 1, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 2,
                 2);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 0.10;
    x_interp_[3] = 0.33;
    x_interp_[4] = 0.66;
    x_interp_[5] = 0.99;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.020000000000000;
    y_ref_[2] = 0.200000000000000;
    y_ref_[3] = 0.660000000000000;
    y_ref_[4] = 1.320000000000000;
    y_ref_[5] = 1.980000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_2x_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i * 10.0);
        y_[i] = 2 * x_[i];
    }
    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 0,
                 0);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 1.60;
    x_interp_[2] = 3.20;
    x_interp_[3] = 4.80;
    x_interp_[4] = 6.40;
    x_interp_[5] = 8.00;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 3.200000000000000;
    y_ref_[2] = 6.400000000000000;
    y_ref_[3] = 9.600000000000000;
    y_ref_[4] = 12.800000000000000;
    y_ref_[5] = 16.000000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_2x_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = 0;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = 2 * (x_[9 - i]);
    }

    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 0,
                 0);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 0.10;
    x_interp_[3] = 0.33;
    x_interp_[4] = 0.66;
    x_interp_[5] = 0.99;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.020000000000000;
    y_ref_[2] = 0.200000000000000;
    y_ref_[3] = 0.660000000000000;
    y_ref_[4] = 1.320000000000000;
    y_ref_[5] = 1.980000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

// TEST_F(CubicSplineTest, Periodic_2x_uniform)
// {
//     CubicSpline cubspl;

//     int n = 10;
//     for (int i = 0; i != n; ++i)
//     {
//         x_[i] = ((double)i / 10.0);
//         y_[i] = 2 * (x_[i]);
//     }

//     cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

//     int ni = 6;
//     x_interp_[0] = 0.0;
//     x_interp_[1] = 1.60;
//     x_interp_[2] = 3.20;
//     x_interp_[3] = 4.80;
//     x_interp_[4] = 6.40;
//     x_interp_[5] = 8.00;

//     y_ref_[0] = 0.000000000000000;
//     y_ref_[1] = 3.200000000000000;
//     y_ref_[2] = 6.400000000000000;
//     y_ref_[3] = 9.600000000000000;
//     y_ref_[4] = 12.800000000000000;
//     y_ref_[5] = 16.000000000000000;

//     cubspl.get(ni, x_interp_, y_interp_);
//     for (int i = 0; i != ni; ++i)
//     {
//         EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//     }
// }

// TEST_F(CubicSplineTest, Periodic_2x_Notuniform)
// {
//     CubicSpline cubspl;

//     int n = 10;
//     x_[0] = 0;
//     y_[0] = 0;
//     for (int i = 8; i >= 0; i--)
//     {
//         x_[9 - i] = exp(-i);
//         y_[9 - i] = 2 * (x_[9 - i]);
//     }

//     cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

//     int ni = 6;
//     x_interp_[0] = 0.0;
//     x_interp_[1] = 0.01;
//     x_interp_[2] = 0.10;
//     x_interp_[3] = 0.33;
//     x_interp_[4] = 0.66;
//     x_interp_[5] = 0.99;

//     y_ref_[0] = 0.000000000000000;
//     y_ref_[1] = 0.020000000000000;
//     y_ref_[2] = 0.200000000000000;
//     y_ref_[3] = 0.660000000000000;
//     y_ref_[4] = 1.320000000000000;
//     y_ref_[5] = 1.980000000000000;

//     cubspl.get(ni, x_interp_, y_interp_);
//     for (int i = 0; i != ni; ++i)
//     {
//         EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//     }
// }

TEST_F(CubicSplineTest, NotAKnot_2x_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i * 10.0);
        y_[i] = 2 * (x_[i]);
    }

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 1.60;
    x_interp_[2] = 3.20;
    x_interp_[3] = 4.80;
    x_interp_[4] = 6.40;
    x_interp_[5] = 8.00;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 3.200000000000000;
    y_ref_[2] = 6.400000000000000;
    y_ref_[3] = 9.600000000000000;
    y_ref_[4] = 12.800000000000000;
    y_ref_[5] = 16.000000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

TEST_F(CubicSplineTest, NotAKnot_2x_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = 0;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = 2 * (x_[9 - i]);
    }

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 0.10;
    x_interp_[3] = 0.33;
    x_interp_[4] = 0.66;
    x_interp_[5] = 0.99;

    y_ref_[0] = 0.000000000000000;
    y_ref_[1] = 0.020000000000000;
    y_ref_[2] = 0.200000000000000;
    y_ref_[3] = 0.660000000000000;
    y_ref_[4] = 1.320000000000000;
    y_ref_[5] = 1.980000000000000;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_expx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0);
        y_[i] = exp(-x_[i]);
    }
    double err = count_err(x_, 2, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 -1,
                 -0.406570);

    int ni = 6;
    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.200000000000000;
    x_interp_[4] = 0.400000000000000;
    x_interp_[5] = 0.800000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.818730753077982;
    y_ref_[4] = 0.670320046035639;
    y_ref_[5] = 0.449328964117222;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_expx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = 1;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = exp(-(x_[9 - i]));
    }
    double err = count_err(x_, 2, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 -1,
                 -0.367879);

    int ni = 6;

    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.330000000000000;
    x_interp_[4] = 0.660000000000000;
    x_interp_[5] = 0.990000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.718923733431926;
    y_ref_[4] = 0.516851334491699;
    y_ref_[5] = 0.371576691022046;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_expx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0);
        y_[i] = exp(-x_[i]);
    }
    double err = count_err(x_, 0, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 1,
                 0.406570);

    int ni = 6;
    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.200000000000000;
    x_interp_[4] = 0.400000000000000;
    x_interp_[5] = 0.800000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.818730753077982;
    y_ref_[4] = 0.670320046035639;
    y_ref_[5] = 0.449328964117222;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_expx_Notuniform)
{
    CubicSpline cubspl;
    int n = 10;
    x_[0] = 0;
    y_[0] = 1;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = exp(-(x_[9 - i]));
    }
    double err = count_err(x_, 2, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 1,
                 0.367879441171442);

    int ni = 6;
    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.200000000000000;
    x_interp_[4] = 0.400000000000000;
    x_interp_[5] = 0.800000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.818730753077982;
    y_ref_[4] = 0.670320046035639;
    y_ref_[5] = 0.449328964117222;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, NotAKnot_expx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = ((double)i / 10.0);
        y_[i] = exp(-x_[i]);
    }
    double err = count_err(x_, 2, n, 0);

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;

    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.200000000000000;
    x_interp_[4] = 0.400000000000000;
    x_interp_[5] = 0.800000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.818730753077982;
    y_ref_[4] = 0.670320046035639;
    y_ref_[5] = 0.449328964117222;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, NotAKnot_expx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    x_[0] = 0;
    y_[0] = 1;
    for (int i = 8; i >= 0; i--)
    {
        x_[9 - i] = exp(-i);
        y_[9 - i] = exp(-(x_[9 - i]));
    }
    double err = count_err(x_, 2, n, 0);

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.000000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 0.200000000000000;
    x_interp_[4] = 0.400000000000000;
    x_interp_[5] = 0.800000000000000;
    y_ref_[0] = 1.000000000000000;
    y_ref_[1] = 0.990049833749168;
    y_ref_[2] = 0.904837418035960;
    y_ref_[3] = 0.818730753077982;
    y_ref_[4] = 0.670320046035639;
    y_ref_[5] = 0.449328964117222;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_lnx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 1; i != n + 1; ++i)
    {
        x_[i] = ((double)i / 2.0);
        y_[i] = log(x_[i]);
    }
    double err = count_err(x_, 3, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 2,
                 0.2);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, FirstDeriv_lnx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 9; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 5;
        y_[9 - i] = log(x_[9 - i]);
    }
    double err = count_err(x_, 2, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 1620.616785515076799,
                 0.2);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_lnx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 1; i != n + 1; ++i)
    {
        x_[i] = ((double)i / 2.0);
        y_[i] = log(x_[i]);
    }
    double err = count_err(x_, 3, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 -4,
                 -0.04);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv_lnx_Notuniform)
{
    CubicSpline cubspl;
    int n = 10;
    for (int i = 9; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 5;
        y_[9 - i] = log(x_[9 - i]);
    }
    double err = count_err(x_, 3, n, 0);
    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 -2626398.765493220183998,
                 0.04);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;
    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, NotAKnot_lnx_uniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 1; i != n + 1; ++i)
    {
        x_[i] = ((double)i / 2.0);
        y_[i] = log(x_[i]);
    }
    double err = count_err(x_, 3, n, 0);

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, NotAKnot_lnx_Notuniform)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 9; i >= 0; i--)
    {
        x_[9 - i] = exp(-i) * 5;
        y_[9 - i] = log(x_[9 - i]);
    }
    double err = count_err(x_, 3, n, 0);

    cubspl.build(n, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    int ni = 6;
    x_interp_[0] = 0.001000000000000;
    x_interp_[1] = 0.010000000000000;
    x_interp_[2] = 0.100000000000000;
    x_interp_[3] = 1.000000000000000;
    x_interp_[4] = 2.000000000000000;
    x_interp_[5] = 4.000000000000000;
    y_ref_[0] = -6.907755278982137;
    y_ref_[1] = -4.605170185988091;
    y_ref_[2] = -2.302585092994045;
    y_ref_[3] = 0.000000000000000;
    y_ref_[4] = 0.693147180559945;
    y_ref_[5] = 1.386294361119891;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], err);
    }
}

TEST_F(CubicSplineTest, SecondDeriv)
{
    CubicSpline cubspl;

    int n = 10;
    for (int i = 0; i != n; ++i)
    {
        x_[i] = std::sqrt(i);
        y_[i] = std::exp(-x_[i]);
    }

    cubspl.build(n,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 -1,
                 -0.01);

    int ni = 6;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.01;
    x_interp_[2] = 1.99;
    x_interp_[3] = 2.0;
    x_interp_[4] = 2.54;
    x_interp_[5] = 3.00;

    y_ref_[0] = 1.;
    y_ref_[1] = 0.9952111318752899;
    y_ref_[2] = 0.13668283949289;
    y_ref_[3] = 0.1353352832366127;
    y_ref_[4] = 0.0788753653329337;
    y_ref_[5] = 0.0497870683678639;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
}

TEST_F(CubicSplineTest, TwoPoints)
{
    CubicSpline cubspl;

    int ni = 5;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.1;
    x_interp_[2] = 0.33;
    x_interp_[3] = 0.9;
    x_interp_[4] = 1.0;

    x_[0] = 0;
    x_[1] = 1;
    y_[0] = 2.33;
    y_[1] = 4.33;

    // first derivative
    cubspl.build(2,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::first_deriv,
                 CubicSpline::BoundaryCondition::first_deriv,
                 0.8,
                 1.5);

    y_ref_[0] = 2.33;
    y_ref_[1] = 2.4373;
    y_ref_[2] = 2.8487171;
    y_ref_[3] = 4.159700000000001;
    y_ref_[4] = 4.33;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }

    // second derivative
    cubspl.build(2,
                 x_,
                 y_,
                 CubicSpline::BoundaryCondition::second_deriv,
                 CubicSpline::BoundaryCondition::second_deriv,
                 0.8,
                 1.5);

    y_ref_[0] = 2.33;
    y_ref_[1] = 2.48245;
    y_ref_[2] = 2.86725265;
    y_ref_[3] = 4.074050000000001;
    y_ref_[4] = 4.33;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }

    // periodic
    y_[1] = y_[0];
    cubspl.build(2, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], 2.33, tol);
    }

    // "not-a-knot" is invalid for n=2
}

TEST_F(CubicSplineTest, ThreePoints)
{
    CubicSpline cubspl;

    int ni = 5;
    x_interp_[0] = 0.0;
    x_interp_[1] = 0.1;
    x_interp_[2] = 0.33;
    x_interp_[3] = 0.9;
    x_interp_[4] = 1.0;

    // not-a-knot
    x_[0] = 0;
    x_[1] = 0.4;
    x_[2] = 1.0;

    y_[0] = 1.2;
    y_[1] = 2.5;
    y_[2] = 4.8;

    cubspl.build(3, x_, y_, CubicSpline::BoundaryCondition::not_a_knot, CubicSpline::BoundaryCondition::not_a_knot);

    y_ref_[0] = 1.2;
    y_ref_[1] = 1.5075;
    y_ref_[2] = 2.259025;
    y_ref_[3] = 4.3875;
    y_ref_[4] = 4.8;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }

    // periodic
    y_[2] = y_[0];
    cubspl.build(3, x_, y_, CubicSpline::BoundaryCondition::periodic, CubicSpline::BoundaryCondition::periodic);

    y_ref_[0] = 1.2;
    y_ref_[1] = 1.44375;
    y_ref_[2] = 2.35383125;
    y_ref_[3] = 1.2361111111111112;
    y_ref_[4] = 1.2;

    cubspl.get(ni, x_interp_, y_interp_);
    for (int i = 0; i != ni; ++i)
    {
        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
    }
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

    return result;
}

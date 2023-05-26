#include "module_base/cubic_spline.h"

#include <cmath>

#include "gtest/gtest.h"
#include "module_base/constants.h"

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

    double tol = 1e-15; ///< tolerance for element-wise numerical error
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

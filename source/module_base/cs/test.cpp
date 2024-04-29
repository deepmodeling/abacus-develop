#include "cubic_spline.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>

#include "gtest/gtest.h"

using BoundaryCondition = CubicSpline::BoundaryCondition;
using BoundaryType = CubicSpline::BoundaryType;

/**
 * @brief Unit test of class CubicSpline
 *
 * Tested functions include:
 *
 *  - build
 *      Constructs a cubic spline interpolant from
 *      a set of data points and boundary conditions.
 *
 *  - eval
 *      Evaluates a single interpolant at multiple places.
 *
 *  - add
 *      Adds an interpolant that shares the same knots.
 *
 *  - multi_eval
 *      Evaluates multiple interpolants at a single place.
 *
 *  - reserve
 *      Reserves memory for multiple interpolants.
 *
 *  - heap_usage
 *      Returns the heap usage of the object.
 *
 */
class CubicSplineTest : public ::testing::Test
{
protected:

    CubicSplineTest():
        n_max_(1000),
        spline_(3 * n_max_),
        x_(spline_.data()),
        y_(x_ + n_max_),
        dy_(y_ + n_max_),
        n_interp_max_(1000),
        interp_(7 * n_interp_max_),
        x_interp_(interp_.data()),
        y_interp_(x_interp_ + n_interp_max_),
        dy_interp_(y_interp_ + n_interp_max_),
        d2y_interp_(dy_interp_ + n_interp_max_),
        y_ref_(d2y_interp_ + n_interp_max_),
        dy_ref_(y_ref_ + n_interp_max_),
        d2y_ref_(dy_ref_ + n_interp_max_)
    {}

    ~CubicSplineTest() = default;

    /// maximum number of knots
    int n_max_;

    /// buffer for a cubic spline (x_, y_ & dy_)
    std::vector<double> spline_;

    /// knots (x-coordinates of data points)
    double* x_;

    /// values at knots (y-coordinates of data points)
    double* y_;

    /// derivatives at knots (computed when building a cubic spline)
    double* dy_;

    /// maximum number of places to evaluate an interpolant
    int n_interp_max_;

    /// buffer for interpolant evaluation
    std::vector<double> interp_;

    /// places to evaluate an interpolant
    double* x_interp_;

    /// values and derivatives of the interpolant at x_interp_
    double* y_interp_;
    double* dy_interp_;
    double* d2y_interp_;

    /// reference values and derivatives
    double* y_ref_;
    double* dy_ref_;
    double* d2y_ref_;

    /// theoretical error bound for complete cubic spline
    double error_bound(
        int n,
        const double* x,
        const std::function<double(double)>& f,
        int d = 0
    ) const;

    ///
    void read(
        const std::string& fname,
        int& n,
        double* x,
        double* y,
        BoundaryCondition& bc_start,
        BoundaryCondition& bc_end,
        int& n_interp,
        double* x_interp,
        double* y_interp,
        double* dy_interp,
        double* d2y_interp
    ) const;

    /// tolerance for (accumulated) rounding errors & cross-check
    constexpr static double tol_ = 1e-12; 
};


double CubicSplineTest::error_bound(
    int n,
    const double* x,
    const std::function<double(double)>& d4f,
    int d
) const
{
    std::vector<double> buffer(n);

    std::adjacent_difference(x, x + n, buffer.begin());
    double max_dx = *std::max_element(buffer.begin() + 1, buffer.end());

    auto d4f_abs = [&d4f](double x) { return std::abs(d4f(x)); };
    std::transform(x, x + n, buffer.begin(), d4f_abs);
    double max_d4f = *std::max_element(buffer.begin(), buffer.end());

    // See Carl de Boor, "A Practical Guide to Splines", Chapter V.
    switch (d)
    {
        case 0:
            return 5.0 / 384.0 * std::pow(max_dx, 4) * max_d4f;
        case 1:
            return 1.0 / 24.0 * std::pow(max_dx, 3) * max_d4f;
        case 2:
            return 3.0 / 8.0 * std::pow(max_dx, 2) * max_d4f;
        default:
            assert(false); // should not reach here
    }
}


void CubicSplineTest::read(
    const std::string& fname,
    int& n,
    double* x,
    double* y,
    BoundaryCondition& bc_start,
    BoundaryCondition& bc_end,
    int& n_interp,
    double* x_interp,
    double* y_interp,
    double* dy_interp,
    double* d2y_interp
) const
{
    std::ifstream ifs(fname);
    assert(ifs.is_open());

    std::string line, bc1, bc2;

    // boundary condition
    std::getline(ifs, line);
    std::stringstream ss(line);
    ss >> bc1 >> bc2;

    auto bc_parse = [](const std::string& bc)
    {
        if (bc == "periodic")
        {
            return BoundaryCondition(BoundaryType::periodic);
        }
        if (bc == "not_a_knot")
        {
            return BoundaryCondition(BoundaryType::not_a_knot);
        }
        if (bc.find("first_deriv") != std::string::npos)
        {
            return BoundaryCondition(BoundaryType::first_deriv,
                                     std::stod(bc.substr(12, std::string::npos)));
        }
        if (bc.find("second_deriv") != std::string::npos)
        {
            return BoundaryCondition(BoundaryType::second_deriv,
                                     std::stod(bc.substr(13, std::string::npos)));
        }
        else
        {
            assert(false);
        }
    };

    bc_start = bc_parse(bc1);
    bc_end = bc_parse(bc2);

    double* data[6] = {x, y, x_interp, y_interp, dy_interp, d2y_interp};
    for (int i = 0; i < 6; ++i)
    {
        std::getline(ifs, line);
        std::stringstream ss(line);
        data[i] = std::copy(std::istream_iterator<double>(ss),
                            std::istream_iterator<double>(), data[i]);
    }
    n = data[0] - x;
    n_interp = data[2] - x_interp;
}



//TEST_F(CubicSplineTest, NotAKnot)
//{
//    int n = 10;
//    for (int i = 0; i != n; ++i)
//    {
//        x_[i] = std::sqrt(i);
//        y_[i] = std::sin(x_[i]);
//    }
//
//    CubicSpline cubspl(n, x_, y_);
//
//    int n_interp = 6;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.01;
//    x_interp_[2] = 1.57;
//    x_interp_[3] = 2.00;
//    x_interp_[4] = 2.99;
//    x_interp_[5] = 3.00;
//
//    y_ref_[0] = 0.;
//    y_ref_[1] = 0.0105903444284005;
//    y_ref_[2] = 1.0000633795463434;
//    y_ref_[3] = 0.9092974268256817;
//    y_ref_[4] = 0.1510153464180796;
//    y_ref_[5] = 0.1411200080598672;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//
//    // static member function
//    CubicSpline::build(n, x_, y_, {}, {}, dy_);
//    CubicSpline::eval(n, x_, y_, dy_, n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//}
//
//
//TEST_F(CubicSplineTest, Periodic)
//{
//    int n = 10;
//    double pi = std::acos(-1.0);
//    for (int i = 0; i != n; ++i)
//    {
//        x_[i] = i * 2 * pi / (n - 1);
//        y_[i] = std::cos(x_[i]);
//    }
//
//    CubicSpline cubspl(n,
//                       x_,
//                       y_,
//                       CubicSpline::BoundaryType::periodic,
//                       CubicSpline::BoundaryType::periodic);
//
//    int n_interp = 5;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.5 * pi;
//    x_interp_[2] = pi;
//    x_interp_[3] = 1.5 * pi;
//    x_interp_[4] = 2 * pi;
//
//    y_ref_[0] = 1.0000000000000000e+00;
//    y_ref_[1] = 1.4356324132368183e-04;
//    y_ref_[2] = -9.9930291851807085e-01;
//    y_ref_[3] = 1.4356324132349871e-04;
//    y_ref_[4] = 1.0000000000000000e+00;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//}
//
//
//TEST_F(CubicSplineTest, FirstDeriv)
//{
//    int n = 10;
//    for (int i = 0; i != n; ++i)
//    {
//        x_[i] = std::sqrt(i);
//        y_[i] = std::exp(-x_[i]);
//    }
//
//    CubicSpline cubspl(n,
//                       x_,
//                       y_,
//                       {CubicSpline::BoundaryType::first_deriv, -3.0},
//                       {CubicSpline::BoundaryType::first_deriv, -0.5});
//
//    int n_interp = 6;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.01;
//    x_interp_[2] = 1.99;
//    x_interp_[3] = 2.0;
//    x_interp_[4] = 2.54;
//    x_interp_[5] = 3.00;
//
//    y_ref_[0] = 1.;
//    y_ref_[1] = 0.9704131180863818;
//    y_ref_[2] = 0.1367376505691157;
//    y_ref_[3] = 0.1353352832366127;
//    y_ref_[4] = 0.0798871927951471;
//    y_ref_[5] = 0.0497870683678639;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//}
//
//
//TEST_F(CubicSplineTest, SecondDeriv)
//{
//    int n = 10;
//    for (int i = 0; i != n; ++i)
//    {
//        x_[i] = std::sqrt(i);
//        y_[i] = std::exp(-x_[i]);
//    }
//
//    CubicSpline cubspl(n,
//                       x_,
//                       y_,
//                       {CubicSpline::BoundaryType::second_deriv, -1},
//                       {CubicSpline::BoundaryType::second_deriv, -0.01});
//
//    int ni = 6;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.01;
//    x_interp_[2] = 1.99;
//    x_interp_[3] = 2.0;
//    x_interp_[4] = 2.54;
//    x_interp_[5] = 3.00;
//
//    y_ref_[0] = 1.;
//    y_ref_[1] = 0.9952111318752899;
//    y_ref_[2] = 0.13668283949289;
//    y_ref_[3] = 0.1353352832366127;
//    y_ref_[4] = 0.0788753653329337;
//    y_ref_[5] = 0.0497870683678639;
//
//    cubspl.eval(ni, x_interp_, y_interp_);
//    for (int i = 0; i != ni; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//}
//
//
//TEST_F(CubicSplineTest, TwoPoints)
//{
//    int n_interp = 5;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.1;
//    x_interp_[2] = 0.33;
//    x_interp_[3] = 0.9;
//    x_interp_[4] = 1.0;
//
//    x_[0] = 0;
//    x_[1] = 1;
//    y_[0] = 2.33;
//    y_[1] = 4.33;
//
//    // first derivative
//    CubicSpline cubspl(2,
//                       x_,
//                       y_,
//                       {CubicSpline::BoundaryType::first_deriv, 0.8},
//                       {CubicSpline::BoundaryType::first_deriv, 1.5});
//
//    y_ref_[0] = 2.33;
//    y_ref_[1] = 2.4373;
//    y_ref_[2] = 2.8487171;
//    y_ref_[3] = 4.159700000000001;
//    y_ref_[4] = 4.33;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//
//    // second derivative
//    cubspl = CubicSpline(2,
//                         x_,
//                         y_,
//                         {CubicSpline::BoundaryType::second_deriv, 0.8},
//                         {CubicSpline::BoundaryType::second_deriv, 1.5});
//
//    y_ref_[0] = 2.33;
//    y_ref_[1] = 2.48245;
//    y_ref_[2] = 2.86725265;
//    y_ref_[3] = 4.074050000000001;
//    y_ref_[4] = 4.33;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//
//    // periodic
//    y_[1] = y_[0];
//    cubspl = CubicSpline(2,
//                         x_,
//                         y_,
//                         CubicSpline::BoundaryType::periodic,
//                         CubicSpline::BoundaryType::periodic);
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], 2.33, tol);
//    }
//
//    // "not-a-knot" is invalid for n=2
//}
//
//
//TEST_F(CubicSplineTest, ThreePoints)
//{
//    int n_interp = 5;
//    x_interp_[0] = 0.0;
//    x_interp_[1] = 0.1;
//    x_interp_[2] = 0.33;
//    x_interp_[3] = 0.9;
//    x_interp_[4] = 1.0;
//
//    // not-a-knot
//    x_[0] = 0;
//    x_[1] = 0.4;
//    x_[2] = 1.0;
//
//    y_[0] = 1.2;
//    y_[1] = 2.5;
//    y_[2] = 4.8;
//
//    CubicSpline cubspl(3, x_, y_);
//
//    y_ref_[0] = 1.2;
//    y_ref_[1] = 1.5075;
//    y_ref_[2] = 2.259025;
//    y_ref_[3] = 4.3875;
//    y_ref_[4] = 4.8;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//
//    // periodic
//    y_[2] = y_[0];
//    cubspl = CubicSpline(3,
//                         x_,
//                         y_,
//                         CubicSpline::BoundaryType::periodic,
//                         CubicSpline::BoundaryType::periodic);
//
//    y_ref_[0] = 1.2;
//    y_ref_[1] = 1.44375;
//    y_ref_[2] = 2.35383125;
//    y_ref_[3] = 1.2361111111111112;
//    y_ref_[4] = 1.2;
//
//    cubspl.eval(n_interp, x_interp_, y_interp_);
//    for (int i = 0; i != n_interp; ++i)
//    {
//        EXPECT_NEAR(y_interp_[i], y_ref_[i], tol);
//    }
//}


TEST_F(CubicSplineTest, ErrorBound)
{
    // error bound formula used in this test correspond to the complete cubic
    // spline interpolant (exact first_deriv boundary conditions at both ends)

    // tested functions and their derivatives (up to fourth-order)
    // f[i][j] is the j-th derivative of the i-th function
    std::vector<std::vector<std::function<double(double)>>> f = {
        {
            [](double x) { return std::sin(x); },
            [](double x) { return std::cos(x); },
            [](double x) { return -std::sin(x); },
            [](double x) { return -std::cos(x); },
            [](double x) { return std::sin(x); },
        },
        {
            [](double x) { return std::exp(-x); },
            [](double x) { return -std::exp(-x); },
            [](double x) { return std::exp(-x); },
            [](double x) { return -std::exp(-x); },
            [](double x) { return std::exp(-x); },
        },
        {
            [](double x) { return std::log(x); },
            [](double x) { return 1.0 / x; },
            [](double x) { return -1.0 / (x * x); },
            [](double x) { return 2.0 / (x * x * x); },
            [](double x) { return -6.0 / (x * x * x * x); },
        },
        // NOTE: functions with vanishing 4-th derivative should in principle
        // be interpolated exactly by a cubic spline. However, in practice,
        // the presence of floating-point rounding errors would lead to some
        // discrepancy between the interpolant and the original function.
        // Such error is difficult to estimate and the tolerance for comparison
        // is simply set to 1e-12 below.
        {
            [](double x) { return x * x * x; },
            [](double x) { return 3.0 * x * x; },
            [](double x) { return 6.0 * x; },
            [](double  ) { return 6.0; },
            [](double  ) { return 0.0; },
        },
        {
            [](double x) { return x * x; },
            [](double x) { return 2.0 * x; },
            [](double  ) { return 2.0; },
            [](double  ) { return 0.0; },
            [](double  ) { return 0.0; },
        },
        {
            [](double x) { return 2.0 * x; },
            [](double  ) { return 2.0; },
            [](double  ) { return 0.0; },
            [](double  ) { return 0.0; },
            [](double  ) { return 0.0; },
        }
    };

    // knots (logspace)
    int n = 100;
    double xmin = 0.1;
    double xmax = 10;

    double rho0 = std::log(xmin);
    double drho = (std::log(xmax) - rho0) / (n - 1);
    std::for_each(x_, x_ + n, [&](double& x) { x = std::exp(rho0 + (&x - x_) * drho); });

    // places to evaluate the interpolant
    int n_interp = 777;
    double dx_interp = (xmax - xmin) / (n_interp - 1);
    std::for_each(x_interp_, x_interp_ + n_interp,
        [&](double& x) { x = (&x - x_interp_) * dx_interp + xmin; });

    // make sure x_interp is inside the range of x
    x_interp_[0] += tol_;
    x_interp_[n_interp - 1] -= tol_;

    for (size_t i = 0; i < f.size(); ++i)
    {
        std::transform(x_, x_ + n, y_, f[i][0]);

        // complete cubic spline (exact first_deriv boundary conditions at both ends)
        CubicSpline::build(
            n,
            x_,
            y_,
            {BoundaryType::first_deriv, f[i][1](x_[0])},
            {BoundaryType::first_deriv, f[i][1](x_[n - 1])},
            dy_
        );

        CubicSpline::eval(
            n,
            x_,
            y_,
            dy_,
            n_interp,
            x_interp_,
            y_interp_,
            dy_interp_,
            d2y_interp_
        );

        double* y_diff = y_interp_;
        for (int d = 0; d <= 0; ++d)
        {
            std::transform(x_interp_, x_interp_ + n_interp, y_ref_, f[i][d]);
            std::transform(y_diff, y_diff + n_interp, y_ref_, y_diff,
                [](double y, double y_ref) { return std::abs(y - y_ref); });

            double err_bound = error_bound(n, x_, f[i][4], d);
            err_bound = std::max(err_bound, tol_);
            EXPECT_TRUE(std::all_of(y_diff, y_diff + n_interp,
                [err_bound](double diff) { return diff < err_bound; }));

            y_diff += n_interp_max_;
        }
    }
}


TEST_F(CubicSplineTest, Reserve)
{
    int n_spline = 20;
    int n = 1000;
    double x0 = 0.0, dx = 0.01;
    for (int i = 0; i < n; ++i)
    {
        x_[i] = x0 + i * dx;
        y_[i] = std::sin(x_[i]);
    }

    CubicSpline cubspl(n, x0, dx, y_);
    cubspl.reserve(n_spline);
    EXPECT_EQ(cubspl.heap_usage(), n_spline * 2 * n * sizeof(double));

    cubspl = CubicSpline(n, x_, y_);
    cubspl.reserve(n_spline);
    EXPECT_EQ(cubspl.heap_usage(), (1 + n_spline * 2) * n * sizeof(double));
}


TEST_F(CubicSplineTest, CrossCheck)
{
    std::vector<std::string> fnames = {
        "not_a_knot.dat",
    };

    int n = 0, n_interp = 0;
    BoundaryCondition bc_start, bc_end;

    for (auto& fname : fnames)
    {
        read(fname, n, x_, y_, bc_start, bc_end, n_interp, x_interp_, y_ref_, dy_ref_, d2y_ref_);
        CubicSpline cubspl(n, x_, y_, bc_start, bc_end);
        cubspl.eval(n_interp, x_interp_, y_interp_, dy_interp_, d2y_interp_);

        double* interp[] = {y_interp_, dy_interp_, d2y_interp_};
        double* ref[] = {y_ref_, dy_ref_, d2y_ref_};
        for (int i = 0; i < 1; ++i)
        {
            std::transform(interp[i], interp[i] + n_interp, ref[i], interp[i], std::minus<double>());
            EXPECT_TRUE(std::all_of(interp[i], interp[i] + n_interp,
                [](double diff) { return std::abs(diff) < tol_; }));
        }
    }
}


int main()
{
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}


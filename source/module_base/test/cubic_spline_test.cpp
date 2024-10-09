#include "module_base/cubic_spline.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>

#include "gtest/gtest.h"

using ModuleBase::CubicSpline;
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
 *  - xmin, xmax
 *      Returns the first and last knots.
 *
 */
class CubicSplineTest : public ::testing::Test
{
protected:

    CubicSplineTest();
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

    /// y/dy/d2y tolerance for cross-check
    const double tol_[3] = {1e-14, 1e-13, 1e-12};


    /**
     * @brief Sample functions & derivatives in error bound check.
     *
     * @note Functions with vanishing 4-th derivative in an interval
     * should in principle be interpolated exactly by a cubic spline.
     * However, the presence of floating-point rounding errors would
     * lead to some discrepancy between the interpolant and the original
     * function. Such error is not covered by the error bound formula.
     * Functions here should not include those kind of functions.
     *
     */
    std::vector<std::vector<std::function<double(double)>>> f_;

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
};


CubicSplineTest::CubicSplineTest():
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
    d2y_ref_(dy_ref_ + n_interp_max_),
    f_{
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
    }
{}


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

    // read boundary conditions
    std::getline(ifs, line);
    std::stringstream ss(line);
    ss >> bc1 >> bc2;

    auto bc_parse = [](const std::string& bc)
    {
        if (bc == "periodic")
        {
            return BoundaryCondition(BoundaryType::periodic);
        }
        if (bc == "not-a-knot")
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


TEST_F(CubicSplineTest, MultiEval)
{
    int n = 100;
    double xmin = 0.1;
    double xmax = 10;
    double dx = (xmax - xmin) / (n - 1);

    std::for_each(x_, x_ + n, [&](double& x) { x = xmin + (&x - x_) * dx; });

    // empty interpolant with specified knots
    CubicSpline cubspl(n, xmin, dx);
    cubspl.reserve(f_.size());

    for (size_t i = 0; i < f_.size(); ++i)
    {
        std::transform(x_, x_ + n, y_, [this, i](double x) { return f_[i][0](x); });
        cubspl.add(y_, {BoundaryType::first_deriv, f_[i][1](xmin)},
                       {BoundaryType::first_deriv, f_[i][1](xmax)});
    }

    EXPECT_EQ(cubspl.n_spline(), f_.size());

    std::vector<std::vector<double>> err_bound(f_.size(), std::vector<double>(3));
    for (size_t i = 0; i < f_.size(); ++i)
    {
        err_bound[i][0] = error_bound(n, x_, f_[i][4], 0);
        err_bound[i][1] = error_bound(n, x_, f_[i][4], 1);
        err_bound[i][2] = error_bound(n, x_, f_[i][4], 2);
    }

    int n_interp = 1000;
    double dx_interp = (xmax - xmin) / (n_interp - 1);
    std::for_each(x_interp_, x_interp_ + n_interp,
        [&](double& x) { x = (&x - x_interp_) * dx_interp + xmin; });

    for (int p = 0; p < n_interp; ++p)
    {
        cubspl.multi_eval(x_interp_[p], y_interp_, dy_interp_, d2y_interp_);
        double ytmp, dytmp, d2ytmp;
        for (size_t i = 0; i < f_.size(); ++i)
        {
            EXPECT_LT(std::abs(y_interp_[i] - f_[i][0](x_interp_[p])), err_bound[i][0]);
            EXPECT_LT(std::abs(dy_interp_[i] - f_[i][1](x_interp_[p])), err_bound[i][1]);
            EXPECT_LT(std::abs(d2y_interp_[i] - f_[i][2](x_interp_[p])), err_bound[i][2]);

            cubspl.eval(1, &x_interp_[p], &ytmp, &dytmp, &d2ytmp, i);
            EXPECT_NEAR(ytmp, y_interp_[i], tol_[0]);
            EXPECT_NEAR(dytmp, dy_interp_[i], tol_[1]);
            EXPECT_NEAR(d2ytmp, d2y_interp_[i], tol_[2]);
        }
    }
}


TEST_F(CubicSplineTest, ErrorBound)
{
    // Error bound formula used in this test correspond to the complete cubic
    // spline interpolant (exact first_deriv boundary conditions at both ends).
    // This does not apply to other types of boundary conditions.

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
    x_interp_[0] += tol_[0];
    x_interp_[n_interp - 1] -= tol_[0];

    for (size_t i = 0; i < f_.size(); ++i)
    {
        std::transform(x_, x_ + n, y_, f_[i][0]);

        // complete cubic spline (exact first_deriv boundary conditions at both ends)
        CubicSpline::build(
            n, x_, y_,
            {BoundaryType::first_deriv, f_[i][1](x_[0])},
            {BoundaryType::first_deriv, f_[i][1](x_[n - 1])},
            dy_
        );

        CubicSpline::eval(
            n, x_, y_, dy_,
            n_interp, x_interp_, y_interp_, dy_interp_, d2y_interp_
        );

        double* diff[3] = {y_interp_, dy_interp_, d2y_interp_};
        for (int d = 0; d < 3; ++d)
        {
            std::transform(x_interp_, x_interp_ + n_interp, diff[d], diff[d],
                [&](double x, double y) { return std::abs(y - f_[i][d](x)); });

            double err_bound = error_bound(n, x_, f_[i][4], d);
            EXPECT_TRUE(std::all_of(diff[d], diff[d] + n_interp,
                [err_bound](double diff) { return diff < err_bound; }));
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


TEST_F(CubicSplineTest, MinMax)
{
    int n = 1000;
    double x0 = 0.0, dx = 0.01;
    for (int i = 0; i < n; ++i)
    {
        x_[i] = x0 + i * dx;
        y_[i] = std::sin(x_[i]);
    }

    CubicSpline cubspl(n, x_, y_);
    EXPECT_EQ(cubspl.xmin(), x_[0]);
    EXPECT_EQ(cubspl.xmax(), x_[n - 1]);

    int m = 300;
    cubspl = CubicSpline(m, x0, dx, y_);
    EXPECT_EQ(cubspl.xmin(), x0);
    EXPECT_EQ(cubspl.xmax(), x0 + (m - 1) * dx);

}


TEST_F(CubicSplineTest, CrossCheck)
{
    std::vector<std::string> fnames = {
        "./data/sin_not_a_knot.dat",
        "./data/cos_periodic.dat",
        "./data/exp_first_deriv.dat",
        "./data/log_second_deriv.dat",
        "./data/sqrt_mix_bc.dat",
        "./data/two_points_periodic.dat",
        "./data/two_points_first_deriv.dat",
        "./data/two_points_second_deriv.dat",
        "./data/three_points_not_a_knot.dat",
    };

    int n = 0, n_interp = 0;
    BoundaryCondition bc_start, bc_end;

    for (const auto& fname : fnames)
    {
        read(fname, n, x_, y_, bc_start, bc_end,
             n_interp, x_interp_, y_ref_, dy_ref_, d2y_ref_);
        CubicSpline cubspl(n, x_, y_, bc_start, bc_end);
        cubspl.eval(n_interp, x_interp_, y_interp_, dy_interp_, d2y_interp_);

        double* diff[] = {y_interp_, dy_interp_, d2y_interp_};
        double* ref[] = {y_ref_, dy_ref_, d2y_ref_};
        for (int d = 0; d < 3; ++d)
        {
            std::transform(diff[d], diff[d] + n_interp, ref[d], diff[d], std::minus<double>());
            EXPECT_TRUE(std::all_of(diff[d], diff[d] + n_interp,
                [this, d](double diff) { return std::abs(diff) < tol_[d]; }));
        }
    }
}


int main()
{
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}



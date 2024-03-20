#include "cubic_spline.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>

namespace ModuleBase
{

CubicSpline::CubicSpline(CubicSpline const& other)
{
    n_ = other.n_;
    is_uniform_ = other.is_uniform_;
    uniform_thr_ = other.uniform_thr_;

    x_ = new double[n_];
    y_ = new double[n_];
    s_ = new double[n_];

    std::memcpy(x_, other.x_, n_ * sizeof(double));
    std::memcpy(y_, other.y_, n_ * sizeof(double));
    std::memcpy(s_, other.s_, n_ * sizeof(double));
}

CubicSpline& CubicSpline::operator=(CubicSpline const& other)
{
    if (this != &other)
    {
        cleanup();
        n_ = other.n_;
        is_uniform_ = other.is_uniform_;
        uniform_thr_ = other.uniform_thr_;

        x_ = new double[n_];
        y_ = new double[n_];
        s_ = new double[n_];

        std::memcpy(x_, other.x_, n_ * sizeof(double));
        std::memcpy(y_, other.y_, n_ * sizeof(double));
        std::memcpy(s_, other.s_, n_ * sizeof(double));
    }

    return *this;
}

void CubicSpline::cleanup()
{
    delete[] x_;
    delete[] y_;
    delete[] s_;

    x_ = nullptr;
    y_ = nullptr;
    s_ = nullptr;
}

void CubicSpline::check_build(const int n,
                              const double* const x,
                              const double* const y,
                              BoundaryCondition bc_start,
                              BoundaryCondition bc_end)
{
    assert(n > 1);

    // periodic boundary condition must apply to both ends
    assert((bc_start == BoundaryCondition::periodic) == (bc_end == BoundaryCondition::periodic));

    // y[0] must equal y[n-1] for periodic boundary condition
    assert(bc_start != BoundaryCondition::periodic || y[0] == y[n - 1]);

    // if one of the boundary condition is not-a-knot, n must be larger than 2
    assert((bc_start != BoundaryCondition::not_a_knot && bc_end != BoundaryCondition::not_a_knot) || n > 2);

    for (int i = 0; i != n - 1; ++i)
    {
        assert(x[i + 1] > x[i] && "CubicSpline: x must be strictly increasing");
    }
}

void CubicSpline::check_interp(const int n,
                               const double* const x,
                               const double* const y,
                               const double* const s,
                               const int n_interp,
                               const double* const x_interp,
                               double* const y_interp,
                               double* const dy_interp)
{
    assert(n > 1 && x && y && s);     // make sure the interpolant exists
    assert(n_interp > 0 && x_interp); // make sure the interpolation points exist
    assert(y_interp || dy_interp);    // make sure at least one of y or dy is not null

    // check that x_interp is in the range of the interpolant
    assert(std::all_of(x_interp, x_interp + n_interp,
                [n, x](const double xi) { return xi >= x[0] && xi <= x[n - 1]; }));
}

void CubicSpline::build(const int n,
                        const double* const x,
                        const double* const y,
                        double* const s,
                        BoundaryCondition bc_start,
                        BoundaryCondition bc_end,
                        const double deriv_start,
                        const double deriv_end)
{
    check_build(n, x, y, bc_start, bc_end);

    if (n == 2 && bc_start == BoundaryCondition::periodic)
    { // in this case the polynomial is a constant
        s[0] = s[1] = 0.0;
    }
    else if (n == 3 && bc_start == BoundaryCondition::not_a_knot && bc_end == BoundaryCondition::not_a_knot)
    { // in this case two conditions coincide; simply build a parabola that passes through the three data points
        double idx10 = 1. / (x[1] - x[0]);
        double idx21 = 1. / (x[2] - x[1]);
        double idx20 = 1. / (x[2] - x[0]);

        s[0] = -y[0] * (idx10 + idx20) + y[1] * (idx21 + idx10) + y[2] * (idx20 - idx21);
        s[1] = -y[1] * (-idx10 + idx21) + y[0] * (idx20 - idx10) + y[2] * (idx21 - idx20);
        s[2] = s[1] + 2.0 * (-y[1] * idx10 + y[2] * idx20) + 2.0 * y[0] * idx10 * idx20 * (x[2] - x[1]);
    }
    else
    {
        double* dx = new double[n - 1];
        double* dd = new double[n - 1];

        // tridiagonal linear system (cyclic if using the periodic boundary condition)
        double* diag = new double[n];
        double* subdiag = new double[n - 1];
        double* supdiag = new double[n - 1];

        for (int i = 0; i != n - 1; ++i)
        {
            dx[i] = x[i + 1] - x[i];
            dd[i] = (y[i + 1] - y[i]) / dx[i];
        }

        // common part of the tridiagonal linear system
        for (int i = 1; i != n - 1; ++i)
        {
            diag[i] = 2.0 * (dx[i - 1] + dx[i]);
            supdiag[i] = dx[i - 1];
            subdiag[i - 1] = dx[i];
            s[i] = 3.0 * (dd[i - 1] * dx[i] + dd[i] * dx[i - 1]);
        }

        if (bc_start == BoundaryCondition::periodic)
        {
            // exclude s[n-1] and solve a a cyclic tridiagonal linear system of size n-1
            diag[0] = 2.0 * (dx[n - 2] + dx[0]);
            supdiag[0] = dx[n - 2];
            subdiag[n - 2] = dx[0];
            s[0] = 3.0 * (dd[0] * dx[n - 2] + dd[n - 2] * dx[0]);
            solve_cyctri(n - 1, diag, supdiag, subdiag, s);
            s[n - 1] = s[0];
        }
        else
        {
            switch (bc_start)
            {
            case BoundaryCondition::first_deriv:
                diag[0] = 1.0 * dx[0];
                supdiag[0] = 0.0;
                s[0] = deriv_start * dx[0];
                break;
            case BoundaryCondition::second_deriv:
                diag[0] = 2.0 * dx[0];
                supdiag[0] = 1.0 * dx[0];
                s[0] = (3.0 * dd[0] - 0.5 * deriv_start * dx[0]) * dx[0];
                break;
            default: // BoundaryCondition::not_a_knot
                diag[0] = dx[1];
                supdiag[0] = x[2] - x[0];
                s[0] = (dd[0] * dx[1] * (dx[0] + 2 * (x[2] - x[0])) + dd[1] * dx[0] * dx[0]) / (x[2] - x[0]);
            }

            switch (bc_end)
            {
            case BoundaryCondition::first_deriv:
                diag[n - 1] = 1.0 * dx[n - 2];
                subdiag[n - 2] = 0.0;
                s[n - 1] = deriv_end * dx[n - 2];
                break;
            case BoundaryCondition::second_deriv:
                diag[n - 1] = 2.0 * dx[n - 2];
                subdiag[n - 2] = 1.0 * dx[n - 2];
                s[n - 1] = (3.0 * dd[n - 2] + 0.5 * deriv_end * dx[n - 2]) * dx[n - 2];
                break;
            default: // BoundaryCondition::not_a_knot
                diag[n - 1] = dx[n - 3];
                subdiag[n - 2] = x[n - 1] - x[n - 3];
                s[n - 1] = (dd[n - 2] * dx[n - 3] * (dx[n - 2] + 2 * (x[n - 1] - x[n - 3]))
                            + dd[n - 3] * dx[n - 2] * dx[n - 2])
                           / (x[n - 1] - x[n - 3]);
            }

            int NRHS = 1;
            int LDB = n;
            int INFO = 0;
            int N = n;

            dgtsv_(&N, &NRHS, subdiag, diag, supdiag, s, &LDB, &INFO);
        }

        delete[] diag;
        delete[] subdiag;
        delete[] supdiag;
        delete[] dx;
        delete[] dd;
    }
}

void CubicSpline::build(const int n,
                        const double* const x,
                        const double* const y,
                        BoundaryCondition bc_start,
                        BoundaryCondition bc_end,
                        const double deriv_start,
                        const double deriv_end)
{
    cleanup();

    n_ = n;
    x_ = new double[n];
    y_ = new double[n];
    s_ = new double[n];

    std::memcpy(x_, x, sizeof(double) * n);
    std::memcpy(y_, y, sizeof(double) * n);
    build(n_, x_, y_, s_, bc_start, bc_end, deriv_start, deriv_end);

    double dx_avg = (x[n - 1] - x[0]) / (n - 1);
    is_uniform_ = std::all_of(x, x + n,
            [dx_avg, &x](const double& xi) -> bool { return std::abs(xi - (&xi - x) * dx_avg) < 1e-15; });
}

void CubicSpline::eval(const int n,
                       const double* const x,
                       const double* const y,
                       const double* const s,
                       const int n_interp,
                       const double* const x_interp,
                       double* const y_interp,
                       double* const dy_interp)
{
    check_interp(n, x, y, s, n_interp, x_interp, y_interp, dy_interp);
    _eval(x, y, s, n_interp, x_interp, y_interp, dy_interp, _gen_search(n, x));
}

void CubicSpline::eval(const int n_interp,
                       const double* const x_interp,
                       double* const y_interp,
                       double* const dy_interp)
{
    check_interp(n_, x_, y_, s_, n_interp, x_interp, y_interp, dy_interp);
    _eval(x_, y_, s_, n_interp, x_interp, y_interp, dy_interp, _gen_search(n_, x_, is_uniform_));
}

std::function<int(double)> CubicSpline::_gen_search(const int n, const double* const x, int is_uniform)
{
    if (is_uniform != 0 && is_uniform != 1)
    {
        double dx_avg = (x[n - 1] - x[0]) / (n - 1);
        is_uniform = std::all_of(x, x + n,
                [dx_avg, &x] (const double& xi) { return std::abs(xi - (&xi - x) * dx_avg) < 1e-15; });
    }

    if (is_uniform)
    {
        double dx = x[1] - x[0];
        return [dx, n, x](double xi) -> int { return xi == x[n - 1] ? n - 2 : xi / dx; };
    }
    else
    {
        return [n, x](double xi) -> int { return (std::upper_bound(x, x + n, xi) - x) - 1; };
    }
}

void CubicSpline::_eval(const double* const x,
                        const double* const y,
                        const double* const s,
                        const int n_interp,
                        const double* const x_interp,
                        double* const y_interp,
                        double* const dy_interp,
                        std::function<int(double)> search)
{
    void (*poly_eval)(double, double, double, double, double, double*, double*) =
        y_interp ? (dy_interp ? _poly_eval<true, true>: _poly_eval<true, false>) : _poly_eval<false, true>;

    for (int i = 0; i != n_interp; ++i)
    {
        int p = search(x_interp[i]);
        double w = x_interp[i] - x[p];

        double dx = x[p + 1] - x[p];
        double dd = (y[p + 1] - y[p]) / dx;

        double c0 = y[p];
        double c1 = s[p];
        double c3 = (s[p] + s[p + 1] - 2.0 * dd) / (dx * dx);
        double c2 = (dd - s[p]) / dx - c3 * dx;

        poly_eval(w, c0, c1, c2, c3, y_interp + i, dy_interp + i);
    }
}

void CubicSpline::solve_cyctri(const int n,
                               double* const diag,
                               double* const supdiag,
                               double* const subdiag,
                               double* const b)
{
    // This function uses the Sherman-Morrison formula to convert a cyclic tridiagonal linear system
    // into a tridiagonal one.

    // NOTE all diags will be overwritten in this function!

    // flexible non-zero parameters that can affect the condition number of the tridiagonal linear system below
    // may have some smart choice, set to 1 for now
    double alpha = 1.0;
    double beta = 1.0;

    double* bp = new double[2 * n];
    std::memcpy(bp, b, sizeof(double) * n);
    bp[n] = 1. / alpha;
    bp[2 * n - 1] = 1. / beta;
    for (int i = n + 1; i != 2 * n - 1; ++i)
    {
        bp[i] = 0.0;
    }

    diag[0] -= supdiag[n - 1] * beta / alpha;
    diag[n - 1] -= subdiag[n - 1] * alpha / beta;

    int nrhs = 2;
    int info = 0;
    int N = n;
    int ldb = n;
    dgtsv_(&N, &nrhs, subdiag, diag, supdiag, bp, &ldb, &info);

    double fac = (beta * supdiag[n - 1] * bp[0] + alpha * subdiag[n - 1] * bp[n - 1])
                 / (1. + beta * supdiag[n - 1] * bp[n] + alpha * subdiag[n - 1] * bp[2 * n - 1]);

    std::memcpy(b, bp, sizeof(double) * n);
    for (int i = 0; i != n; ++i)
    {
        b[i] -= fac * bp[n + i];
    }

    delete[] bp;
}

} // namespace ModuleBase

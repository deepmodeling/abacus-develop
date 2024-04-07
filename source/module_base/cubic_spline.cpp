#include "cubic_spline.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>
#include <chrono>

namespace ModuleBase
{
    std::chrono::duration<double> CubicSpline::time_c01;
    std::chrono::duration<double> CubicSpline::time_c23;
    std::chrono::duration<double> CubicSpline::time_poly;
    std::chrono::duration<double> CubicSpline::time_seg;
    std::chrono::duration<double> CubicSpline::time_dd;
    std::chrono::duration<double> CubicSpline::time_alloc;

CubicSpline::CubicSpline(
    const int n,
    const double* const x,
    const double* const y,
    BoundaryCondition bc_start,
    BoundaryCondition bc_end,
    const double deriv_start,
    const double deriv_end
)
{
    x_ = std::vector<double>(x, x + n);
    y_ = std::vector<double>(y, y + n);
    s_.resize(n);
    build(n, x, y, s_.data(), -1, bc_start, bc_end, deriv_start, deriv_end);
}

CubicSpline::CubicSpline(
    const int n,
    const double x0,
    const double dx,
    const double* const y,
    BoundaryCondition bc_start,
    BoundaryCondition bc_end,
    const double deriv_start,
    const double deriv_end
)
{
    x0_ = x0;
    dx_ = dx;
    y_ = std::vector<double>(y, y + n);
    s_.resize(n);
    build(n, nullptr, y, s_.data(), dx, bc_start, bc_end, deriv_start, deriv_end);
}

void CubicSpline::_validate_build(
    const int n,
    const double* const x,
    const double* const y,
    const double dx,
    BoundaryCondition bc_start,
    BoundaryCondition bc_end
)
{
//#ifdef __DEBUG
    assert(n > 1);

    // if periodic boundary condition is specified, it must be applied to both ends
    assert((bc_start == BoundaryCondition::periodic) == (bc_end == BoundaryCondition::periodic));

    // y[0] must equal y[n-1] for periodic boundary condition
    assert(bc_start != BoundaryCondition::periodic || y[0] == y[n - 1]);

    // not-a-knot boundary condition requires the existence of "internal" knot, n must be at least 3
    assert((bc_start != BoundaryCondition::not_a_knot && bc_end != BoundaryCondition::not_a_knot) || n > 2);

    // knots must be STRICTLY increasing
    assert((x && std::is_sorted(x, x + n, std::less_equal<double>())) || dx > 0.0);
//#endif
}

void CubicSpline::_validate_eval(
    const int n,
    const double* const x,
    const double* const y,
    const double* const s,
    const int n_interp,
    const double* const x_interp,
    const double x0,
    const double dx
)
{
//#ifdef __DEBUG
    assert(n > 1 && y && s);
    assert((x && std::is_sorted(x, x + n, std::less_equal<double>())) || dx > 0.0);

    assert(n_interp > 0 && x_interp);

    double xmin = x ? x[0] : x0;
    double xmax = x ? x[n - 1] : x0 + (n - 1) * dx;
    assert(std::all_of(x_interp, x_interp + n_interp,
                       [xmin, xmax](const double x_i) { return xmin <= x_i && x_i <= xmax; }));
//#endif
}

void CubicSpline::build(
    const int n,
    const double* const x,
    const double* const y,
    double* const s,
    const double dx_uniform,
    BoundaryCondition bc_start,
    BoundaryCondition bc_end,
    const double deriv_start,
    const double deriv_end
)
{
    _validate_build(n, x, y, dx_uniform, bc_start, bc_end);

    if (n == 2 && bc_start == BoundaryCondition::periodic)
    { // in this case the polynomial is a constant
        s[0] = s[1] = 0.0;
    }
    else if (n == 3 && bc_start == BoundaryCondition::not_a_knot && bc_end == BoundaryCondition::not_a_knot)
    { // in this case two conditions coincide; simply build a parabola that passes through the three data points
        double idx10, idx21, idx20, dx21;
        if (x)
        {
            idx10 = 1. / (x[1] - x[0]);
            idx21 = 1. / (x[2] - x[1]);
            idx20 = 1. / (x[2] - x[0]);
            dx21 = x[2] - x[1];
        }
        else
        {
            idx10 = idx21 = 1. / dx_uniform;
            idx20 = 0.5 / dx_uniform;
            dx21 = dx_uniform;
        }

        s[0] = -y[0] * (idx10 + idx20) + y[1] * (idx21 + idx10) + y[2] * (idx20 - idx21);
        s[1] = -y[1] * (-idx10 + idx21) + y[0] * (idx20 - idx10) + y[2] * (idx21 - idx20);
        s[2] = s[1] + 2.0 * (-y[1] * idx10 + y[2] * idx20) + 2.0 * y[0] * idx10 * idx20 * dx21;
    }
    else
    {
        double* buffer = new double[5 * n];

        double* dx = buffer;
        double* dd = buffer + n; // divided differences

        std::adjacent_difference(y, y + n, dd);
        dd += 1; // the first element computed by adjacent_difference is not a difference

        if (x)
        {
            std::adjacent_difference(x, x + n, dx);
            dx += 1;
        }
        else
        {
            std::fill(dx, dx + n - 1, dx_uniform);
        }

        std::transform(dd, dd + n - 1, dx, dd, std::divides<double>());

        // tridiagonal linear system (cyclic tridiagonal if periodic boundary condition)
        double* diag = buffer + 2 * n;
        double* subdiag = buffer + 3 * n;
        double* supdiag = buffer + 4 * n;

        // common part of the tridiagonal linear system

        std::memcpy(subdiag, dx + 1, sizeof(double) * (n - 2));
        std::memcpy(supdiag + 1, dx, sizeof(double) * (n - 2));

        for (int i = 1; i != n - 1; ++i)
        {
            diag[i] = 2.0 * (dx[i - 1] + dx[i]);
            s[i] = 3.0 * (dd[i - 1] * dx[i] + dd[i] * dx[i - 1]);
        }

        // below for the boundary-condition-specific part

        if (bc_start == BoundaryCondition::periodic)
        {
            // exclude s[n-1] and solve a a cyclic tridiagonal linear system of size n-1
            diag[0] = 2.0 * (dx[n - 2] + dx[0]);
            supdiag[0] = dx[n - 2];
            subdiag[n - 2] = dx[0];
            s[0] = 3.0 * (dd[0] * dx[n - 2] + dd[n - 2] * dx[0]);
            _solve_cyctri(n - 1, diag, supdiag, subdiag, s);
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
                supdiag[0] = dx[0] + dx[1];
                s[0] = (dd[0] * dx[1] * (dx[0] + 2 * (dx[0] + dx[1])) + dd[1] * dx[0] * dx[0]) / (dx[0] + dx[1]);
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
                subdiag[n - 2] = dx[n - 3] + dx[n - 2];
                s[n - 1] = (dd[n - 2] * dx[n - 3] * (dx[n - 2] + 2 * (dx[n - 3] + dx[n - 2]))
                            + dd[n - 3] * dx[n - 2] * dx[n - 2])
                           / (dx[n - 3] + dx[n - 2]);
            }

            int NRHS = 1;
            int LDB = n;
            int INFO = 0;
            int N = n;

            dgtsv_(&N, &NRHS, subdiag, diag, supdiag, s, &LDB, &INFO);
        }

        delete[] buffer;
    }
}

void CubicSpline::eval(
    const int n,
    const double* const x,
    const double* const y,
    const double* const s,
    const int n_interp,
    const double* const x_interp,
    double* const y_interp,
    double* const s_interp,
    const double x0,
    const double dx
)
{
    _validate_eval(n, x, y, s, n_interp, x_interp, x0, dx);

    iclock::time_point start = iclock::now();

    // segind[p] will store the index of the polynomial segment that contains x_interp[p]
    int* segind = new int[n_interp];

    double* buffer = new double[6 * n_interp];
    double* dd = buffer; // divided difference
    double* w  = buffer +     n_interp;
    double* c0 = buffer + 2 * n_interp;
    double* c1 = buffer + 3 * n_interp;
    double* c2 = buffer + 4 * n_interp;
    double* c3 = buffer + 5 * n_interp;

    time_alloc += iclock::now() - start;
    start = iclock::now();

    if (x)
    {
        std::transform(x_interp, x_interp + n_interp, segind,
                       [n, x](double x_i) { return _index(n, x, x_i); });

        std::transform(segind, segind + n_interp, x_interp, w,
                       [x](int p, double x_i) { return x_i - x[p]; });

        std::transform(segind, segind + n_interp, dd,
                       [x, y](int p) { return (y[p + 1] - y[p]) / (x[p + 1] - x[p]); });

        std::transform(segind, segind + n_interp, dd, c2,
                       [x, s](int p, double d) { return (3.0 * d - 2.0 * s[p] - s[p + 1]) / (x[p + 1] - x[p]); });

        std::transform(segind, segind + n_interp, dd, c3,
                       [x, s](int p, double d) { return (s[p] + s[p + 1] - 2.0 * d) / std::pow(x[p + 1] - x[p], 2); });
    }
    else
    {
        std::transform(x_interp, x_interp + n_interp, segind,
                       [n, x0, dx](double x_i) { return _index(n, x0, dx, x_i); });

        std::transform(segind, segind + n_interp, x_interp, w,
                       [x0, dx](int p, double x_i) { return x_i - x0 - p * dx; });
        time_seg += iclock::now() - start;
        start = iclock::now();

        std::transform(segind, segind + n_interp, dd,
                       [dx, y](int p) { return (y[p + 1] - y[p]) / dx; });
        time_dd += iclock::now() - start;
        start = iclock::now();

        std::transform(segind, segind + n_interp, dd, c2,
                       [dx, s](int p, double d) { return (3.0 * d - 2.0 * s[p] - s[p + 1]) / dx; });

        std::transform(segind, segind + n_interp, dd, c3,
                       [dx, s](int p, double d) { return (s[p] + s[p + 1] - 2.0 * d) / (dx * dx); });
        time_c23 += iclock::now() - start;
        start = iclock::now();
    }

    std::transform(segind, segind + n_interp, c0, [y](int p) { return y[p]; });
    std::transform(segind, segind + n_interp, c1, [s](int p) { return s[p]; });

    time_c01 += iclock::now() - start;
    start = iclock::now();

    if (y_interp)
    {
        for (int i = 0; i < n_interp; ++i)
        {
            y_interp[i] = ((c3[i] * w[i] + c2[i]) * w[i] + c1[i]) * w[i] + c0[i];
        }
    }

    if (s_interp)
    {
        for (int i = 0; i < n_interp; ++i)
        {
            s_interp[i] = (3.0 * c3[i] * w[i] + 2.0 * c2[i]) * w[i] + c1[i];
        }
    }
    time_poly += iclock::now() - start;

    delete[] segind;
    delete[] buffer;
}

void CubicSpline::eval(
    const int n_interp,
    const double* const x_interp,
    double* const y_interp,
    double* const s_interp
) const
{
    _validate_eval(y_.size(), x_.data(), y_.data(), s_.data(), n_interp, x_interp, x0_, dx_);
    eval(y_.size(), x_.data(), y_.data(), s_.data(), n_interp, x_interp, y_interp, s_interp, x0_, dx_);
}

void CubicSpline::_solve_cyctri(
    const int n,
    double* const diag,
    double* const supdiag,
    double* const subdiag,
    double* const b
)
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

#include "cubic_spline.h"

#include <cassert>
#include <cmath>
#include <cstring>

namespace ModuleBase
{

CubicSpline::~CubicSpline()
{
    cleanup();
}

void CubicSpline::cleanup()
{
    delete[] x_;
    delete[] c0_;
    delete[] c1_;
    delete[] c2_;
    delete[] c3_;

    x_ = nullptr;
    c0_ = nullptr;
    c1_ = nullptr;
    c2_ = nullptr;
    c3_ = nullptr;
}

void CubicSpline::sanity_check(const int n,
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

void CubicSpline::build(const int n,
                        const double* const x,
                        const double* const y,
                        BoundaryCondition bc_start,
                        BoundaryCondition bc_end,
                        const double deriv_start,
                        const double deriv_end)
{
    sanity_check(n, x, y, bc_start, bc_end);

    cleanup();

    n_ = n;
    x_ = new double[n];
    std::memcpy(x_, x, sizeof(double) * n);

    c0_ = new double[n_ - 1];
    c1_ = new double[n_ - 1];
    c2_ = new double[n_ - 1];
    c3_ = new double[n_ - 1];

    if (n == 2 && bc_start == BoundaryCondition::periodic)
    {
        // in this case the polynomial is a constant
        c0_[0] = y[0];
        c1_[0] = c2_[0] = c3_[0] = 0.0;
    }
    else if (n == 3 && bc_start == BoundaryCondition::not_a_knot && bc_end == BoundaryCondition::not_a_knot)
    {
        // in this case two conditions coincide; simply build a parabola that passes through the three data points
        double idx10 = 1. / (x[1] - x[0]);
        double idx21 = 1. / (x[2] - x[1]);
        double idx20 = 1. / (x[2] - x[0]);

        c0_[0] = y[0];
        c1_[0] = -y[0] * (idx10 + idx20) + y[1] * (idx21 + idx10) + y[2] * (idx20 - idx21);
        c2_[0] = y[0] * idx10 * idx20 - y[1] * idx10 * idx21 + y[2] * idx20 * idx21;
        c3_[0] = 0.0;

        c0_[1] = y[1];
        c1_[1] = -y[1] * (-idx10 + idx21) + y[0] * (idx20 - idx10) + y[2] * (idx21 - idx20);
        c2_[1] = -y[1] * idx10 * idx21 + y[0] * idx10 * idx20 + y[2] * idx21 * idx20;
        c3_[1] = 0.0;
    }
    else
    {
        double* dx = new double[n - 1];
        double* dd = new double[n - 1];

        // tridiagonal linear system (cyclic if using the periodic boundary condition)
        double* diag = new double[n];
        double* subdiag = new double[n - 1];
        double* supdiag = new double[n - 1];
        double* s = new double[n];

        is_uniform_ = true;
        double dx_avg = (x[n - 1] - x[0]) / (n - 1);

        for (int i = 0; i != n - 1; ++i)
        {
            dx[i] = x[i + 1] - x[i];
            dd[i] = (y[i + 1] - y[i]) / dx[i];

            if (std::abs(dx[i] - dx_avg) > uniform_thr_)
            {
                is_uniform_ = false;
            }
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
            ;
            solve_cyctri(n - 1, diag, supdiag, subdiag, s);
            s[n - 1] = s[0];
        }
        else
        {

            switch (bc_start)
            {
            case BoundaryCondition::not_a_knot:
                diag[0] = dx[1];
                supdiag[0] = x[2] - x[0];
                s[0] = (dd[0] * dx[1] * (dx[0] + 2 * (x[2] - x[0])) + dd[1] * dx[0] * dx[0]) / (x[2] - x[0]);
                break;
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
            default:
                assert(false && "CubicSpline: BoundaryCondition error!"); // not supposed to happend
            }

            switch (bc_end)
            {
            case BoundaryCondition::not_a_knot:
                diag[n - 1] = dx[n - 3];
                subdiag[n - 2] = x[n - 1] - x[n - 3];
                s[n - 1] = (dd[n - 2] * dx[n - 3] * (dx[n - 2] + 2 * (x[n - 1] - x[n - 3]))
                            + dd[n - 3] * dx[n - 2] * dx[n - 2])
                           / (x[n - 1] - x[n - 3]);
                break;
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
            default:
                assert(false && "CubicSpline: BoundaryCondition error!"); // not supposed to happend
            }

            int NRHS = 1;
            int LDB = n;
            int INFO = 0;
            int N = n;

            dgtsv_(&N, &NRHS, subdiag, diag, supdiag, s, &LDB, &INFO);
        }

        std::memcpy(c0_, y, sizeof(double) * (n - 1));
        std::memcpy(c1_, s, sizeof(double) * (n - 1));
        for (int i = 0; i != n - 1; ++i)
        {
            c3_[i] = (s[i] + s[i + 1] - 2 * dd[i]) / (dx[i] * dx[i]);
            c2_[i] = (dd[i] - s[i]) / dx[i] - c3_[i] * dx[i];
        }

        delete[] diag;
        delete[] subdiag;
        delete[] supdiag;
        delete[] dx;
        delete[] dd;
        delete[] s;
    }
}

void CubicSpline::get(const int n, const double* const x, double* const y)
{
    assert(x_ && c0_ && c1_ && c2_ && c3_); // make sure the interpolant exists

    int p = 0;
    double w = 0.0;

    for (int i = 0; i != n; ++i)
    {
        assert(x[i] >= x_[0] && x[i] <= x_[n_ - 1]);
    }

    if (is_uniform_)
    {
        double dx = x_[1] - x_[0];
        for (int i = 0; i != n; ++i)
        {
            p = x[i] == x_[n_ - 1] ? n_ - 2 : x[i] / dx;
            w = x[i] - x_[p];
            y[i] = ((c3_[p] * w + c2_[p]) * w + c1_[p]) * w + c0_[p];
        }
    }
    else
    {
        for (int i = 0; i != n; ++i)
        {
            int q = n_ - 1;
            while (q - p > 1)
            {
                int m = (p + q) / 2;
                if (x_[m] > x[i])
                {
                    q = m;
                }
                else
                {
                    p = m;
                }
            }
            w = x[i] - x_[p];
            y[i] = ((c3_[p] * w + c2_[p]) * w + c1_[p]) * w + c0_[p];
        }
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

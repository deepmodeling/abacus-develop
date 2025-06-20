#include "cubic_spline.h"

#include <cassert>
#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

using ModuleBase::CubicSpline;

extern "C"
{
    // solve a tridiagonal linear system
    void dgtsv_(int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO);
};


CubicSpline::BoundaryCondition::BoundaryCondition(BoundaryType type)
    : type(type)
{
    assert(type == BoundaryType::periodic || type == BoundaryType::not_a_knot);
}


CubicSpline::BoundaryCondition::BoundaryCondition(BoundaryType type, double val)
    : type(type), val(val)
{
    assert(type == BoundaryType::first_deriv || type == BoundaryType::second_deriv);
}


CubicSpline::CubicSpline(
    int n,
    const double* x,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end
): n_spline_(1), n_(n), xmin_(x[0]), xmax_(x[n - 1]), x_(x, x + n), y_(2 * n)
{
    std::copy(y, y + n, y_.begin());
    build(n, x, y, bc_start, bc_end, &y_[n]);
}


CubicSpline::CubicSpline(
    int n,
    double x0,
    double dx,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end
): n_spline_(1), n_(n), xmin_(x0), xmax_(x0 + (n - 1) * dx), dx_(dx), y_(2 * n)
{
    std::copy(y, y + n, y_.begin());
    build(n, dx, y, bc_start, bc_end, &y_[n]);
}


CubicSpline::CubicSpline(int n, const double* x)
    : n_spline_(0), n_(n), xmin_(x[0]), xmax_(x[n - 1]), x_(x, x + n)
{
}


CubicSpline::CubicSpline(int n, double x0, double dx)
    : n_spline_(0), n_(n), xmin_(x0), xmax_(x0 + (n - 1) * dx), dx_(dx)
{
}


void CubicSpline::add(
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end
)
{
    int offset = n_spline_ * 2 * n_;
    y_.resize(offset + 2 * n_);

    std::copy(y, y + n_, &y_[offset]);
    double* dy = &y_[offset + n_];
    if (x_.empty()) // evenly spaced knots
    {
        build(n_, dx_, y, bc_start, bc_end, dy);
    }
    else
    {
        build(n_, x_.data(), y, bc_start, bc_end, dy);
    }
    ++n_spline_;
}


void CubicSpline::eval(
    int n_interp,
    const double* x_interp,
    double* y_interp,
    double* dy_interp,
    double* d2y_interp,
    int i_spline
) const
{
    assert(0 <= i_spline && i_spline < n_spline_);

    const double* y = &y_[i_spline * 2 * n_];
    const double* dy = y + n_;
    if (x_.empty()) // evenly spaced knots
    {
        eval(n_, xmin_, dx_, y, dy, n_interp, x_interp, y_interp, dy_interp, d2y_interp);
    }
    else
    {
        eval(n_, x_.data(), y, dy, n_interp, x_interp, y_interp, dy_interp, d2y_interp);
    }
}


void CubicSpline::multi_eval(
    int n_spline,
    const int* i_spline,
    double x_interp,
    double* y_interp,
    double* dy_interp,
    double* d2y_interp
) const
{
    assert(std::all_of(i_spline, i_spline + n_spline,
        [this](int i) { return 0 <= i && i < n_spline_; }));
    _validate_eval(n_, {xmin_, dx_}, x_.empty() ? nullptr : x_.data(),
                   y_.data(), &y_[n_], 1, &x_interp);

    int p = 0;
    double dx = 0.0, r = 0.0; 
    if (x_.empty()) // evenly spaced knots
    {
        p = _index(n_, xmin_, dx_, x_interp);
        dx = dx_;
        r = (x_interp - xmin_) / dx - p;
    }
    else
    {
        p = _index(n_, x_.data(), x_interp);
        dx = x_[p + 1] - x_[p];
        r = (x_interp - x_[p]) / dx;
    }

    const double r2 = r * r;
    const double r3 = r2 * r;
    double wy0 = 0.0, wy1 = 0.0, ws0 = 0.0, ws1 = 0.0;
    int offset = 0;

    if (y_interp)
    {
        wy1 = 3.0 * r2 - 2.0 * r3;
        wy0 = 1.0 - wy1;
        ws0 = (r - 2.0 * r2 + r3) * dx;
        ws1 = (r3 - r2) * dx;
        for (int i = 0; i < n_spline; ++i)
        {
            offset = i_spline[i] * 2 * n_ + p;
            y_interp[i] = wy0 * y_[offset] + wy1 * y_[offset + 1]
                + ws0 * y_[offset + n_] + ws1 * y_[offset + n_ + 1];
        }
    }

    if (dy_interp)
    {
        wy1 = 6.0 * (r - r2) / dx; // wy0 = -wy1
        ws0 = 3.0 * r2 - 4.0 * r + 1.0;
        ws1 = 3.0 * r2 - 2.0 * r;
        for (int i = 0; i < n_spline; ++i)
        {
            offset = i_spline[i] * 2 * n_ + p;
            dy_interp[i] = wy1 * (y_[offset + 1] - y_[offset])
                + ws0 * y_[offset + n_] + ws1 * y_[offset + n_ + 1];
        }
    }

    if (d2y_interp)
    {
        wy1 = (6.0 - 12.0 * r) / (dx * dx); // wy0 = -wy1
        ws0 = (6.0 * r - 4.0) / dx;
        ws1 = (6.0 * r - 2.0) / dx;
        for (int i = 0; i < n_spline; ++i)
        {
            offset = i_spline[i] * 2 * n_ + p;
            d2y_interp[i] = wy1 * (y_[offset + 1] - y_[offset])
                + ws0 * y_[offset + n_] + ws1 * y_[offset + n_ + 1];
        }
    }
}


void CubicSpline::multi_eval(
    double x,
    double* y,
    double* dy,
    double* d2y
) const
{
    std::vector<int> i_spline(n_spline_);
    std::iota(i_spline.begin(), i_spline.end(), 0);
    multi_eval(i_spline.size(), i_spline.data(), x, y, dy, d2y);
}


void CubicSpline::build(
    int n,
    const double* x,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end,
    double* dy
)
{
    std::vector<double> dx(n);
    std::adjacent_difference(x, x + n, dx.begin());
    _build(n, &dx[1], y, bc_start, bc_end, dy);
}


void CubicSpline::build(
    int n,
    double dx,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end,
    double* dy
)
{
    std::vector<double> dx_(n - 1, dx);
    _build(n, dx_.data(), y, bc_start, bc_end, dy);
}


void CubicSpline::eval(
    int n,
    const double* x,
    const double* y,
    const double* dy,
    int n_interp,
    const double* x_interp,
    double* y_interp,
    double* dy_interp,
    double* d2y_interp
)
{
    _validate_eval(n, {}, x, y, dy, n_interp, x_interp);

    // indices of the polynomial segments that contain x_interp
    std::vector<int> _ind(n_interp);
    std::transform(x_interp, x_interp + n_interp, _ind.begin(),
        [n, x](double x_i) { return _index(n, x, x_i); });

    std::vector<double> buffer(n_interp * 5);
    double* _w = buffer.data();
    double* _c0 = _w + n_interp;
    double* _c1 = _c0 + n_interp;
    double* _c2 = _c1 + n_interp;
    double* _c3 = _c2 + n_interp;

    for (int i = 0; i < n_interp; ++i)
    {
        int p = _ind[i];
        double dx = x[p + 1] - x[p];
        double inv_dx = 1.0 / dx;
        double dd = (y[p + 1] - y[p]) * inv_dx;
        _w[i] = x_interp[i] - x[p];
        _c0[i] = y[p];
        _c1[i] = dy[p];
        _c3[i] = (_c1[i] + dy[p + 1] - 2.0 * dd) * inv_dx * inv_dx;
        _c2[i] = (dd - _c1[i]) * inv_dx - _c3[i] * dx;
    }

    _cubic(n_interp, _w, _c0, _c1, _c2, _c3, y_interp, dy_interp, d2y_interp);
}


void CubicSpline::eval(
    int n,
    double x0,
    double dx,
    const double* y,
    const double* dy,
    int n_interp,
    const double* x_interp,
    double* y_interp,
    double* dy_interp,
    double* d2y_interp
)
{
    _validate_eval(n, {x0, dx}, nullptr, y, dy, n_interp, x_interp);

    // indices of the polynomial segments that contain x_interp
    std::vector<int> _ind(n_interp);
    std::transform(x_interp, x_interp + n_interp, _ind.begin(),
        [n, x0, dx](double x_i) { return _index(n, x0, dx, x_i); });

    std::vector<double> buffer(n_interp * 5);
    double* _w = buffer.data();
    double* _c0 = _w + n_interp;
    double* _c1 = _c0 + n_interp;
    double* _c2 = _c1 + n_interp;
    double* _c3 = _c2 + n_interp;

    double inv_dx = 1.0 / dx;
    double inv_dx2 = inv_dx * inv_dx;
    for (int i = 0; i < n_interp; ++i)
    {
        int p = _ind[i];
        double dd = (y[p + 1] - y[p]) * inv_dx;
        _w[i] = x_interp[i] - x0 - p * dx;
        _c0[i] = y[p];
        _c1[i] = dy[p];
        _c3[i] = (_c1[i] + dy[p + 1] - 2.0 * dd) * inv_dx2;
        _c2[i] = (dd - _c1[i]) * inv_dx - _c3[i] * dx;
    }

    _cubic(n_interp, _w, _c0, _c1, _c2, _c3, y_interp, dy_interp, d2y_interp);
}


void CubicSpline::_validate_build(
    int n,
    const double* dx,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end
)
{
    assert(n > 1);

    // if periodic boundary condition is specified, it must be applied to both ends
    assert((bc_start.type == BoundaryType::periodic)
            == (bc_end.type == BoundaryType::periodic));

    // y[0] must equal y[n-1] for periodic boundary condition
    assert(bc_start.type != BoundaryType::periodic || y[0] == y[n - 1]);

    // not-a-knot boundary condition requires the existence of "internal" knot
    // so n must be at least 3
    assert((bc_start.type != BoundaryType::not_a_knot && 
            bc_end.type != BoundaryType::not_a_knot) || n > 2);

    // knots must be strictly increasing
    assert(std::all_of(dx, dx + n - 1, [](double d) { return d > 0.0; }));
}


void CubicSpline::_validate_eval(
    int n,
    const double (&u)[2],
    const double* x,
    const double* y,
    const double* dy,
    int n_interp,
    const double* x_interp
)
{
    assert(n > 1 && y && dy);
    assert((x && std::is_sorted(x, x + n, std::less_equal<double>())) || u[1] > 0.0);

    assert((n_interp > 0 && x_interp) || n_interp == 0);

    double xmin = x ? x[0] : u[0];
    double xmax = x ? x[n - 1] : u[0] + (n - 1) * u[1];
    assert(std::all_of(x_interp, x_interp + n_interp,
                       [xmin, xmax](double x_i) { return xmin <= x_i && x_i <= xmax; }));
}


void CubicSpline::_build(
    int n,
    const double* dx,
    const double* y,
    const BoundaryCondition& bc_start,
    const BoundaryCondition& bc_end,
    double* dy
)
{
    _validate_build(n, dx, y, bc_start, bc_end);

    if (n == 2 && bc_start.type == BoundaryType::periodic)
    {
        dy[0] = dy[1] = 0.0; // the only possible solution: constant
    }
    else if (n == 3 && bc_start.type == BoundaryType::not_a_knot
                    && bc_end.type == BoundaryType::not_a_knot)
    {
        // in this case two conditions coincide
        // simply build a parabola that passes through the three data points
        double dd01 = (y[1] - y[0]) / dx[0]; // divided difference f[x0,x1]
        double dd12 = (y[2] - y[1]) / dx[1]; // f[x1,x2]
        double dd012 = (dd12 - dd01) / (dx[0] + dx[1]); // f[x0,x1,x2]

        dy[0] = dd01 - dd012 * dx[0];
        dy[1] = 2.0 * dd01 - dy[0];
        dy[2] = dd01 + dd012 * (dx[0] + 2.0 * dx[1]);
    }
    else
    {
        std::vector<double> buffer(4 * n);

        double* dd = buffer.data(); // divided differences
        std::adjacent_difference(y, y + n, dd);
        dd += 1; // the first element computed by adjacent_difference is not a difference
        std::transform(dd, dd + n - 1, dx, dd, std::divides<double>());

        // tridiagonal linear system (cyclic tridiagonal if periodic boundary condition)
        double* d = buffer.data() + n; // main diagonal
        double* l = buffer.data() + 2 * n; // subdiagonal
        double* u = buffer.data() + 3 * n; // superdiagonal

        //***********************************************
        // common part of the tridiagonal linear system
        //***********************************************
        std::copy(dx + 1, dx + n - 1, l);
        std::copy(dx, dx + n - 2, u + 1);

        for (int i = 1; i != n - 1; ++i)
        {
            d[i] = 2.0 * (dx[i - 1] + dx[i]);
            dy[i] = 3.0 * (dd[i - 1] * dx[i] + dd[i] * dx[i - 1]);
        }

        //***********************************************
        //          boundary-specific part
        //***********************************************
        if (bc_start.type == BoundaryType::periodic)
        {
            // exclude s[n-1] and solve a a cyclic tridiagonal linear system of size n-1
            d[0] = 2.0 * (dx[n - 2] + dx[0]);
            u[0] = dx[n - 2];
            l[n - 2] = dx[0];
            dy[0] = 3.0 * (dd[0] * dx[n - 2] + dd[n - 2] * dx[0]);
            _solve_cyctri(n - 1, d, u, l, dy);
            dy[n - 1] = dy[0];
        }
        else
        {
            switch (bc_start.type)
            {
            case BoundaryType::first_deriv:
                d[0] = 1.0 * dx[0];
                u[0] = 0.0;
                dy[0] = bc_start.val * dx[0];
                break;
            case BoundaryType::second_deriv:
                d[0] = 2.0 * dx[0];
                u[0] = 1.0 * dx[0];
                dy[0] = (3.0 * dd[0] - 0.5 * bc_start.val * dx[0]) * dx[0];
                break;
            default: // BoundaryCondition::not_a_knot
                d[0] = dx[1];
                u[0] = dx[0] + dx[1];
                dy[0] = (dd[0] * dx[1] * (dx[0] + 2 * u[0]) + dd[1] * dx[0] * dx[0]) / u[0];
            }

            switch (bc_end.type)
            {
            case BoundaryType::first_deriv:
                d[n - 1] = 1.0 * dx[n - 2];
                l[n - 2] = 0.0;
                dy[n - 1] = bc_end.val * dx[n - 2];
                break;
            case BoundaryType::second_deriv:
                d[n - 1] = 2.0 * dx[n - 2];
                l[n - 2] = 1.0 * dx[n - 2];
                dy[n - 1] = (3.0 * dd[n - 2] + 0.5 * bc_end.val * dx[n - 2]) * dx[n - 2];
                break;
            default: // BoundaryCondition::not_a_knot
                d[n - 1] = dx[n - 3];
                l[n - 2] = dx[n - 3] + dx[n - 2];
                dy[n - 1] = (dd[n - 2] * dx[n - 3] * (dx[n - 2] + 2 * l[n - 2])
                             + dd[n - 3] * dx[n - 2] * dx[n - 2]) / l[n - 2];
            }

            int nrhs = 1;
            int ldb = n;
            int info = 0;
            dgtsv_(&n, &nrhs, l, d, u, dy, &ldb, &info);
        }
    }
}


void CubicSpline::_cubic(
    int n,
    const double* w,
    const double* c0,
    const double* c1,
    const double* c2,
    const double* c3,
    double* y,
    double* dy,
    double* d2y
)
{
    if (y)
    {
        for (int i = 0; i < n; ++i)
        {
            y[i] = ((c3[i] * w[i] + c2[i]) * w[i] + c1[i]) * w[i] + c0[i];
        }
    }

    if (dy)
    {
        for (int i = 0; i < n; ++i)
        {
            dy[i] = (3.0 * c3[i] * w[i] + 2.0 * c2[i]) * w[i] + c1[i];
        }
    }

    if (d2y)
    {
        for (int i = 0; i < n; ++i)
        {
            d2y[i] = 6.0 * c3[i] * w[i] + 2.0 * c2[i];
        }
    }
}


int CubicSpline::_index(int n, const double* knots, double x)
{
    int i = (std::upper_bound(knots, knots + n, x) - knots) - 1;
    return i - (i == n - 1);
}


int CubicSpline::_index(int n, double x0, double dx, double x)
{
    int i = (x - x0) / dx;
    return i - (i == n - 1);
}


void CubicSpline::_solve_cyctri(int n, double* d, double* u, double* l, double* b)
{
    // flexible non-zero parameters that can affect the condition number of the
    // tridiagonal linear system
    double alpha = 1.0;
    double beta = -d[0] / u[n - 1];

    std::vector<double> bp(2 * n, 0.0);
    std::copy(b, b + n, bp.begin());
    bp[n] = 1. / alpha;
    bp[2 * n - 1] = 1. / beta;

    d[0] -= u[n - 1] * beta / alpha;
    d[n - 1] -= l[n - 1] * alpha / beta;

    int nrhs = 2;
    int info = 0;
    int ldb = n;
    dgtsv_(&n, &nrhs, l, d, u, bp.data(), &ldb, &info);

    double fac = (beta * u[n - 1] * bp[0] + alpha * l[n - 1] * bp[n - 1])
                 / (1. + beta * u[n - 1] * bp[n] + alpha * l[n - 1] * bp[2 * n - 1]);

    std::transform(bp.begin(), bp.begin() + n, bp.begin() + n, b,
            [fac](double yi, double zi) { return yi - fac * zi; });
}



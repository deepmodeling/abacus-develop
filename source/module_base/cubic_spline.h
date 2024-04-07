#ifndef CUBIC_SPLINE_INTERPOLATION_H_
#define CUBIC_SPLINE_INTERPOLATION_H_

#include "lapack_connector.h"
#include <algorithm>
#include <chrono>
#include <vector>
#include <chrono>

using iclock = std::chrono::high_resolution_clock;

namespace ModuleBase
{

/**
 * @brief Cubic spline interplation.
 *
 * This class interpolates a given set of data points (x[i], y[i]) (i=0,...,n-1)
 * by piecewise cubic polynomials with continuous first and second derivatives
 * at x[i] ("knots").
 *
 * There are two ways to use this class. The first way uses the class as an
 * interpolant object; the second way merely uses static member functions.
 *
 * Usage-1: interpolant object
 *
 *      // build the interpolant object
 *      // n is the number of data points (x[i], y[i]) (i=0,...,n-1)
 *      CubicSpline cubspl(n, x, y);
 *
 *      // evaluates the interpolant at x_interp
 *      cubspl.eval(n_interp, x_interp, y_interp); // values are returned in y_interp
 *
 *      // evaluates both the interpolant and its derivative at x_interp
 *      cubspl.eval(n_interp, x_interp, y_interp, s_interp);
 *
 *      // evaluates the derivative only
 *      cubspl.eval(n_interp, x_interp, nullptr, s_interp);
 *
 *      // for evenly spaced x (x[i] = x0 + i*dx), x0 and dx can be passed instead of x
 *      CubicSpline cubspl2(n, nullptr, y, dx, x0);
 *
 *      // "not-a-knot" boundary condition is used by default at both ends
 *      // other supported boundary conditions include first/second derivatives and
 *      // periodic boundary condition.
 *
 *      // not-a-knot, first/second_deriv can be independently applied to each end
 *      // for example, to specify f''(start) = 1.0 and f'(end) = 3.0,
 *      CubicSpline cubspl3(..., CubicSpline::BoundaryCondition::second_deriv,
 *                               CubicSpline::BoundaryCondition::first_deriv,
 *                               1.0, 3.0);
 *
 *      // periodic boundary condition needs to be specified at both ends
 *      // and y[0] must equal y[n-1]
 *      CubicSpline cubspl4(..., CubicSpline::BoundaryCondition::periodic,
 *                               CubicSpline::BoundaryCondition::periodic);
 *
 * Usage-2: static member functions (no interpolant object)
 *
 *      // step-1: computes the first-derivatives at knots
 *      CubicSpline::compute(n, x, y, s); // s is updated with the first derivatives
 *
 *      // various boundary conditions and explicitly specified evenly spaced
 *      // knots are supported in the same way as the interpolant object
 *      // now n, x, y, s are sufficient to define a piecewise cubic polynomial
 *
 *      // step-2: interpolates with knots, values & derivatives
 *      CubicSpline::interp(n, x, y, s, n_interp, x_interp, y_interp, s_interp);
 *
 *
 */
class CubicSpline
{
public:
    CubicSpline()                   = delete;
    CubicSpline(CubicSpline const&) = default;
    CubicSpline(CubicSpline &&)     = default;

    CubicSpline& operator=(CubicSpline const&)  = default;
    CubicSpline& operator=(CubicSpline &&)      = default;

    ~CubicSpline() = default; 

    /**
     * @brief Boundary conditions of cubic spline interpolations.
     *
     * Available boundary conditions include:
     * - not_a_knot:   the first two pieces at the start or the last two at the end
     *                 are the same polynomial, i.e., x[1] or x[n-2] is not a "knot";
     * - first_deriv:  User-defined first derivative;
     * - second_deriv: User-defined second derivative;
     * - periodic:     The first and second derivatives at two ends are continuous.
     *                 This condition requires that y[0] = y[n-1] and it must be
     *                 applied to both ends.
     */
    enum class BoundaryCondition
    {
        not_a_knot,
        first_deriv,
        second_deriv,
        periodic
    };

    /**
     * @brief Builds the interpolant object.
     *
     * This constructor computes the first derivatives at x ("knots"), and stores, together
     * with x & y, in the object.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points, must be strictly increasing
     * @param[in]   y               y coordiates of data points
     * @param[in]   all_poly_coef   if true, all polynomial coefficients are stored in the object,
     *                              which will accelerates the evaluation of the interpolant by
     *                              about 1/3 but increases the memory usage by about 2/3.
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     * @param[in]   deriv_start     first or second derivative at the start,
     *                              ignored if bc_start is not first_deriv or second_deriv
     * @param[in]   deriv_end       first or second derivative at the end,
     *                              ignored if bc_end is not first_deriv or second_deriv
     */
    CubicSpline(
        const int n,
        const double* const x,
        const double* const y,
        //const bool all_poly_coef = false,
        BoundaryCondition bc_start = BoundaryCondition::not_a_knot,
        BoundaryCondition bc_end = BoundaryCondition::not_a_knot,
        const double deriv_start = 0.0,
        const double deriv_end = 0.0
    );

    /**
     * @brief Builds the interpolant object with evenly-spaced knots.
     *
     * This constructor computes the first derivatives at x[i] = x0 + i*dx, and stores,
     * together with x & y, in the object.
     *
     * @param[in]   n               number of data points
     * @param[in]   x0              starting x coordinate
     * @param[in]   dx              spacing between knots
     * @param[in]   y               y coordiates of data points
     * @param[in]   all_poly_coef   if true, all polynomial coefficients are stored in the object,
     *                              which will accelerates the evaluation of the interpolant by
     *                              about 1/3 but doubles the memory usage.
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     * @param[in]   deriv_start     first or second derivative at the start,
     *                              ignored if bc_start is not first_deriv or second_deriv
     * @param[in]   deriv_end       first or second derivative at the end,
     *                              ignored if bc_end is not first_deriv or second_deriv
     */
    CubicSpline(
        const int n,
        const double x0,
        const double dx,
        const double* const y,
        //const bool all_poly_coef = false,
        BoundaryCondition bc_start = BoundaryCondition::not_a_knot,
        BoundaryCondition bc_end = BoundaryCondition::not_a_knot,
        const double deriv_start = 0.0,
        const double deriv_end = 0.0
    );
    /**
     * @brief Evaluates the interpolant object.
     *
     * This function evaluates the interpolant and its first derivative at x_interp.
     *
     * @param[in]   n_interp        number of points to evaluate the interpolant
     * @param[in]   x_interp        places where the interpolant is evaluated, must be
     *                              within [x_[0], x_[n-1]] or [x0_, x0_+(n-1)*dx_]
     * @param[out]  y_interp        interpolated values
     * @param[out]  Dy_interp       derivatives at x_interp
     *
     * @note If y_interp or s_interp is nullptr, the corresponding calculation is skipped.
     */
    void eval(
        const int n_interp,
        const double* const x_interp,
        double* const y_interp,
        double* const Dy_interp = nullptr//,
        //double* const D2y_interp = nullptr
    ) const;

    /// Returns the object's heap memory usage in bytes.
    //size_t heap_memory_usage() const { return (x_.capacity() + c_.capacity()) * sizeof(double); }

    /**
     * @brief Computes the first derivatives at knots for cubic spline interpolation.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points ("knots"), must be
     *                              strictly increasing. if x is nullptr, knots are
     *                              assumed to be uniformly spaced with spacing dx.
     * @param[in]   y               y coordiates of data points
     * @param[out]  dy              first derivatives at knots
     * @param[out]  d2y             second derivatives at knots
     * @param[out]  d3y             third derivatives at knots
     * @param[in]   dx_uniform      spacing between knots (ignored if x is not nullptr)
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     * @param[in]   deriv_start     first or second derivative at the start,
     *                              ignored if bc_start is not first_deriv or second_deriv
     * @param[in]   deriv_end       first or second derivative at the end,
     *                              ignored if bc_end is not first_deriv or second_deriv
     */
    static void build(
        const int n,
        const double* const x,
        const double* const y,
        double* const dy,
        //double* const d2y = nullptr,
        //double* const d3y = nullptr,
        const double dx_uniform = -1.0,
        BoundaryCondition bc_start = BoundaryCondition::not_a_knot,
        BoundaryCondition bc_end = BoundaryCondition::not_a_knot,
        const double deriv_start = 0.0,
        const double deriv_end = 0.0
    );

    /**
     * @brief Performs a cubic spline interpolation.
     *
     * Given knots, values and first derivatives at knots, this function evaluates
     * values and first derivatives at x_interp with cubic spline interpolation.
     *
     * @param[in]   n               number of knots
     * @param[in]   x               knots. if x is nullptr, knots are assumed to be
     *                              uniformly spaced as x0 + i*dx (i=0,...,n-1).
     * @param[in]   y               values at knots
     * @param[in]   s               first derivatives at knots
     * @param[in]   n_interp        number of points to evaluate the spline polynomial
     * @param[in]   x_interp        places where the spline polynomial is evaluated,
     *                              must be within [x[0], x[n-1]]
     * @param[out]  y_interp        interpolated values at x_interp
     * @param[out]  s_interp        derivatives at x_interp
     * @param[in]   dx              spacing between knots (ignored if x is not nullptr)
     * @param[in]   x0              starting x coordinate (ignored if x is not nullptr)
     *
     * @note If y_interp or s_interp is nullptr, the corresponding calculation is skipped.
     *
     */
    static void eval(
        const int n,
        const double* const x,
        const double* const y,
        const double* const dy,
        //const double* const d2y,
        //const double* const d3y,
        const int n_interp,
        const double* const x_interp,
        double* const y_interp,
        double* const dy_interp = nullptr,
        //double* const d2y_interp = nullptr,
        const double x0 = 0.0,
        const double dx = -1.0
    );

    static std::chrono::duration<double> time_alloc;
    static std::chrono::duration<double> time_seg;
    static std::chrono::duration<double> time_dd;
    static std::chrono::duration<double> time_c23;
    static std::chrono::duration<double> time_c01;
    static std::chrono::duration<double> time_poly;


private:
    /// number of knots
    int n_ = -1;

    /// starting x coordinate (ignored if all knots are explicitly specified in build)
    double x0_ = 0.0;

    /// spacing between knots (ignored if all knots are explicitly specified in build)
    double dx_ = -1.0;

    /// knots (will remain empty if x0 & dx are specified in build)
    std::vector<double> x_;

    /**
     * @brief coefficients of cubic polynomials.
     *
     */
    //std::vector<double> c_;
    std::vector<double> y_;
    std::vector<double> s_;

    /// Asserts that the input arguments are valid for computing cubic spline.
    /// This function is skipped when __DEBUG is not defined.
    static void _validate_build(
        const int n,
        const double* const x,
        const double* const y,
        const double dx,
        BoundaryCondition bc_start,
        BoundaryCondition bc_end
    );

    /// Asserts that the input arguments are valid for an interpolation.
    /// This function is skipped when __DEBUG is not defined.
    static void _validate_eval(
        const int n,
        const double* const x,
        const double* const y,
        const double* const s,
        const int n_interp,
        const double* const x_interp,
        const double x0,
        const double dx
    );

    /**
     * @brief Solves a cyclic tridiagonal linear system.
     *
     * This function solves a cyclic tridiagonal linear system A*x=b where b
     * is a vector and A is given by
     *
     *      D[0]   U[0]                           L[n-1]
     *      L[0]   D[1]   U[1]
     *             L[1]   D[2]   U[2]
     *                    ...    ...      ...
     *                          L[n-3]   D[n-2]   U[n-2]
     *      U[n-1]                       L[n-2]   D[n-1]
     *
     * On finish, b is overwritten by the solution.
     *
     * Sherman-Morrison formula is used to convert the problem into a tridiagonal
     * linear system, after which the problem can be solved by dgtsv efficiently.
     *
     * @param[in]       n       size of the linear system
     * @param[in]       D       main diagonal
     * @param[in]       U       superdiagonal
     * @param[in]       L       subdiagonal
     * @param[in,out]   b       right hand side of the linear system, will be
     *                          overwritten by the solution on finish.
     *
     * @note D, L, U are all overwritten in this function. Use with care!
     *
     */
    static void _solve_cyctri(
        const int n,
        double* const D,
        double* const U,
        double* const L,
        double* const b
    );

    /**
     * @brief Index of the polynomial segment containing x_interp.
     *
     * Assuming a strictly increasing x and x[0] <= x_interp <= x[n-1],
     * this function returns n-2 if x_interp == x[n-1], otherwise an index i
     * such that x[i] <= x_interp < x[i+1].
     *
     * @note No check is performed for the above assumptions. Use with care!
     *
     * @param[in]   n           number of knots
     * @param[in]   x           knots of the interpolant
     * @param[in]   x_interp    place where the interpolant is evaluated
     *
     */
    static int _index(const int n, const double* const x, double x_interp)
    {
        int i = (std::upper_bound(x, x + n, x_interp) - x) - 1;
        return i - (i == n - 1);
    }

    /// Index of the polynomial segment containing x_interp (for evenly spaced knots).
    static int _index(const int n, const double x0, const double dx, double x_interp)
    {
        int i = (x_interp - x0) / dx;
        return i - (i == n - 1);
    }
};

}; // namespace ModuleBase

#endif

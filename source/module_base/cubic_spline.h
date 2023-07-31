#ifndef CUBIC_SPLINE_INTERPOLATOR_H_
#define CUBIC_SPLINE_INTERPOLATOR_H_

#include "lapack_connector.h"

namespace ModuleBase
{

//! A class that performs cubic spline interplations.
/*!
 *  This class interpolates a given set of data points (x[i],y[i]) (i=0,...,n-1)
 *  by piecewise cubic polynomials with continuous first and second derivatives
 *  at x[i].
 *
 *  Usage:
 *
 *      CubicSpline cubspl;
 *
 *      // build the interpolant
 *      // n is the number of data points (x[i],y[i]) (i=0,...,n-1)
 *      // by default "not-a-knot" boundary condition is used at both ends
 *      cubspl.build(n, x, y);
 *
 *      // alternatively one can specify the first or second derivatives at the ends
 *      cubspl.build(n, x, y, CubicSpline::BoundaryCondition::second_deriv,
 *                      CubicSpline::BoundaryCondition::second_deriv, 1.0, -2.0);
 *
 *      // not-a-knot, first_deriv & second_deriv can be independently applied to
 *      // each end; two ends do not necessarily have the same type of condition
 *      cubspl.build(n, x, y, CubicSpline::BoundaryCondition::first_deriv,
 *                      CubicSpline::BoundaryCondition::not-a-knot, 1.0);
 *
 *      // periodic boundary condition needs to be specified at both ends
 *      // and y[0] must equal y[n-1] in this case
 *      cubspl.build(n, x, y, CubicSpline::BoundaryCondition::periodic,
 *                            CubicSpline::BoundaryCondition::periodic );
 *
 *      // calculate the values of interpolant at x_interp[i]
 *      cubspl.interp(n_interp, x_interp, y_interp);
 *
 *      // calculate the values & derivatives of interpolant at x_interp[i]
 *      cubspl.interp(n_interp, x_interp, y_interp, dy_interp);
 *                                                                                  */
class CubicSpline
{
  public:
    CubicSpline(){};
    CubicSpline(CubicSpline const&) = delete;
    CubicSpline& operator=(CubicSpline const&) = delete;

    ~CubicSpline();

    //! Boundary condition of cubic spline interpolations.
    /*!
     *  Available boundary conditions include:
     *  - not_a_knot:   the first two pieces at the start or the last two at the end
     *                  are the same polynomial, i.e., x[1] or x[n-2] is not a "knot";
     *  - first_deriv:  User-defined first derivative;
     *  - second_deriv: User-defined second derivative;
     *  - periodic:     The first and second derivatives at two ends are continuous.
     *                  This condition requires that y[0] = y[n-1] and it must be
     *                  applied to both ends.
     *                                                                              */
    enum class BoundaryCondition
    {
        not_a_knot,
        first_deriv,
        second_deriv,
        periodic
    };

    //! Builds the interpolant.
    /*!
     *  This function computes all cubic spline polynomial coefficients from given
     *  data points and boundary conditions.
     *                                                                              */
    void build(const int n,           //!< [in] number of data points
               const double* const x, //!< [in] x coordiates of data points, must be strictly increasing
               const double* const y, //!< [in] y coordiates of data points
               BoundaryCondition bc_start = BoundaryCondition::not_a_knot, //!< [in] boundary condition at the start
               BoundaryCondition bc_end = BoundaryCondition::not_a_knot,   //!< [in] boundary condition at the end
               const double deriv_start = 0.0, //!< [in] first or second derivative at the start,
                                               //!<      used if bc_start is first_deriv or second_deriv
               const double deriv_end = 0.0    //!< [in] first or second derivative at the end,
                                               //!<      used if bc_end is first_deriv or second_deriv
    );

    //! Calculates the values of interpolant.
    /*!
     *  This function evaluates the interpolant at x_interp[i].
     *  On finish, interpolated values are placed in y_interp,
     *  and the derivatives at x_interp[i] are placed in dy_interp.
     *
     *  If y_interp or dy_interp is nullptr, the corresponding values are
     *  not calculated. y_interp and dy_interp must not be both nullptr.
     *                                                                              */
    void interp(const int n,                       //!< [in]  number of points to evaluate the interpolant
                const double* const x_interp,      //!< [in]  places where the interpolant is evaluated;
                                                   //!<       must be within [x_[0], x_[n-1]]
                double* const y_interp,            //!< [out] interpolated values
                double* const dy_interp = nullptr  //!< [out] derivatives at x_interp
    );

    /// knots of the interpolant
    const double* x() const { return x_; }

    /// values at knots
    const double* y() const { return y_; }

    /// first-order derivatives at knots
    const double* s() const { return s_; }

  private:
    //! number of data points
    int n_ = 0;

    //! knots (x coordinates of data points)
    double* x_ = nullptr;

    //! values at knots
    double* y_ = nullptr;

    //! first-order derivatives at knots
    double* s_ = nullptr;

    //! A flag that tells whether the knots are evenly spaced.
    bool is_uniform_ = false;

    //! Numerical threshold for determining whether the knots are evenly spaced.
    /*!
     *  The knots are considered uniform (evenly spaced) if for every i from 0 to n-2
     *
     *          abs( (x[i+1]-x[i]) - (x[n-1]-x[0])/(n-1) ) < uniform_thr_
     *                                                                              */
    double uniform_thr_ = 1e-15;

    //! Checks whether the input arguments of build() are valid
    void sanity_check(const int n,
                      const double* const x,
                      const double* const y,
                      BoundaryCondition bc_start,
                      BoundaryCondition bc_end);

    //! Solves a cyclic tridiagonal linear system.
    /*!
     *  This function solves a cyclic tridiagonal linear system A*x=b where b
     *  is a vector and A is given by
     *
     *      D[0]   U[0]                           L[n-1]
     *      L[0]   D[1]   U[1]
     *             L[1]   D[2]   U[2]
     *                    ...    ...      ...
     *                          L[n-3]   D[n-2]   U[n-2]
     *      U[n-1]                       L[n-2]   D[n-1]
     *
     *  On finish, b is overwritten by the solution.
     *
     *  Sherman-Morrison formula is used to convert the problem into a tridiagonal
     *  linear system, after which the problem can be solved by dgtsv efficiently.
     *
     *  @note D, L, U are all modified in this function, so use with care!
     *                                                                              */
    void solve_cyctri(const int n,     //!< [in] size of the linear system
                      double* const D, //!< [in] main diagonal
                      double* const U, //!< [in] superdiagonal
                      double* const L, //!< [in] subdiagonal
                      double* const b  //!< [in,out] right hand side of the linear system;
                                       //!< will be overwritten by the solution on finish.
    );

    //! Wipes off the interpolant (if any) and deallocates memories.
    void cleanup();

    static void _eval_y(double w, double c0, double c1, double c2, double c3, double* y, double*)
    {
        *y = ((c3 * w + c2) * w + c1) * w + c0;
    }

    static void _eval_dy(double w, double, double c1, double c2, double c3, double*, double* dy)
    {
        *dy = (3.0 * c3 * w + 2.0 * c2) * w + c1;
    }

    static void _eval_y_dy(double w, double c0, double c1, double c2, double c3, double* y, double* dy)
    {
        *y = ((c3 * w + c2) * w + c1) * w + c0;
        *dy = (3.0 * c3 * w + 2.0 * c2) * w + c1;
    }
};

}; // namespace ModuleBase

#endif

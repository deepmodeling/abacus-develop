#ifndef CUBIC_SPLINE_INTERPOLATOR_H_
#define CUBIC_SPLINE_INTERPOLATOR_H_

#include "lapack_connector.h"

namespace ModuleBase
{

/*!
 * @brief A class that performs cubic spline interplations.
 *
 * This class interpolates a given set of data points (x[i],y[i]) (i=0,...,n-1)
 * by piecewise cubic polynomials with continuous first and second derivatives
 * at x[i].
 *
 * There are two ways to use this class. The first way treats the class as an
 * interpolant object, and the second way uses the static member functions.
 *
 * Usage-1: interpolant object
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
 *      cubspl.eval(n_interp, x_interp, y_interp);
 *
 *      // calculate the values & derivatives of interpolant at x_interp[i]
 *      cubspl.eval(n_interp, x_interp, y_interp, dy_interp);
 *
 * Usage-2: static member functions
 *
 *      // gets the first-derivatives (s) at knots
 *      // may build with various boundary conditions as above
 *      CubicSpline::build(n, x, y, s);
 *
 *      // evaluates the interpolant with knots, values & derivatives
 *      CubicSpline::eval(n, x, y, s, n_interp, x_interp, y_interp, dy_interp);
 *
 *                                                                                 */
class CubicSpline
{
  public:
    CubicSpline(){};
    CubicSpline(CubicSpline const&);
    CubicSpline& operator=(CubicSpline const&);

    ~CubicSpline() { cleanup(); }

    /*!
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
     *                                                                             */
    enum class BoundaryCondition
    {
        not_a_knot,
        first_deriv,
        second_deriv,
        periodic
    };

    /*!
     * @brief Builds the interpolant.
     *
     * By calling this function, the class object computes and stores the first-order
     * derivatives at knots by from given data points and boundary conditions, which
     * makes the object an interpolant.
     *                                                                              */
    void build(const int n,           //!< [in] number of data points
               const double* const x, //!< [in] x coordiates of data points, must be strictly increasing
               const double* const y, //!< [in] y coordiates of data points
               BoundaryCondition bc_start = BoundaryCondition::not_a_knot, //!< [in] boundary condition at start
               BoundaryCondition bc_end = BoundaryCondition::not_a_knot,   //!< [in] boundary condition at end
               const double deriv_start = 0.0, //!< [in] first or second derivative at the start,
                                               //!<      used if bc_start is first_deriv or second_deriv
               const double deriv_end = 0.0    //!< [in] first or second derivative at the end,
                                               //!<      used if bc_end is first_deriv or second_deriv
    );

    /*!
     * @brief Evaluates the interpolant.
     *
     * This function evaluates the interpolant at x_interp[i].
     * On finish, interpolated values are placed in y_interp,
     * and the derivatives at x_interp[i] are placed in dy_interp.
     *
     * If y_interp or dy_interp is nullptr, the corresponding values are
     * not calculated. They must not be nullptr at the same time.
     *
     * @note the interpolant must be built before calling this function.
     *                                                                              */
    void eval(const int n,                      //!< [in]  number of points to evaluate the interpolant
              const double* const x_interp,     //!< [in]  places where the interpolant is evaluated;
                                                //!<       must be within [x_[0], x_[n-1]]
              double* const y_interp,           //!< [out] interpolated values
              double* const dy_interp = nullptr //!< [out] derivatives at x_interp
    );

    /// knots of the interpolant
    const double* x() const { return x_; }

    /// values at knots
    const double* y() const { return y_; }

    /// first-order derivatives at knots
    const double* s() const { return s_; }

    static void build(const int n,           //!< [in] number of data points
                      const double* const x, //!< [in] x coordiates of data points, must be strictly increasing
                      const double* const y, //!< [in] y coordiates of data points
                      double* const s,       //!< [out] first-order derivatives at knots
                      BoundaryCondition bc_start = BoundaryCondition::not_a_knot, //!< [in] boundary condition at start
                      BoundaryCondition bc_end = BoundaryCondition::not_a_knot,   //!< [in] boundary condition at end
                      const double deriv_start = 0.0, //!< [in] first or second derivative at the start,
                                                      //!<      used if bc_start is first_deriv or second_deriv
                      const double deriv_end = 0.0    //!< [in] first or second derivative at the end,
                                                      //!<      used if bc_end is first_deriv or second_deriv
    );

    static void eval(const int n,                      //!< [in]  number of knots
                     const double* const x,            //!< [in]  knots of the interpolant
                     const double* const y,            //!< [in]  values at knots
                     const double* const s,            //!< [in]  first-order derivatives at knots
                     const int n_interp,               //!< [in]  number of points to evaluate the interpolant
                     const double* const x_interp,     //!< [in]  places where the interpolant is evaluated;
                                                       //!<       must be within [x_[0], x_[n-1]]
                     double* const y_interp,           //!< [out] interpolated values
                     double* const dy_interp = nullptr //!< [out] derivatives at x_interp
    );

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

    /// Checks whether the input arguments are valid for building a cubic spline.
    static void check_build(const int n,
                            const double* const x,
                            const double* const y,
                            BoundaryCondition bc_start,
                            BoundaryCondition bc_end);

    /// Checks whether the input arguments are valid for evaluating a cubic spline.
    static void check_interp(const int n,
                             const double* const x,
                             const double* const y,
                             const double* const s,
                             const int n_interp,
                             const double* const x_interp,
                             double* const y_interp,
                             double* const dy_interp);

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
    static void solve_cyctri(const int n,     //!< [in] size of the linear system
                             double* const D, //!< [in] main diagonal
                             double* const U, //!< [in] superdiagonal
                             double* const L, //!< [in] subdiagonal
                             double* const b  //!< [in,out] right hand side of the linear system;
                                              //!< will be overwritten by the solution on finish.
    );

    //! Wipes off the interpolant (if any) and deallocates memories.
    void cleanup();

    /// Evaluates a cubic polynomial and its derivative.
    template <bool EvalY, bool EvalDy>
    static void _poly_eval(double w, double c0, double c1, double c2, double c3, double* y, double* dy)
    {
        if (EvalY)
        {
            *y = ((c3 * w + c2) * w + c1) * w + c0;
        }

        if (EvalDy)
        {
            *dy = (3.0 * c3 * w + 2.0 * c2) * w + c1;
        }
    }

    /*!
     * @brief Generates a function that returns the index of the left knot of the
     *        spline polynomial to be evaluated.
     *
     * This function takes the knots of the interpolant and returns a function that
     * takes a value x_interp and returns the index of the left knot of the spline
     * polynomial. If "is_uniform" is not 0/1, this function will checks whether
     * the knots are evenly spaced. The returned function makes use of "is_uniform"
     * to speed up the search.
     *                                                                              */
    static std::function<int(double)> _gen_search(const int n, const double* const x, int is_uniform = -1);

    /*!
     * @brief Evaluates a cubic spline with given knots, values and derivatives.
     *                                                                              */
    static void _eval(const double* const x,            //!< [in]  knots of the interpolant
                      const double* const y,            //!< [in]  values at knots
                      const double* const s,            //!< [in]  first-order derivatives at knots
                      const int n_interp,               //!< [in]  number of points to evaluate the interpolant
                      const double* const x_interp,     //!< [in]  places where the interpolant is evaluated;
                                                        //!<       must be within [x_[0], x_[n-1]]
                      double* const y_interp,           //!< [out] interpolated values
                      double* const dy_interp,          //!< [out] derivatives at x_interp
                      std::function<int(double)> search //!< [in]  a function that returns the index of the left
                                                        //         knot of the spline polynomial to be evaluated
    );
};

}; // namespace ModuleBase

#endif

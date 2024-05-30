#ifndef CUBIC_SPLINE_INTERPOLATION_H_
#define CUBIC_SPLINE_INTERPOLATION_H_

#include <vector>
#include <cstddef>

namespace ModuleBase
{

/**
 * @brief Cubic spline interplation.
 *
 * Interpolating a set of data points (x[i], y[i]) (i=0,...,n-1) by piecewise
 * cubic polynomials
 *
 *      p_i(x) = c0[i] + c1[i]*(x-x[i]) + c2[i]*(x-x[i])^2 + c3[i]*(x-x[i])^3
 *
 * with continuous first and second derivatives. (p_i(x) is defined on the
 * interval [x[i], x[i+1]])
 *
 * There are two ways to use this class. The first way treats class objects as
 * interpolants; the second way uses static member functions.
 *
 * Usage-1: object as interpolant
 *
 *      //---------------------------------------------------------------------
 *      //                          basic usage
 *      //---------------------------------------------------------------------
 *      // build the interpolant object
 *      // n is the number of data points (x[i], y[i]) (i=0,...,n-1)
 *      CubicSpline cubspl(n, x, y);
 *
 *      // evaluates the interpolant at multiple places (x_interp)
 *      // n_interp is the number of places to evaluate the interpolant
 *      cubspl.eval(n_interp, x_interp, y_interp); // values are returned in y_interp
 *
 *      // evaluates both the values and first derivatives at x_interp
 *      cubspl.eval(n_interp, x_interp, y_interp, dy_interp);
 *
 *      // evaluates the second derivative only
 *      cubspl.eval(n_interp, x_interp, nullptr, nullptr, d2y_interp);
 *
 *      //---------------------------------------------------------------------
 *      //                      evenly spaced knots
 *      //---------------------------------------------------------------------
 *      // Interpolants with evenly spaced knots can be built by a different
 *      // constructor, which allows faster evaluation due to quicker index lookup.
 *
 *      // build an interpolant with evenly spaced knots x[i] = x0 + i*dx
 *      CubicSpline cubspl(n, x0, dx, y);
 *
 *      //---------------------------------------------------------------------
 *      //                      boundary conditions
 *      //---------------------------------------------------------------------
 *      // By default "not-a-knot" boundary condition is applied to both ends.
 *      // Other supported boundary conditions include first/second derivatives
 *      // and periodic boundary condition.
 *
 *      // build an interpolant with f''(start) = 1.0 and f'(end) = 3.0
 *      CubicSpline cubspl(n, x, y,
 *                         {CubicSpline::BoundaryType::second_deriv, 1.0},
 *                         {CubicSpline::BoundaryType::first_deriv, 3.0});
 *
 *      // build an interpolant with periodic boundary condition
 *      CubicSpline cubspl(n, x, y, // y[0] must equal y[n-1]
 *                         {CubicSpline::BoundaryType::periodic},
 *                         {CubicSpline::BoundaryType::periodic});
 *
 *      //---------------------------------------------------------------------
 *      //                      multiple interpolants
 *      //---------------------------------------------------------------------
 *      // Once an object is constructed, more interpolants that share the same
 *      // knots can be added.
 *      // Such interpolants can be evaluated simultaneously at a single place.
 *
 *      // build an object with 5 interpolants
 *      CubicSpline cubspl5(n, x, y);
 *      cubspl5.reserve(5); // reduce memory reallocations & data copies
 *      cubspl5.add(y2);
 *      cubspl5.add(y3, {CubicSpline::BoundaryType::first_deriv, 1.0}, {});
 *      cubspl5.add(y4, {}, {CubicSpline::BoundaryType::second_deriv, 2.0});
 *      cubspl5.add(y5);
 *
 *      // evaluates the five interpolants simultaneously at a single place
 *      cubspl5.multi_eval(x_interp, y_interp)
 *
 *      // evaluate the first and third interpolants at a single place
 *      std::vector<int> ind = {0, 2};
 *      cubspl5.multi_eval(ind.size(), ind.data(), x_interp, y_interp)
 *
 *      // evaluate the last interpolant (i_spline == 4) at multiple places
 *      cubspl5.eval(n_interp, x_interp, y_interp, nullptr, nullptr, 4)
 *
 *
 * Usage-2: static member functions
 *
 *      // step-1: computes the first-derivatives at knots
 *      // boundary conditions defaulted to "not-a-knot"
 *      CubicSpline::build(n, x, y, {}, {}, dy);
 *
 *      // Various boundary conditions and evenly spaced knots are supported
 *      // in the same way as the interpolant object.
 *
 *      // step-2: computes the interpolated values & derivatives
 *      CubicSpline::eval(n, x, y, dy, n_interp, x_interp, y_interp, dy_interp);
 *
 *      // Simultaneous evaluation of multiple interpolants are not supported
 *      // for static functions.
 *
 */
class CubicSpline
{    
    //*****************************************************************
    //                      boundary condition
    //*****************************************************************

public:

    /**
     * @brief Types of cubic spline boundary conditions.
     *
     * Supported types include:
     * - not_a_knot     The first or last two pieces are the same polynomial,
     *                  i.e., x[1] or x[n-2] is not a "knot". This does not
     *                  rely on any prior knowledge of the original function
     *                  and is the default option.
     * - first_deriv    user-defined first derivative
     * - second_deriv   user-defined second derivative
     * - periodic       the first and second derivatives at two ends are continuous.
     *                  This condition requires that y[0] = y[n-1] and it must be
     *                  applied to both ends
     */
    enum class BoundaryType
    {
        not_a_knot,
        first_deriv,
        second_deriv,
        periodic
    };


    /**
     * @brief Boundary condition for cubic spline interpolation.
     *
     * An object of this struct represents an actual boundary condition at one end,
     * which contains a type and possibly a value (first/second_deriv only).
     *
     */
    struct BoundaryCondition
    {
        // for not_a_knot and periodic
        BoundaryCondition(BoundaryType type = BoundaryType::not_a_knot);

        // for first/second_deriv
        BoundaryCondition(BoundaryType type, double val);

        BoundaryType type;
        double val = 0.0;
    };


    //*****************************************************************
    //                      interpolant object
    //*****************************************************************

public:

    CubicSpline()                   = delete;
    CubicSpline(CubicSpline const&) = default;
    CubicSpline(CubicSpline &&)     = default;

    CubicSpline& operator=(CubicSpline const&)  = default;
    CubicSpline& operator=(CubicSpline &&)      = default;

    ~CubicSpline() = default; 


    /**
     * @brief Builds an interpolant object.
     *
     * Constructing a cubic spline interpolant from a set of data points
     * (x[i], y[i]) (i=0,1,...,n-1) and boundary conditions.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points
     *                              ("knots", must be strictly increasing)
     * @param[in]   y               y coordinates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     *
     */
    CubicSpline(
        int n,
        const double* x,
        const double* y,
        const BoundaryCondition& bc_start = {},
        const BoundaryCondition& bc_end = {}
    );


    /**
     * @brief Builds an interpolant object with evenly-spaced knots.
     *
     * Constructing a cubic spline interpolant from a set of data points
     * (x0+i*dx, y[i]) (i=0,1,...,n-1) and boundary conditions.
     *
     * @param[in]   n               number of data points
     * @param[in]   x0              x coordinate of the first data point (first knot)
     * @param[in]   dx              spacing between knots (must be positive)
     * @param[in]   y               y coordinates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     *
     */
    CubicSpline(
        int n,
        double x0,
        double dx,
        const double* y,
        const BoundaryCondition& bc_start = {},
        const BoundaryCondition& bc_end = {}
    );


    /**
     * @brief Add an interpolant that shares the same knots.
     *
     * An object of this class can hold multiple interpolants with the same knots.
     * Once constructed, more interpolants sharing the same knots can be added by
     * this function. Multiple interpolants can be evaluated simultaneously at a
     * single place by multi_eval.
     *
     * @param[in]   y               y coordinates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     *
     */
    void add(
        const double* y,
        const BoundaryCondition& bc_start = {},
        const BoundaryCondition& bc_end = {}
    );


    /**
     * @brief Evaluates a single interpolant at multiple places.
     *
     * @param[in]   n_interp        number of places to evaluate the interpolant
     * @param[in]   x_interp        places where an interpolant is evaluated
     *                              (must be within the range of knots)
     * @param[out]  y_interp        values at x_interp
     * @param[out]  dy_interp       first derivatives at x_interp
     * @param[out]  d2y_interp      second derivatives at x_interp
     * @param[in]   i_spline        index of the interpolant to evaluate
     *
     * @note pass nullptr to any of the output would suppress the corresponding calculation
     *
     */
    void eval(
        int n_interp,
        const double* x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr,
        int i_spline = 0
    ) const;


    /**
     * @brief Evaluates multiple interpolants at a single place.
     *
     * @param[in]   n_spline        number of interpolants to evaluate
     * @param[in]   i_spline        indices of interpolants to evaluate
     * @param[in]   x_interp        place where interpolants are evaluated
     *                              (must be within the range of knots)
     * @param[out]  y_interp        values at x_interp
     * @param[out]  dy_interp       first derivatives at x_interp
     * @param[out]  d2y_interp      second derivatives at x_interp
     *
     * @note pass nullptr to any of the output would suppress the corresponding calculation
     *
     */
    void multi_eval(
        int n_spline,
        const int* i_spline,
        double x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr
    ) const;


    /**
     * @brief Evaluates all interpolants at a single place.
     *
     * @param[in]   x_interp        place where interpolants are evaluated
     *                              (must be within the range of knots)
     * @param[out]  y_interp        values at x_interp
     * @param[out]  dy_interp       first derivatives at x_interp
     * @param[out]  d2y_interp      second derivatives at x_interp
     *
     * @note pass nullptr to any of the output would suppress the corresponding calculation
     *
     */
    void multi_eval(
        double x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr
    ) const;


    /**
     * @brief Reserves memory for holding more interpolants.
     *
     * By default this class does not reserve memory for multiple interpolants.
     * Without reservation, whenever a new interpolant is added, memory has to
     * be reallocated and old data copied, which could waste a lot of time if
     * there's a large number of interpolants to add.
     *
     * This function help avoid repetitive memory reallocations and data copies
     * by a one-shot reservation.
     *
     * @param[in]   n_spline        expected total number of interpolants
     *
     */
    void reserve(int n_spline) { y_.reserve(n_spline * n_ * 2); }


    /// heap memory usage in bytes
    size_t heap_usage() const { return (x_.capacity() + y_.capacity()) * sizeof(double); }

    /// first knot
    double xmin() const { return xmin_; }

    /// last knot
    double xmax() const { return xmax_; }


private:

    /// number of cubic spline interpolants
    int n_spline_ = 0;

    /// number of knots
    int n_ = 0;

    /// first knot
    double xmin_ = 0.0;

    /// last knot
    double xmax_ = 0.0;

    /// spacing between knots (only used for evenly-spaced knots)
    double dx_ = 0.0;

    /// knots of the spline polynomial (remains empty for evenly-spaced knots)
    std::vector<double> x_;

    /// values and first derivatives at knots
    std::vector<double> y_;


    //*****************************************************************
    //                      static functions
    //*****************************************************************

public:

    /**
     * @brief Computes the first derivatives at knots for cubic spline
     * interpolation.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points
     *                              ("knots", must be strictly increasing)
     * @param[in]   y               y coordinates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     * @param[out]  dy              first derivatives at knots
     *
     */
    static void build(
        int n,
        const double* x,
        const double* y,
        const BoundaryCondition& bc_start,
        const BoundaryCondition& bc_end,
        double* dy
    );


    /**
     * @brief Computes the first derivatives at evenly-spaced knots for
     * cubic spline interpolation.
     *
     * @param[in]   n               number of data points
     * @param[in]   dx              spacing between knots (must be positive)
     * @param[in]   y               y coordinates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     * @param[out]  dy              first derivatives at knots
     *
     */
    static void build(
        int n,
        double dx,
        const double* y,
        const BoundaryCondition& bc_start,
        const BoundaryCondition& bc_end,
        double* dy
    );


    /**
     * @brief Evaluates a cubic spline polynomial at multiple places.
     *
     * @param[in]   n               number of knots
     * @param[in]   x               knots (must be strictly increasing)
     * @param[in]   y               values at knots
     * @param[in]   dy              first derivatives at knots
     * @param[in]   n_interp        number of places to evaluate the interpolant
     * @param[in]   x_interp        places where the interpolant is evaluated
     *                              (must be within the range of knots)
     * @param[out]  y_interp        values at x_interp
     * @param[out]  dy_interp       first derivatives at x_interp
     * @param[out]  d2y_interp      second derivatives at x_interp
     *
     * @note pass nullptr to any of the output would suppress the corresponding calculation
     *
     */
    static void eval(
        int n,
        const double* x,
        const double* y,
        const double* dy,
        int n_interp,
        const double* x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr
    );


    /**
     * @brief Evaluates a cubic spline polynomial with evenly spaced knots.
     *
     * @param[in]   n               number of knots
     * @param[in]   x0              first knot
     * @param[in]   dx              spacing between knots
     * @param[in]   y               values at knots
     * @param[in]   dy              first derivatives at knots
     * @param[in]   n_interp        number of places to evaluate the interpolant
     * @param[in]   x_interp        places where the interpolant is evaluated
     *                              (must be within the range of knots)
     * @param[out]  y_interp        values at x_interp
     * @param[out]  dy_interp       first derivatives at x_interp
     * @param[out]  d2y_interp      second derivatives at x_interp
     *
     * @note pass nullptr to any of the output would suppress the corresponding calculation
     *
     */
    static void eval(
        int n,
        double x0,
        double dx,
        const double* y,
        const double* dy,
        int n_interp,
        const double* x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr
    );


private:

    /// Computational routine for building cubic spline interpolant
    static void _build(
        int n,
        const double* dx,
        const double* y,
        const BoundaryCondition& bc_start,
        const BoundaryCondition& bc_end,
        double* dy
    );


    /**
     * @brief Segment index lookup.
     *
     * Given a strictly increasing array x and a target within the range of
     * x, this function returns an index i such that x[i] <= target < x[i+1]
     * if target != x[n-1], or n-2 if t == x[n-1].
     *
     */
    static inline int _index(int n, const double* x, double target);


    /// Segment index lookup (evenly spaced knots).
    static inline int _index(int n, double x0, double dx, double target);


    /// Evaluates a batch of cubic polynomials.
    static inline void _cubic(
        int n,
        const double* w,
        const double* c0,
        const double* c1,
        const double* c2,
        const double* c3,
        double* y,
        double* dy,
        double* d2y
    );


    /// Asserts that the input arguments are valid for constructing a cubic spline.
    static void _validate_build(
        int n,
        const double* dx,
        const double* y,
        const BoundaryCondition& bc_start,
        const BoundaryCondition& bc_end
    );


    /// Asserts that the input arguments are valid for interpolating a cubic spline.
    static void _validate_eval(
        int n,
        const double (&u)[2],
        const double* x,
        const double* y,
        const double* dy,
        int n_interp,
        const double* x_interp
    );


    /**
     * @brief Solves a cyclic tridiagonal linear system.
     *
     * A cyclic tridiagonal linear system A*x=b where b is a vector and
     *
     *        --                                             --   
     *        |  d[0]   u[0]                           l[n-1] |
     *        |  l[0]   d[1]   u[1]                           |
     *   A =  |         l[1]   d[2]   u[2]                    |
     *        |                ...    ...      ...            |
     *        |                      l[n-3]   d[n-2]   u[n-2] |
     *        |  u[n-1]                       l[n-2]   d[n-1] |
     *        --                                             --
     *
     * is transformed to a tridiagonal linear system by the Sherman-Morrison
     * formula, and then solved by dgtsv.
     *
     * @param[in]       n       size of the linear system
     * @param[in]       d       main diagonal
     * @param[in]       u       superdiagonal
     * @param[in]       l       subdiagonal
     * @param[in,out]   b       right hand side of the linear system; will be
     *                          overwritten by the solution on finish.
     *
     * @note d, l, u are all overwritten in this function.
     *
     */
    static void _solve_cyctri(int n, double* d, double* u, double* l, double* b);
};

} // namespace ModuleBase

#endif

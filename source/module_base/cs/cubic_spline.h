#ifndef CUBIC_SPLINE_INTERPOLATION_H_
#define CUBIC_SPLINE_INTERPOLATION_H_

#include <vector>

/**
 * @brief Cubic spline interplation.
 *
 * Interpolating a given set of data points (x[i], y[i]) (i=0,...,n-1) by piecewise
 * cubic polynomials with continuous first and second derivatives:
 *
 *      p_i(x) = c0[i] + c1[i]*(x-x[i]) + c2[i]*(x-x[i])^2 + c3[i]*(x-x[i])^3
 *
 * where p_i(x) is defined on [x[i], x[i+1]].
 *
 * There are two ways to use this class. The first way uses the class as an
 * interpolant object; the second way uses static member functions.
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
 *      cubspl.eval(n_interp, x_interp, y_interp, dy_interp);
 *
 *      // evaluates the derivative only
 *      cubspl.eval(n_interp, x_interp, nullptr, dy_interp);
 *
 *      // Build an interpolant with evenly spaced x (x[i] = x0 + i*dx)
 *      CubicSpline cubspl2(n, x0, dx, y);
 *
 *      // "not-a-knot" boundary condition is used by default at both ends.
 *      // Other supported boundary conditions include first/second derivatives and
 *      // periodic boundary condition.
 *
 *      // f''(start) = 1.0 and f'(end) = 3.0
 *      CubicSpline cubspl3(n, x0, dx, y,
 *                          {CubicSpline::BoundaryType::second_deriv, 1.0},
 *                          {CubicSpline::BoundaryType::first_deriv, 3.0});
 *
 *      // applying periodic boundary condition requires that
 *      // 1. it be specified at both ends;
 *      // 2. y[0] must equal y[n-1]
 *      CubicSpline cubspl4(n, x0, dx, y,
 *                          {CubicSpline::BoundaryType::periodic},
 *                          {CubicSpline::BoundaryType::periodic});
 *
 *      // interpolant objects support holding multiple interpolants with the same knots
 *      cubspl.add(y2);
 *      cubspl.add(y3);
 *
 *      // such multiple interpolants can be evaluated simultaneously at a single place
 *
 *      // evaluate all interpolants and their first derivatives at x_interp
 *      cubspl.multi_eval(x_interp, y, dy)
 *
 *      // evaluate the first and third interpolants at x_interp
 *      std::vector<int> ind = {0, 2};
 *      cubspl.multi_eval(ind.size(), ind.data(), x_interp, y)
 *
 *
 * Usage-2: static member functions
 *
 *      // step-1: computes the first-derivatives at knots
 *      // PS: boundary conditions defaulted to "not-a-knot"
 *      CubicSpline::build(n, x, y, {}, {}, dy);
 *
 *      // Various boundary conditions and explicitly specified evenly spaced
 *      // knots are supported in the same way as the interpolant object.
 *
 *      // step-2: interpolates with knots, values & derivatives
 *      CubicSpline::eval(n, x, y, dy, n_interp, x_interp, y_interp, dy_interp);
 *
 *      // Interpolating multiple interpolants are not supported for static functions.
 *
 */
class CubicSpline
{

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


    /**
     * @brief Builds the interpolant object.
     *
     * Constructing a cubic spline interpolant from the given set of data points
     * (x[i], y[i]) (i=0,1,...,n-1) and boundary conditions.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points, must be strictly increasing
     * @param[in]   y               y coordiates of data points
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
     * @brief Builds the interpolant object with evenly-spaced knots.
     *
     * Constructing a cubic spline interpolant from the given set of data points
     * (x0+i*dx, y[i]) (i=0,1,...,n-1) and boundary conditions.
     *
     * @param[in]   n               number of data points
     * @param[in]   x0              starting x coordinate
     * @param[in]   dx              spacing between knots
     * @param[in]   y               y coordiates of data points
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
     * @brief Adds more interpolants of the same knots to an existing object.
     *
     * An object of this class supports holding multiple interpolants with the same knots.
     * Once the first interpolant is contructed by the constructor, more interpolants can
     * be added by this function.
     *
     * @param[in]   y               y coordiates of data points
     * @param[in]   bc_start        boundary condition at start
     * @param[in]   bc_end          boundary condition at end
     *
     * @note whether the new interpolant stores full coefficients or not is determined by
     * the constructor.
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
     * @param[in]   x_interp        places where the interpolant is evaluated; must be
     *                              within [x_[0], x_[n-1]] or [x0_, x0_+(n-1)*dx_]
     * @param[out]  y_interp        interpolated values
     * @param[out]  dy_interp       first derivatives at x
     * @param[out]  d2y_interp      second derivatives at x
     * @param[in]   ind             index of the interpolant to evaluate
     *
     * @note the corresponding calculation is skipped if either y, dy or d2y is nullptr.
     *
     */
    void eval(
        int n_interp,
        const double* x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr,
        int ind = 0
    ) const;


    /**
     * @brief Evaluates selected multiple interpolants at a single place.
     *
     * @param[in]   n_spline        number of interpolants to evaluate
     * @param[in]   ind             indices of interpolants to evaluate
     * @param[in]   x_interp        place where interpolants are evaluated; must be
     *                              within [x_[0], x_[n-1]] or [x0_, x0_+(n-1)*dx_]
     * @param[out]  y_interp        interpolated values
     * @param[out]  dy_interp       first derivatives at x
     * @param[out]  d2y_interp      second derivatives at x
     *
     * @note the corresponding calculation is skipped if either y, dy or d2y is nullptr.
     *
     */
    void multi_eval(
        int n_spline,
        const int* ind,
        double x_interp,
        double* y_interp,
        double* dy_interp = nullptr,
        double* d2y_interp = nullptr
    ) const;


    /**
     * @brief Evaluates all interpolants at a single place.
     *
     * @param[in]   x               place where interpolants are evaluated; must be
     *                              within [x_[0], x_[n-1]] or [x0_, x0_+(n-1)*dx_]
     * @param[out]  y               interpolated values
     * @param[out]  dy              first derivatives at x
     * @param[out]  d2y             second derivatives at x
     *
     * @note the corresponding calculation is skipped if either y, dy or d2y is nullptr.
     *
     */
    void multi_eval(
        double x,
        double* y,
        double* dy = nullptr,
        double* d2y = nullptr
    ) const;


    /**
     * @brief Reserves memory for multiple interpolants.
     *
     * This class does not reserve memory for multiple interpolants by default.
     * Memory is reallocated and old data copied whenever a new interpolant is added,
     * which could be slow if the number of interpolants to be added is large.
     *
     * Memory reallocations and data copy can be avoided by reserving a sufficient
     * amount of memory in advance, which can be done by this function.
     *
     * @param[in]   n               expected total number of interpolants
     *
     */
    void reserve(int n) { c_.reserve(n * n_ * 2); }


    /// heap memory usage in bytes
    size_t heap_usage() const { return (x_.capacity() + y_.capacity()) * sizeof(double); }


private:

    /// number of knots
    int n_ = 0;

    /// starting x coordinate
    double x0_ = 0.0;

    /// spacing between knots (used for evenly-space knots only)
    double dx_ = 0.0;

    /// knots of the spline polynomial (remains empty for evenly-spaced knots)
    std::vector<double> x_;

    /// values and first derivatives at knots (all values come first)
    std::vector<double> y_;


    //*****************************************************************
    //                      static functions
    //*****************************************************************

public:

    /**
     * @brief Computes the coefficients of piecewise cubic polynomials in a cubic spline.
     *
     * @param[in]   n               number of data points
     * @param[in]   x               x coordinates of data points ("knots"),
     *                              must be strictly increasing
     * @param[in]   y               y coordiates of data points
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
     * @brief Computes the coefficients of piecewise cubic polynomials in a cubic spline
     * with evenly spaced knots.
     *
     * @param[in]   n               number of data points
     * @param[in]   dx              spacing between knots
     * @param[in]   y               y coordiates of data points
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
     * @brief Evaluates a cubic spline polynomial.
     *
     * @param[in]   n               number of knots
     * @param[in]   x               knots
     * @param[in]   y               y coordiates of data points
     * @param[in]   dy              first derivatives at knots
     * @param[in]   n_interp        number of points to evaluate the spline polynomial
     * @param[in]   x_interp        places where the spline polynomial is evaluated,
     *                              must be within [x[0], x[n-1]]
     * @param[out]  y_interp        interpolated values at x_interp
     * @param[out]  dy_interp       interpolated first derivatives at x_interp
     * @param[out]  d2y_interp      interpolated second derivatives at x_interp
     *
     * @note pass nullptr to any of the y's would suppress the corresponding calculation
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
     * @param[in]   x0              starting x coordinate
     * @param[in]   dx              spacing between knots
     * @param[in]   y               values at knots
     * @param[in]   dy              first derivatives at knots
     * @param[in]   n_interp        number of points to evaluate the spline polynomial
     * @param[in]   x_interp        places where the spline polynomial is evaluated,
     *                              must be within [x[0], x[n-1]]
     * @param[out]  y_interp        interpolated values at x_interp
     * @param[out]  dy_interp       interpolated first derivatives at x_interp
     * @param[out]  d2y_interp      interpolated second derivatives at x_interp
     *
     * @note pass nullptr to any of the y's would suppress the corresponding calculation
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
     * @brief Index of the polynomial piece containing x.
     *
     * Given strictly increasing knots and knots[0] <= x <= knots[n-1],
     * this function returns an index i such that knots[i] <= x < knots[i+1]
     * if x != knots[n-1], or n-2 if x == knots[n-1].
     *
     * @param[in]   n           number of knots
     * @param[in]   knots       x coordinates of interpolant's data points
     * @param[in]   x           place where the interpolant is evaluated
     *
     * @note No sanity check is performed in this function
     *
     */
    static int _index(int n, const double* knots, double x)
    {
        int i = (std::upper_bound(knots, knots + n, x) - knots) - 1;
        return i - (i == n - 1);
    }

    /// Index of the polynomial segment containing x (evenly spaced knots)
    static int _index(int n, double x0, double dx, double x)
    {
        int i = (x - x0) / dx;
        return i - (i == n - 1);
    }


    /// Evaluates a batch of cubic polynomials.
    static inline void _cubic_eval(
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
        const double* x,
        double x0,
        double dx,
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
     *        |  d[0]   u[0]                            l[0]  |
     *        |  l[1]   d[1]   u[1]                           |
     *   A =  |         l[2]   d[2]   u[2]                    |
     *        |                ...    ...      ...            |
     *        |                      l[n-2]   d[n-2]   u[n-2] |
     *        |  u[n-1]                       l[n-1]   d[n-1] |
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

#endif

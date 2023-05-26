#ifndef CUBIC_SPLINE_INTERPOLATOR_H_
#define CUBIC_SPLINE_INTERPOLATOR_H_

#include "lapack_connector.h"

namespace ModuleBase
{

class CubicSpline
{
  public:
    CubicSpline() {};
    CubicSpline(CubicSpline const&) = delete;
    CubicSpline& operator=(CubicSpline const&) = delete;

    ~CubicSpline();

    enum class BoundaryCondition { not_a_knot, first_deriv, second_deriv, periodic };

    void build(const int n, 
            const double* const x, 
            const double* const y, 
            BoundaryCondition bc_start = BoundaryCondition::not_a_knot, 
            BoundaryCondition bc_end = BoundaryCondition::not_a_knot,
            const double deriv_start = 0.0, 
            const double deriv_end = 0.0
    );

    void interpolate(const int n, const double* const x, double* const y);

  private:
    int n_ = 0;
    double* x_ = nullptr; // knots

    // polynomial coefficients
    // The i-th piece polynomial is given by (i=0,1,...)
    // P[i](x) = c0[i] + c1[i]*(x-x[i]) + c2[i]*(x-x[i])^2 + c3[i]*(x-x[i])^3
    double* c0_ = nullptr;
    double* c1_ = nullptr;
    double* c2_ = nullptr;
    double* c3_ = nullptr;

    bool is_uniform_ = false;
    double uniform_thr_ = 1e-15;

    void sanity_check(const int n, const double* const x, const double* const y, BoundaryCondition bc_start, BoundaryCondition bc_end);

    void solve_cyctri(int n, double* const diag, double* const supdiag, double* const subdiag, double* const b);

    void cleanup();
};

}; // namespace ModuleBase

#endif

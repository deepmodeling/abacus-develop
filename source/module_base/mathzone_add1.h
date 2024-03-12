#ifndef MATHZONE_ADD1_H
#define MATHZONE_ADD1_H

#include <cassert>
#include <cmath>

namespace ModuleBase
{

/************************************************************************
LiaoChen add @ 2010/03/09 to add efficient functions in LCAO calculation
************************************************************************/
// Only used in module_basis/module_ao
class Mathzone_Add1
{
  public:
    Mathzone_Add1();
    ~Mathzone_Add1();

    static double dualfac(const int& l);
    static double factorial(const int& l);
    /**
     * @brief calculate second derivatives for cubic
     * spline interpolation
     *
     * @param[in] rad   x before interpolation
     * @param[in] rad_f f(x) before interpolation
     * @param[in] mesh  number of x before interpolation
     * @param[in] yp1   f'(0) boundary condition
     * @param[in] ypn   f'(n) boundary condition
     * @param[out] y2   f''(x)
     */
    static void SplineD2(const double* rad,
                         const double* rad_f,
                         const int& mesh,
                         const double& yp1,
                         const double& ypn,
                         double* y2);

    /**
     * @brief cubic spline interpolation
     *
     * @param[in] rad   x before interpolation
     * @param[in] rad_f f(x) before inpterpolation
     * @param[in] y2    f''(x) before interpolation
     * @param[in] mesh  number of x before interpolation
     * @param[in] r     x after interpolation
     * @param[in] rsize number of x after interpolation
     * @param[out] y    f(x) after interpolation
     * @param[out] dy   f'(x) after interpolation
     * @param[out] ddy   f''(x) after interpolation
     */
    static void Cubic_Spline_Interpolation(const double* const rad,
                                           const double* const rad_f,
                                           const double* const y2,
                                           const int& mesh,
                                           const double* const r,
                                           const int& rsize,
                                           double* const y,
                                           double* const dy);

    /**
     * @brief "spline like interpolation" of a uniform
     * funcation of r
     *
     * @param[in] rad_f f(x) before interpolation
     * @param[in] msh   number of x known
     * @param[in] dr    uniform distance of x
     * @param R         f(R) is to be calculated
     * @return double   f(R)
     */
    static double Uni_RadialF(const double* rad_f, const int& msh, const double& dr, const double& R);

    static void Uni_Deriv_Phi(const double* radf, const int& mesh, const double& dr, const int& nd, double* phind);

  private:
    const static int sph_lmax = 20;
    static double** c_ln_c;
    static double** c_ln_s;
};

} // namespace ModuleBase

#endif

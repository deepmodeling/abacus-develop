#ifndef GRID_RADIAL_H
#define GRID_RADIAL_H

#include <vector>

namespace Grid {
namespace Radial {

/**
 * @brief Radial quadratures.
 *
 * This namespace contains functions that generate grids and weights
 * for numerical integration
 *
 *      / inf     2
 *      |     dr r  g(r) ~ \sum_i w[i] g(r[i])
 *      /  0
 *
 */

/**
 * Baker, J., Andzelm, J., Scheiner, A., & Delley, B. (1994).
 * The effect of grid quality and weight derivatives in
 * density functional calculations.
 * The Journal of chemical physics, 101(10), 8894-8902.
 *
 * Zhang, I. Y., Ren, X., Rinke, P., Blum, V., & Scheffler, M. (2013).
 * Numeric atom-centered-orbital basis sets with valence-correlation
 * consistency from H to Ar.
 * New Journal of Physics, 15(12), 123033.
 *
 * @note nbase is the number of points of the "base" grid, i.e.,
 * before applying the "radial multiplier" introduced by Zhang et al.
 * The true number of grid points is (nbase+1) * mult - 1.
 */
void baker(int nbase, double R, double* r, double* w, int mult = 1);
void baker(int nbase, double R, std::vector<double>& r,
           std::vector<double>& w, int mult = 1);


/**
 * Murray, C. W., Handy, N. C., & Laming, G. J. (1993).
 * Quadrature schemes for integrals of density functional theory.
 * Molecular Physics, 78(4), 997-1014.
 */
void murray(int n, double R, double* r, double* w);


/**
 * Treutler, O., & Ahlrichs, R. (1995).
 * Efficient molecular numerical integration schemes.
 * The Journal of Chemical Physics, 102(1), 346-354.
 *
 * @note The M4 grid reduces to the M3 grid when alpha = 0.
 */
void treutler_m4(int n, double R, double* r, double* w, double alpha = 0.6);


/**
 * Mura, M. E., & Knowles, P. J. (1996).
 * Improved radial grids for quadrature in molecular
 * density‐functional calculations.
 * The Journal of chemical physics, 104(24), 9848-9858.
 */
void mura(int n, double R, double* r, double* w);

} // end of namespace Radial
} // end of namespace Grid

#endif

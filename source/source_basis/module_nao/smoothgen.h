#ifndef SMOOTHGEN_H_
#define SMOOTHGEN_H_

#include <vector>

/**
 * @brief Smoothly truncates a radial function at a cutoff radius.
 *
 * Starting with a radial function chi on grid r, this function generates
 *
 *     alpha(r) = N * chi(r) * (1 - exp(-(r - rcut)^2 / (2 sigma^2)))
 *
 * on the truncated grid, i.e. the first lower_bound(rcut) + 1 points of r
 * (all points below rcut plus the first point at or beyond rcut), so that
 * alpha vanishes at rcut. The factor N normalizes alpha such that
 * \int r^2 alpha(r)^2 dr = 1. The smoothness parameter sigma is optimized
 * by a bracketed bisection search to minimize the squared radial derivative
 * difference \int r^2 (alpha'(r) - chi'(r))^2 dr, with derivatives evaluated
 * by cubic spline interpolation.
 *
 * Requires r to be sorted in ascending order and rcut < r[nr - 1].
 *
 * @param[in]   nr      number of grid points
 * @param[in]   r       radial grid
 * @param[in]   chi     radial function
 * @param[in]   rcut    cutoff radius
 * @param[out]  alpha   smoothed radial function on the truncated grid,
 *                      resized to lower_bound(rcut) + 1 points
 */
void smoothgen(
    const int nr,
    const double* r,
    const double* chi,
    const double rcut,
    std::vector<double>& alpha);

#endif

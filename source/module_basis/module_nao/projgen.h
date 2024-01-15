#ifndef MODULE_NAO_PROJGEN_H_
#define MODULE_NAO_PROJGEN_H_

#include <vector>

/**
 * @brief Generates a projector's radial function for DFT+U
 *
 * Starting with a numerical radial function chi on grid r with angular momentum l,
 * given a smaller cutoff radius rcut and the number of spherical Bessel components (j_l),
 * this function generates a new radial function alpha of the same l on the truncated grid
 * which satisfies the following conditions:
 *
 * * alpha = \sum_p j_l(theta[p]*r/rcut) * c[p] where theta[p] is the p-th zero of j_l;
 * * \int_0^rcut alpha(r) r^2 dr = 1 (normalization);
 * * \int_0^rcut alpha(r) chi(r) r^2 dr is maximized;
 *
 * @param[in]   l       angular momentum
 * @param[in]   nr      number of grid points
 * @param[in]   r       radial grid
 * @param[in]   chi     radial function
 * @param[in]   rcut    cutoff radius of the projector
 * @param[in]   nbes    number of spherical Bessel components
 * @param[out]  alpha   new radial function of the projector
 *
 */
void projgen(
    const int l, 
    const int nr, 
    const double* r, 
    const double* chi, 
    const double rcut, 
    const int nbes, 
    std::vector<double>& alpha);

#endif
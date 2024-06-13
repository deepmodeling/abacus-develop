#ifndef NUMERICAL_BASIS_JYJY_H
#define NUMERICAL_BASIS_JYJY_H

#include "module_base/complexarray.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"

#include <tuple>
#include <vector>

namespace NumericalBasis
{
std::vector<std::tuple<int, int, int, int>> indexgen(const std::vector<int>& natom, const std::vector<int>& lmax);

/**
 * @brief <jy|op|jy> overlap matrix (two-center integration)
 *
 *
 * @param[in]   type                'S' (op = 1) or 'T' (kinetic, op = -\nabla^2)
 * @param[in]   lmax                maximum angular momentum
 * @param[in]   nbes                number of Bessel functions
 * @param[in]   rcut                cutoff radius
 * @param[in]   tau_cart            atomic positions (in Bohr)
 * @param[in]   latvec              lattice vectors (in Bohr)
 * @param[in]   mu_index            composite index
 *
 */
ModuleBase::ComplexArray cal_overlap_Sq(const char type, // 'S' or 'T'
                                        const int lmax, const int nbes, const double rcut,
                                        const std::vector<std::vector<ModuleBase::Vector3<double>>>& tau_cart,
                                        const ModuleBase::Matrix3& latvec,
                                        const std::vector<std::tuple<int, int, int, int>>& mu_index);

/**
 * @brief Searching for all relative position vectors for periodic images
 * within a cutoff radius.
 *
 * Given an initial relative position vector d0 and a searching radius r,
 * this function returns all d such that
 *
 *      d = d0 + n0*a0 + n1*a1 + n2*a2   and |d| < r
 *
 * where n0, n1, n2 are integers and a0, a1, a2 are the lattice vectors.
 *
 */
std::vector<ModuleBase::Vector3<double>> neighbor_vec(const ModuleBase::Vector3<double>& d0,
                                                      const ModuleBase::Matrix3& latvec, const double r);

} // namespace NumericalBasis

#endif

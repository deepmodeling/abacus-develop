#ifndef NUMERICAL_BASIS_JYJY_H
#define NUMERICAL_BASIS_JYJY_H

#include <tuple>
#include <vector>
#include "module_base/vector3.h"
#include "module_base/complexarray.h"

namespace NumericalBasis
{
std::vector<std::tuple<int, int, int, int>> indexgen(const std::vector<int>& natom,
                                                     const std::vector<int>& lmax);

/**
 * @brief <jy|op|jy> overlap matrix (two-center integration)
 *
 *
 * @param[in]   type                'S' (op = 1) or 'T' (kinetic, op = -\nabla^2)
 * @param[in]   lmax                maximum angular momentum
 * @param[in]   nbes                number of Bessel functions
 * @param[in]   rcut                cutoff radius
 * @param[in]   sigma               smoothing parameter
 * @param[in]   tau_cart            atomic positions (in Bohr)
 * @param[in]   mu_index            composite index
 *
 */
ModuleBase::ComplexArray cal_overlap_Sq(
    const char type, // 'S' or 'T'
    const int lmax,
    const int nbes,
    const double rcut,
    const std::vector<std::vector<ModuleBase::Vector3<double>>>& tau_cart,
    const std::vector<std::tuple<int, int, int, int>>& mu_index
);

}

#endif

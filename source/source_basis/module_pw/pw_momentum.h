#ifndef MODULE_PW_MOMENTUM_H
#define MODULE_PW_MOMENTUM_H

#include "source_basis/module_pw/pw_basis_k.h"

namespace ModulePW
{
/** @brief Return the kinetic correction 2*p*A + A^2 in Ry for a Hartree vector potential.
 *  @note Evaluate geometry in double precision; callers convert the result to their native precision.
 */
inline double kinetic_shift(const PW_Basis_K& basis, const int ik, const int ig,
                            const ModuleBase::Vector3<double>& A_ha)
{
    const ModuleBase::Vector3<double> p = basis.getgpluskcar(ik, ig) * basis.tpiba;
    return 2.0 * (p.x * A_ha.x + p.y * A_ha.y + p.z * A_ha.z) + A_ha.norm2();
}

/** @brief Return the shifted kinetic energy in Ry, retaining the cached unshifted squared momentum. */
inline double shifted_kinetic(const PW_Basis_K& basis, const int ik, const int ig,
                              const ModuleBase::Vector3<double>& A_ha)
{
    return basis.getgk2(ik, ig) * basis.tpiba2 + kinetic_shift(basis, ik, ig, A_ha);
}
} // namespace ModulePW

#endif

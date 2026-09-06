/// @file dftu_nao_str_r.h
/// @brief DFT+U stress calculation in real space (r-space)
///
/// This file provides the real-space implementation of DFT+U stress contribution
/// from a single atom pair (I,J,R). It is independent of k-point sampling because
/// the real-space density matrix (DMR) already contains the Brillouin-zone integration.
///
/// Naming convention: _r suffix denotes real-space implementation,
/// corresponding to _k suffix for k-space (legacy) implementation.
///
/// The stress formula for atom pair (I,J,R) is:
///
///   sigma_{alpha,beta} += -(1/Omega) * sum_{m,m'} V_U_{mm'}(I) * DMR_{mu,nu}(J1,J2,R) * [
///       d<phi_{mu,0}|chi_m(I)>/d tau_{J1,alpha} * <chi_m'(I)|phi_{nu,R}> * R_{J1,beta}
///     + <phi_{mu,0}|chi_m(I)> * d<chi_m'(I)|phi_{nu,R}>/d tau_{J2,alpha} * R_{J2,beta}
///   ]
///
/// where R_{J1} and R_{J2} are the position vectors of atoms J1 and J2 relative to
/// the on-site atom I, and Omega is the unit cell volume.

#ifndef DFTU_NAO_STR_R_H
#define DFTU_NAO_STR_R_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/vector3.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <unordered_map>
#include <vector>

namespace hamilt
{

// Forward declarations to avoid circular dependency with dftu_lcao_op.h
template <typename TK, typename TR>
class OperatorLCAO;

template <typename T>
class DFTU;

/**
 * @brief Compute DFT+U stress contribution from a single atom pair (I,J,R) in real space
 *
 * The stress tensor component (alpha,beta) receives contributions from both the derivative
 * of <phi|chi> with respect to the strain and the explicit R_beta factor:
 *
 *   sigma_{alpha,beta} += -(1/Omega) * sum_{m,m'} V_U_{mm'}(I) * DMR_{mu,nu}(J1,J2,R) * [
 *       d<phi_{mu,0}|chi_m(I)>/d tau_{J1,alpha} * <chi_m'(I)|phi_{nu,R}> * R_{J1,beta}
 *     + <phi_{mu,0}|chi_m(I)> * d<chi_m'(I)|phi_{nu,R}>/d tau_{J2,alpha} * R_{J2,beta}
 *   ]
 *
 * The strain derivative of the two-center integral is approximated by:
 *   d<phi|chi>/d epsilon_{alpha,beta} ~ d<phi|chi>/d tau_alpha * R_beta
 *
 * @note The stress is accumulated in Voigt notation (6 components: xx, xy, xz, yy, yz, zz)
 *       and symmetrized at the caller level. The final scaling by lat0/Omega is applied
 *       after MPI reduction.
 *
 * @param iat1        [in] global atom index of J1
 * @param iat2        [in] global atom index of J2
 * @param pv          [in] Parallel_Orbitals for basis index mapping
 * @param nlm1_all    [in] pre-computed <phi|chi> and derivatives for atom J1
 * @param nlm2_all    [in] pre-computed <phi|chi> and derivatives for atom J2
 * @param pot_onsite  [in] flattened V_U matrix
 * @param dmR_pointer [in] pointer to DMR matrix blocks for each spin
 * @param nspin       [in] number of spin channels
 * @param dis1        [in] position vector of J1 relative to I
 * @param dis2        [in] position vector of J2 relative to I
 * @param stress      [out] stress accumulator (6 components in Voigt notation)
 */
template <typename TK, typename TR>
void cal_str_IJR_nao_r(const DFTU<OperatorLCAO<TK, TR>>* dftu_op,
                      const int& iat1,
                      const int& iat2,
                      const Parallel_Orbitals* pv,
                      const std::unordered_map<int, std::vector<double>>& nlm1_all,
                      const std::unordered_map<int, std::vector<double>>& nlm2_all,
                      const std::vector<double>& pot_onsite_in,
                      const hamilt::BaseMatrix<double>** dmR_pointer,
                      const int nspin,
                      const ModuleBase::Vector3<double>& dis1,
                      const ModuleBase::Vector3<double>& dis2,
                      double* stress);

} // namespace hamilt

#endif // DFTU_NAO_STR_R_H

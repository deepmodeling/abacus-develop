/// @file dftu_nao_for_r.h
/// @brief DFT+U force calculation in real space (r-space)
///
/// This file provides the real-space implementation of DFT+U force contribution
/// from a single atom pair (I,J,R). It is independent of k-point sampling because
/// the real-space density matrix (DMR) already contains the Brillouin-zone integration.
///
/// Naming convention: _r suffix denotes real-space implementation,
/// corresponding to _k suffix for k-space (legacy) implementation.
///
/// The force formula for atom pair (I,J,R) is:
///
///   F_{J1} += sum_{m,m'} V_U_{mm'}(I) * <phi_{mu,0}|chi_m(I)>
///             * d<chi_m'(I)|phi_{nu,R}>/d tau_{J1} * DMR_{mu,nu}(J1,J2,R)
///
///   F_{J2} -= sum_{m,m'} V_U_{mm'}(I) * d<phi_{mu,0}|chi_m(I)>/d tau_{J2}
///             * <chi_m'(I)|phi_{nu,R}> * DMR_{mu,nu}(J1,J2,R)
///
/// where V_U_{mm'}(I) = U_eff * (delta_{mm'}/2 - n_{m'm}(I)) is the on-site
/// Hubbard potential, and the two-center integrals <phi|chi> are pre-computed
/// by TwoCenterIntegrator::snap().

#ifndef DFTU_NAO_FOR_R_H
#define DFTU_NAO_FOR_R_H

#include "source_basis/module_ao/parallel_orbitals.h"
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
 * @brief Compute DFT+U force contribution from a single atom pair (I,J,R) in real space
 *
 * For a given on-site atom I with correlated orbital l, and a pair of basis atoms (J1, J2)
 * separated by lattice vector R, the force contribution is:
 *
 *   F_{J1} += sum_{m,m'} V_U_{mm'}(I) * <phi_{mu,0}|chi_m(I)>
 *             * d<chi_m'(I)|phi_{nu,R}>/d tau_{J1} * DMR_{mu,nu}(J1,J2,R)
 *
 *   F_{J2} -= sum_{m,m'} V_U_{mm'}(I) * d<phi_{mu,0}|chi_m(I)>/d tau_{J2}
 *             * <chi_m'(I)|phi_{nu,R}> * DMR_{mu,nu}(J1,J2,R)
 *
 * where mu runs over orbitals on J1, nu over orbitals on J2, and m,m' over the 2l+1
 * magnetic quantum numbers of the correlated orbital.
 *
 * The two-center integrals <phi|chi> and their derivatives are pre-computed by
 * TwoCenterIntegrator::snap() and stored in nlm_tot.
 *
 * @param iat1        [in] global atom index of J1
 * @param iat2        [in] global atom index of J2
 * @param pv          [in] Parallel_Orbitals for basis index mapping
 * @param nlm1_all    [in] pre-computed <phi|chi> and derivatives for atom J1
 * @param nlm2_all    [in] pre-computed <phi|chi> and derivatives for atom J2
 * @param pot_onsite  [in] flattened V_U matrix: [m1*m_size + m2 + is*m_size2]
 * @param dmR_pointer [in] pointer to DMR matrix blocks for each spin
 * @param nspin       [in] number of spin channels (1, 2, or 4)
 * @param force1      [out] force accumulator for atom J1 (3 components)
 * @param force2      [out] force accumulator for atom J2 (3 components)
 *
 * @note For nspin=1, the force is scaled by 2.0 at the caller level to account
 *       for spin degeneracy. For nspin=2, spin-up and spin-down are summed explicitly.
 *       For nspin=4 (non-collinear), the spinor structure is handled via npol=2 indexing.
 */
template <typename TK, typename TR>
void cal_for_IJR_nao_r(const DFTU<OperatorLCAO<TK, TR>>* dftu_op,
                     const int& iat1,
                     const int& iat2,
                     const Parallel_Orbitals* pv,
                     const std::unordered_map<int, std::vector<double>>& nlm1_all,
                     const std::unordered_map<int, std::vector<double>>& nlm2_all,
                     const std::vector<double>& pot_onsite_in,
                     const hamilt::BaseMatrix<double>** dmR_pointer,
                     const int nspin,
                     double* force1,
                     double* force2);

} // namespace hamilt

#endif // DFTU_NAO_FOR_R_H

/// @file dftu_nao_fs_r.h
/// @brief DFT+U force and stress unified entry in real space (r-space)
///
/// This file provides the unified entry for DFT+U force/stress using the
/// real-space density matrix (DMR). It is independent of k-point sampling
/// because DMR already contains the Brillouin-zone integration.
///
/// Naming convention: _r suffix denotes real-space implementation,
/// corresponding to _k suffix for k-space (legacy) implementation.
///
/// The DFT+U force on atom J is derived from the Hubbard correction energy:
///
///   E_U = (U_eff/2) * sum_{I,m,m',sigma} [ n^sigma_{mm'}(I) * (delta_{mm'} - n^sigma_{m'm}(I)) ]
///
/// where n^sigma_{mm'}(I) is the on-site occupation matrix for correlated orbital l on atom I.
///
/// The force on atom J is:
///
///   F_J = -dE_U/d tau_J
///       = -sum_{I,R} sum_{m,m'} V_U_{mm'}(I) * [
///             sum_{mu,nu} DMR_{mu,nu}(I,R) * d<phi_{mu,0}|chi_m(I)>/d tau_J * <chi_m'(I)|phi_{nu,R}>
///           ]
///
/// For stress, the derivative is with respect to strain tensor epsilon_{alpha,beta}:
///
///   sigma_{alpha,beta} = -(1/Omega) * dE_U/d epsilon_{alpha,beta}
///       = -(1/Omega) * sum_{I,R} sum_{m,m'} V_U_{mm'}(I) * [
///             sum_{mu,nu} DMR_{mu,nu}(I,R) * (
///                 d<phi_{mu,0}|chi_m(I)>/d epsilon_{alpha,beta} * <chi_m'(I)|phi_{nu,R}> * R_beta
///               + <phi_{mu,0}|chi_m(I)> * d<chi_m'(I)|phi_{nu,R}>/d epsilon_{alpha,beta} * R_beta
///             )
///           ]

#ifndef DFTU_NAO_FS_R_H
#define DFTU_NAO_FS_R_H

#include "source_base/matrix.h"

namespace hamilt
{

// Forward declarations to avoid circular dependency with dftu_lcao_op.h
template <typename TK, typename TR>
class OperatorLCAO;

template <typename T>
class DFTU;

/**
 * @brief Calculate DFT+U force and stress in real space (unified for gamma-only and multik)
 *
 * This is the unified entry for DFT+U force/stress calculation. It loops over all
 * on-site atoms with correlated orbitals, computes the two-center integrals <phi|chi>
 * via TwoCenterIntegrator, and accumulates force/stress contributions from all
 * atom pairs (I,J,R) using OpenMP parallelization.
 *
 * @note This implementation uses the real-space density matrix DMR (HContainer<double>)
 *       and two-center integrals <phi|chi> computed by TwoCenterIntegrator. It is
 *       independent of k-point sampling because DMR already contains the BZ integration.
 *
 * @param dftu_op     [in] pointer to the DFTU operator object (for accessing ucell, dftu, intor_)
 * @param cal_force   [in] whether to compute force
 * @param cal_stress  [in] whether to compute stress
 * @param force       [out] force matrix (nat, 3), accumulated
 * @param stress      [out] stress matrix (3, 3), accumulated
 *
 * @warning The density matrix must be set via Plus_U::set_dmr() before calling this.
 *          If get_dmr(0) returns nullptr, the function aborts with WARNING_QUIT.
 */
template <typename TK, typename TR>
void cal_fs_nao_r(DFTU<OperatorLCAO<TK, TR>>* dftu_op,
                        const bool cal_force,
                        const bool cal_stress,
                        ModuleBase::matrix& force,
                        ModuleBase::matrix& stress);

} // namespace hamilt

#endif // DFTU_NAO_FS_R_H

#ifndef DFTU_NAO_FS_R_H
#define DFTU_NAO_FS_R_H

/// @file dftu_nao_fs_r.h
/// @brief Unified entry for DFT+U force and stress in real space (r-space)
///
/// Loops over atom pairs (I,J,R) and dispatches to the per-pair kernels
/// cal_for_IJR_nao_r / cal_str_IJR_nao_r. Independent of k-point sampling
/// because the real-space density matrix (DMR) already contains the
/// Brillouin-zone integration.
///
/// Naming convention: _r suffix denotes the real-space implementation,
/// corresponding to _k for the k-space (legacy) one.
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

#include "source_base/matrix.h"

#include <vector>

class UnitCell;
class Plus_U_Base;
class TwoCenterIntegrator;
class AdjacentAtomInfo;

namespace hamilt
{
// Forward declarations to avoid circular dependency with dftu_lcao_op.h
template <typename TK, typename TR>
class OperatorLCAO;

template <typename T>
class DFTU_onsite;

template <typename T>
class HContainer;
} // namespace hamilt

namespace DFTU_LCAO
{

/**
 * @brief Post-process force: MPI Allreduce and spin-degeneracy scaling.
 *
 * @param force  [in,out] force matrix (nat, 3), reduced and scaled
 * @param nspin  [in] number of spin channels (1, 2, or 4)
 */
void reduce_force_impl(ModuleBase::matrix& force, int nspin);

/**
 * @brief Post-process stress: MPI Allreduce and renormalization to
 *        full 3x3 tensor from Voigt-like 6-component form.
 *
 * @param ucell      [in] unit cell (for lat0/omega)
 * @param stress_tmp [in] Voigt-like 6-component stress, reduced in-place
 * @param stress     [out] full 3x3 stress tensor
 */
void reduce_stress_impl(const UnitCell* ucell,
                        const std::vector<double>& stress_tmp,
                        ModuleBase::matrix& stress);

/**
 * @brief Non-template core of DFT+U force/stress in real space.
 *
 * All types are concrete and independent of the operator's k-point type (TK)
 * and real-space type (TR). The template wrapper cal_fs_nao_r only validates
 * the density matrix and forwards arguments here.
 *
 * @note Hard to unit-test: requires TwoCenterIntegrator (numerical orbital
 * radial integration) to build nlm via build_nlm. The inner function
 * acc_fs_pairs is directly testable with a mock nlm_tot; build_nlm is not.
 *
 * @param ucell       [in] unit cell
 * @param dftu        [in] DFT+U base object (occupation matrix, U values)
 * @param intor       [in] two-center integrator for <phi|chi> and gradients
 * @param nspin       [in] number of spin channels (1, 2, or 4)
 * @param adjs_all    [in] adjacent atom info for all atoms with plus-U
 * @param dmR         [in] density matrices in real space, size nspin
 * @param cal_force   [in] whether to compute force
 * @param cal_stress  [in] whether to compute stress
 * @param force       [out] force matrix (nat, 3), accumulated
 * @param stress      [out] stress matrix (3, 3), accumulated
 */
void cal_fs_nao_r_impl(const UnitCell* ucell,
                       Plus_U_Base* dftu,
                       const TwoCenterIntegrator* intor,
                       int nspin,
                       const std::vector<AdjacentAtomInfo>& adjs_all,
                       const std::vector<const hamilt::HContainer<double>*>& dmR,
                       bool cal_force,
                       bool cal_stress,
                       ModuleBase::matrix& force,
                       ModuleBase::matrix& stress);

/**
 * @brief Calculate DFT+U force and stress in real space from explicit
 *        environment arguments (non-template overload).
 *
 * This overload does not require a DFTU operator object; it is intended for
 * callers that only need force/stress and already hold the required data.
 *
 * @param ucell       [in] unit cell
 * @param dftu        [in] DFT+U base object (occupation matrix, U values)
 * @param intor       [in] two-center integrator for <phi|chi> and gradients
 * @param nspin       [in] number of spin channels (1, 2, or 4)
 * @param adjs_all    [in] adjacent atom info for all atoms with plus-U
 * @param dmR         [in] density matrices in real space, size nspin
 * @param cal_force   [in] whether to compute force
 * @param cal_stress  [in] whether to compute stress
 * @param force       [out] force matrix (nat, 3), accumulated
 * @param stress      [out] stress matrix (3, 3), accumulated
 */
void cal_fs_nao_r(const UnitCell* ucell,
                  Plus_U_Base* dftu,
                  const TwoCenterIntegrator* intor,
                  int nspin,
                  const std::vector<AdjacentAtomInfo>& adjs_all,
                  const std::vector<const hamilt::HContainer<double>*>& dmR,
                  bool cal_force,
                  bool cal_stress,
                  ModuleBase::matrix& force,
                  ModuleBase::matrix& stress);

} // namespace DFTU_LCAO

#endif // DFTU_NAO_FS_R_H

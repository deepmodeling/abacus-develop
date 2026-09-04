/**
 * @file deltaspin_pw_mi.h
 * @brief PW-basis DeltaSpin computation path as free functions, decoupled from SpinConstrain.
 *
 * @par Purpose
 * These functions implement the PW-basis branch of DeltaSpin that used to be
 * member functions of SpinConstrain<std::complex<double>>, defined out-of-line in
 * source/source_pw/module_pwdft/deltaspin_pw_impl.cpp. They are now free functions
 * in spinconstrain::pw so the PW implementation lives in module_deltaspin instead
 * of source_pw (removing the reverse source_pw -> module_deltaspin implementation
 * dependency), and so dependencies are passed explicitly per the AGENTS.md rule.
 *
 * The PW path is always instantiated on complex wavefunctions
 * (psi::Psi<std::complex<double>>), so these functions are not templated on TK.
 *
 * @par Workflow
 * - cal_mi_pw(): compute atomic magnetic moments Mi from becp = <alpha|psi>.
 * - calculate_delta_hcc(): add the DeltaSpin correction H += becp^dagger * lambda * becp
 *   to a subspace Hamiltonian.
 * - update_psi_charge_pw_cpu()/update_psi_charge_pw_gpu(): subspace diagonalization +
 *   optional full-space refinement to update psi and charge density after the lambda loop.
 */
#ifndef DELTASPIN_PW_MI_H
#define DELTASPIN_PW_MI_H

#include <complex>

#include "source_base/vector3.h"
#include "deltaspin_state.h"

// Forward declarations to keep header dependencies minimal (AGENTS.md rule 3).
namespace psi
{
template <typename T, typename Device>
class Psi;
}
namespace hamilt
{
template <typename T, typename Device>
class Hamilt;
}
namespace elecstate
{
class ElecState;
}
namespace ModulePW
{
class PW_Basis_K;
}

namespace spinconstrain
{
namespace pw
{

class SubspaceCache;

/**
 * @brief Calculate atomic magnetic moments using projector overlap (PW basis).
 *
 * @details For each k-point: tabulate atomic projectors, compute becp = <alpha|psi>
 * via OnsiteProjector::overlap_proj_psi, then decompose becp into magnetic moments.
 * Finally Mi is summed across all MPI k-pool ranks.
 *
 * @param state   Constraint state; state.Mi_ is filled in place.
 * @param psi     PW wavefunctions (psi::Psi<std::complex<double>>*, passed as void*
 *                to keep this header free of the Device template parameter).
 * @param pelec   Electronic state (provides wg weights and k-list spin signs).
 */
void cal_mi_pw(ScState& state,
               void* psi,
               elecstate::ElecState* pelec);

/**
 * @brief Compute DeltaSpin correction to a subspace Hamiltonian.
 *
 * @details Adds H += becp^dagger * ps, where ps = delta_lambda * becp. For npol=2
 * uses the full 2x2 Pauli matrix; for npol=1 uses the z-component with spin_sign.
 *
 * @param state        Constraint state (npol, Mi size, lambda snapshot via cache).
 * @param cache        PW subspace cache; supplies lambda_in_sub() for full_update.
 * @param pelec        Electronic state (for collinear spin_sign lookup).
 * @param h_tmp        Subspace Hamiltonian (nbands x nbands, modified in place).
 * @param becp_k       Projector coefficients for k-point ik.
 * @param delta_lambda Lambda change per atom (or full lambda if full_update).
 * @param nbands       Number of bands.
 * @param nkb          Total number of projectors.
 * @param nh_iat       Number of projectors per atom.
 * @param ik           K-point index (for collinear spin_sign lookup).
 * @param full_update  If true, compute delta = lambda_current - lambda_at_save.
 */
void calculate_delta_hcc(ScState& state,
                         const SubspaceCache& cache,
                         elecstate::ElecState* pelec,
                         std::complex<double>* h_tmp,
                         const std::complex<double>* becp_k,
                         const ModuleBase::Vector3<double>* delta_lambda,
                         const int nbands,
                         const int nkb,
                         const int* nh_iat,
                         const int ik,
                         const bool full_update);

/**
 * @brief CPU implementation of PW wavefunction and charge density update.
 *
 * @par Two-stage process
 * Stage 1 - Subspace diagonalization: apply DeltaSpin correction to the saved
 *   subspace H, then diagonalize to rotate wavefunctions (cheap, nbands x nbands).
 * Stage 2 - Full-space update: if pw_solve, run HSolverPW for iterative refinement;
 *   else update weights from new eigenvalues and call psiToRho().
 *
 * Frees the subspace cache after use (allocated in cal_mw_from_lambda).
 *
 * @param state        Constraint state.
 * @param cache        PW subspace cache (released here).
 * @param psi          PW wavefunctions (psi::Psi<std::complex<double>>*, as void*).
 * @param p_hamilt     Hamiltonian (hamilt::Hamilt<std::complex<double>>*, as void*).
 * @param pelec        Electronic state.
 * @param pw_wfc       PW basis for wavefunction storage.
 * @param delta_lambda Lambda change for incremental H correction.
 * @param pw_solve     If true, run full PW solver; if false, just update weights.
 * @param full_update  If true, apply full lambda (not delta) to H correction.
 */
void update_psi_charge_pw_cpu(ScState& state,
                              SubspaceCache& cache,
                              void* psi,
                              void* p_hamilt,
                              elecstate::ElecState* pelec,
                              ModulePW::PW_Basis_K* pw_wfc,
                              const ModuleBase::Vector3<double>* delta_lambda,
                              bool pw_solve,
                              bool full_update);

#if ((defined __CUDA) || (defined __ROCM))
/**
 * @brief GPU implementation of PW wavefunction and charge density update.
 * @details Same algorithm as update_psi_charge_pw_cpu(), but with GPU memory
 * management (device allocation, host-device synchronization).
 */
void update_psi_charge_pw_gpu(ScState& state,
                              SubspaceCache& cache,
                              void* psi,
                              void* p_hamilt,
                              elecstate::ElecState* pelec,
                              ModulePW::PW_Basis_K* pw_wfc,
                              const ModuleBase::Vector3<double>* delta_lambda,
                              bool pw_solve,
                              bool full_update);
#endif // __CUDA || __ROCM

} // namespace pw
} // namespace spinconstrain

#endif // DELTASPIN_PW_MI_H

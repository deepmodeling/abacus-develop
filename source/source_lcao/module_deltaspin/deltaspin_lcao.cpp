#include "deltaspin_lcao.h"
#include "spin_constrain.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/elecstate.h"

/**
 * @file deltaspin_lcao.cpp
 * @brief Wrapper/facade layer between ESolver and DeltaSpin module.
 *
 * @par Purpose
 * Provides a simplified interface to the ESolver layer, hiding the
 * SpinConstrain Singleton details. The ESolver calls these functions
 * rather than accessing SpinConstrain directly.
 *
 * @par Design rationale
 * - Template functions: Support both TK=double (nspin=2) and TK=complex<double> (nspin=4)
 * - Early returns: If sc_mag_switch is false, all functions return immediately
 *   without any overhead
 * - #ifdef __LCAO: The density matrix pointer is only available in LCAO builds
 *
 * @par Workflow
 * 1. ESolver calls init_deltaspin_lcao() at start of calculation
 * 2. Each SCF iteration:
 *    a. ESolver calls cal_mi_lcao_wrapper() to compute magnetic moments
 *    b. ESolver calls run_deltaspin_lambda_loop_lcao() to optimize lambda
 *    c. If skip_solve=true, ESolver skips the Hamiltonian solve (lambda loop already did it)
 */

namespace ModuleESolver
{

/**
 * @brief Initialize the SpinConstrain singleton with all input parameters.
 *
 * @details Called once at the start of a DeltaSpin calculation. Checks
 * sc_mag_switch first; if disabled, returns immediately without any action.
 *
 * @par Conditional compilation
 * The density matrix pointer (dm) is only available when __LCAO is defined.
 * For non-LCAO builds (PW-only), init_sc() is called without the dm parameter.
 */
template <typename TK>
void init_deltaspin_lcao(const UnitCell& ucell,
                          const Input_para& inp,
                          void* pv,
                          const K_Vectors& kv,
                          void* p_hamilt,
                          void* psi,
                          void* dm,
                          void* pelec)
{
    // Early exit if DeltaSpin is not enabled
    if (!inp.sc_mag_switch)
    {
        return;
    }

    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
#ifdef __LCAO
    // LCAO build: pass density matrix pointer
    sc.init_sc(inp.sc_thr, inp.nsc, inp.nsc_min, inp.alpha_trial,
               inp.sccut, inp.sc_drop_thr, ucell, inp.sc_direction_only,
               static_cast<Parallel_Orbitals*>(pv),
               inp.nspin, kv, p_hamilt, psi,
               static_cast<elecstate::DensityMatrix<TK, double>*>(dm),
               static_cast<elecstate::ElecState*>(pelec));
#else
    // Non-LCAO build: no density matrix
    sc.init_sc(inp.sc_thr, inp.nsc, inp.nsc_min, inp.alpha_trial,
               inp.sccut, inp.sc_drop_thr, ucell, inp.sc_direction_only,
               static_cast<Parallel_Orbitals*>(pv),
               inp.nspin, kv, p_hamilt, psi,
               static_cast<elecstate::ElecState*>(pelec));
#endif
}

/**
 * @brief Wrapper: calculate magnetic moments for current SCF iteration.
 *
 * @details If DeltaSpin is enabled, calls SpinConstrain::cal_mi_lcao().
 * The moments are stored in Mi_ and can be retrieved via get_target_mag().
 */
template <typename TK>
void cal_mi_lcao_wrapper(const int iter, const Input_para& inp)
{
    if (!inp.sc_mag_switch)
    {
        return;
    }

#ifdef __LCAO
    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
    sc.cal_mi_lcao(iter);
#endif
}

/**
 * @brief Wrapper: run the lambda optimization loop.
 *
 * @details Decision logic for when to run the lambda loop:
 *
 *   Case 1: NOT converged AND charge density is close enough (drho < sc_scf_thr)
 *   -> Run lambda loop, mark as converged, skip_solve = true
 *   Rationale: The charge density is stable enough to optimize lambda.
 *   The lambda loop does its own diagonalization, so skip the outer solve.
 *
 *   Case 2: Already converged
 *   -> Still run lambda loop (to refine for the current charge density)
 *   -> skip_solve = true
 *   Rationale: Even if converged, the charge density may have changed
 *   slightly, requiring lambda refinement.
 *
 *   Case 3: NOT converged AND charge density is NOT close enough (drho >= sc_scf_thr)
 *   -> Do nothing, skip_solve = false
 *   Rationale: The charge density is still changing significantly, so
 *   optimizing lambda would be premature. Wait for SCF to stabilize first.
 *
 * @param iter Current SCF iteration number
 * @param drho Charge density convergence criterion (max|drho|)
 * @param inp Input parameters
 * @return true if the ESolver should skip the Hamiltonian solve
 */
template <typename TK>
bool run_deltaspin_lambda_loop_lcao(const int iter,
                                     const double drho,
                                     const Input_para& inp)
{
    bool skip_solve = false;

    if (inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();

        if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
        {
            /// Charge density is stable enough: optimize lambda for the first time
            sc.run_lambda_loop(iter);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            /// Already converged: refine lambda for the current charge density
            sc.run_lambda_loop(iter);
            skip_solve = true;
        }
    }

    return skip_solve;
}

/// Template instantiations for both spin types
template void init_deltaspin_lcao<double>(const UnitCell& ucell,
                                           const Input_para& inp,
                                           void* pv,
                                           const K_Vectors& kv,
                                           void* p_hamilt,
                                           void* psi,
                                           void* dm,
                                           void* pelec);
template void init_deltaspin_lcao<std::complex<double>>(const UnitCell& ucell,
                                                          const Input_para& inp,
                                                          void* pv,
                                                          const K_Vectors& kv,
                                                          void* p_hamilt,
                                                          void* psi,
                                                          void* dm,
                                                          void* pelec);

template void cal_mi_lcao_wrapper<double>(const int iter, const Input_para& inp);
template void cal_mi_lcao_wrapper<std::complex<double>>(const int iter, const Input_para& inp);

template bool run_deltaspin_lambda_loop_lcao<double>(const int iter,
                                                      const double drho,
                                                      const Input_para& inp);
template bool run_deltaspin_lambda_loop_lcao<std::complex<double>>(const int iter,
                                                                      const double drho,
                                                                      const Input_para& inp);

} // namespace ModuleESolver

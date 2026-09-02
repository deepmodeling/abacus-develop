#ifndef LAMBDA_LOOP_HELPER_H
#define LAMBDA_LOOP_HELPER_H

#include <ostream>
#include <vector>

#include "source_base/vector3.h"
#include "spin_constrain.h"

/**
 * @file lambda_loop_helper.h
 * @brief Free-function helpers for the DeltaSpin lambda optimization loop.
 *
 * @par Background
 * Originally these routines were member functions of SpinConstrain<TK>.
 * They have been lifted to free functions in the spinconstrain namespace to
 * shrink the SpinConstrain god class. Each helper takes the SpinConstrain
 * instance as its first parameter and accesses internal state through the
 * public getters (get_Mi, get_sc_lambda, get_constrain, ...).
 *
 * @par Template parameter TK
 * - std::complex<double>: full implementation (nspin=2 and nspin=4)
 * - double: stub specialization (no-ops) provided elsewhere
 */

namespace spinconstrain
{

/**
 * @brief Print final spin and lambda values when the lambda loop terminates.
 *
 * @param sc          SpinConstrain instance (read Mi_ and lambda_)
 * @param ofs_running Log output stream
 */
template <typename TK>
void print_termination(const SpinConstrain<TK>& sc, std::ostream& ofs_running);

/**
 * @brief Check whether RMS error is below the convergence threshold or the
 *        maximum number of inner steps has been reached.
 *
 * @param sc             SpinConstrain instance
 * @param outer_step     Current SCF outer iteration
 * @param i_step         Current inner lambda step
 * @param rms_error      Current RMS error of Mi - M_target
 * @param duration       Wall time for this step (s)
 * @param total_duration Cumulative wall time for the inner loop (s)
 * @param ofs_running    Log output stream
 * @return true if the inner loop should terminate, false otherwise
 */
template <typename TK>
bool check_rms_stop(const SpinConstrain<TK>& sc,
                    int outer_step,
                    int i_step,
                    double rms_error,
                    double duration,
                    double total_duration,
                    std::ostream& ofs_running);

/**
 * @brief Print header at the start of the lambda optimization loop.
 *
 * @param sc          SpinConstrain instance
 * @param ofs_running Log output stream
 */
template <typename TK>
void print_header(const SpinConstrain<TK>& sc, std::ostream& ofs_running);

/**
 * @brief Cap the step size to prevent the optimizer from overshooting.
 *
 * @details If |alpha_trial * max(search)| exceeds restrict_current_, the
 * trial step is reduced so that the maximum lambda change per step is
 * bounded. alpha_trial is modified in place.
 *
 * @param sc           SpinConstrain instance
 * @param search       Current search direction (per atom, 3 components)
 * @param alpha_trial  Trial step size, modified in place if capped
 * @param ofs_running  Log output stream
 */
template <typename TK>
void check_restriction(const SpinConstrain<TK>& sc,
                       const std::vector<ModuleBase::Vector3<double>>& search,
                       double& alpha_trial,
                       std::ostream& ofs_running);

/**
 * @brief Compute the optimal step size via two-point linear interpolation.
 *
 * @par Algorithm
 *   alpha_opt = sum_k / sum_k2 * alpha_trial
 * where
 *   sum_k  = sum((target - spin) . (spin_plus - spin))   over constrained components
 *   sum_k2 = sum(|spin - spin_plus|^2)                   over constrained components
 *
 * @param sc           SpinConstrain instance
 * @param spin         Mi at current lambda
 * @param spin_plus    Mi at trial lambda (current + alpha_trial * search)
 * @param alpha_trial  Current trial step size
 * @return Optimal step size; falls back to alpha_trial when sum_k2 ~ 0
 */
template <typename TK>
double cal_alpha_opt(const SpinConstrain<TK>& sc,
                     std::vector<ModuleBase::Vector3<double>> spin,
                     std::vector<ModuleBase::Vector3<double>> spin_plus,
                     const double alpha_trial);

/**
 * @brief Check whether the magnetic susceptibility gradient dM/dlambda has
 *        decayed below the per-atom-type threshold.
 *
 * @par Algorithm
 * 1. Compute spin_change = new_spin - spin
 * 2. Compute nu_change     = delta_lambda - dnu_last_step
 * 3. Build full gradient matrix dM[ia][ic]/dlambda[ja][jc]
 * 4. Extract diagonal; pick max abs per atom type
 * 5. Return true if max(|diag|) < decay_grad[itype] for any type
 *
 * @param sc             SpinConstrain instance
 * @param new_spin       Mi at current lambda
 * @param spin           Mi at previous lambda
 * @param delta_lambda   Current lambda change
 * @param dnu_last_step  Previous cumulative step
 * @param print          Whether to print detailed gradient info
 * @param ofs_running    Log output stream
 * @return true if gradient decayed below threshold for any atom type
 */
template <typename TK>
bool check_gradient_decay(const SpinConstrain<TK>& sc,
                          std::vector<ModuleBase::Vector3<double>> new_spin,
                          std::vector<ModuleBase::Vector3<double>> spin,
                          std::vector<ModuleBase::Vector3<double>> delta_lambda,
                          std::vector<ModuleBase::Vector3<double>> dnu_last_step,
                          bool print,
                          std::ostream& ofs_running);

/**
 * @brief Print atomic magnetic moments Mi in a formatted table.
 *
 * @par Output format
 * - nspin=2: "Total Magnetism (uB)" with single z-component column
 * - nspin=4: three columns (Mx, My, Mz)
 *
 * @param sc          SpinConstrain instance
 * @param ofs_running Log output stream
 */
template <typename TK>
void print_Mi(const SpinConstrain<TK>& sc, std::ostream& ofs_running);

/**
 * @brief Print the magnetic force (-lambda) per atom in eV/uB.
 *
 * @par Physical meaning
 * Magnetic force = dL/dMi = -lambda. Large |lambda| means the system
 * strongly resists the target moment constraint.
 *
 * @param sc          SpinConstrain instance
 * @param ofs_running Log output stream
 */
template <typename TK>
void print_Mag_Force(const SpinConstrain<TK>& sc, std::ostream& ofs_running);

} // namespace spinconstrain

#endif // LAMBDA_LOOP_HELPER_H

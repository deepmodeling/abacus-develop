#include "basic_funcs.h"
#include "spin_constrain.h"

/**
 * @file lambda_loop_helper.cpp
 * @brief Helper/auxiliary methods for the lambda optimization loop.
 *
 * @par Functions overview
 * - print_termination(): Print final spin and lambda values when loop exits
 * - check_rms_stop(): Check convergence and print step info
 * - print_header(): Print header at loop start
 * - check_restriction(): Cap step size to prevent overshooting
 * - cal_alpha_opt(): Compute optimal step size via linear interpolation
 * - check_gradient_decay(): Check if dM/dlambda has decayed below threshold
 */

/**
 * @brief Print final spin and lambda values when lambda loop terminates.
 *
 * @par Output
 * - "after-optimization spin (uB)": Final magnetic moments Mi for each atom
 * - "after-optimization lambda (eV/uB)": Final Lagrange multipliers for each atom
 * - "Inner optimization for lambda ends.": Termination marker
 *
 * @par Interpretation
 * - Mi close to target_mag: constraint successfully satisfied
 * - Mi far from target_mag: constraint not converged (check RMS error in log)
 * - lambda ≈ 0: system naturally has the target moment
 * - lambda large: system resists the constraint (may indicate unrealistic target)
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::print_termination()
{
    print_2d("after-optimization spin (uB): (print in the inner loop): ", this->Mi_, this->nspin_);
    print_2d("after-optimization lambda (eV/uB): (print in the inner loop): ", this->lambda_, this->nspin_, ModuleBase::Ry_to_eV);
    std::cout << "Inner optimization for lambda ends." << std::endl;
    std::cout << "===============================================================================" << std::endl;
}

/**
 * @brief Check if RMS error is below convergence threshold or max steps reached.
 *
 * @par Output
 * Prints step info: "Step (Outer -- Inner) = X -- Y   RMS = Z   TIME(s) = T"
 *
 * @par Termination messages
 * - "Meet convergence criterion": RMS < current_sc_thr_ (successfully converged)
 * - "Reach maximum number of steps": i_step == nsc_ - 1 (did not converge)
 *
 * @par Return value
 * - true: loop should terminate (either converged or max steps)
 * - false: continue optimization
 *
 * @param outer_step Current SCF outer iteration
 * @param i_step Current inner lambda step
 * @param rms_error Current RMS error of Mi - M_target
 * @param duration Time for this step
 * @param total_duration Cumulative time for inner loop
 */
template <>
bool spinconstrain::SpinConstrain<std::complex<double>>::check_rms_stop(int outer_step,
                                                                                  int i_step,
                                                                                  double rms_error,
                                                                                  double duration,
                                                                                  double total_duration)
{
    std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << std::left << std::setw(5) << i_step + 1
              << "       RMS = " << rms_error << "     TIME(s) = " << std::setw(11) << duration << std::endl;
    if (rms_error < this->current_sc_thr_ || i_step == this->nsc_ - 1)
    {
        if (rms_error < this->current_sc_thr_)
        {
            std::cout << "Meet convergence criterion ( < " << this->current_sc_thr_ << " ), exit.";
            std::cout << "       Total TIME(s) = " << total_duration << std::endl;
        }
        else if (i_step == this->nsc_ - 1)
        {
            std::cout << "Reach maximum number of steps ( " << this->nsc_ << " ), exit.";
            std::cout << "              Total TIME(s) = " << total_duration << std::endl;
        }
        this->print_termination();
        return true;
    }
    return false;
}

/// @brief Print header at start of lambda optimization loop
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::print_header()
{
    std::cout << "===============================================================================" << std::endl;
    std::cout << "Inner optimization for lambda begins ..." << std::endl;
    std::cout << "Covergence criterion for the iteration: " << this->sc_thr_ << std::endl;
}

/**
 * @brief Cap step size to prevent overshooting in lambda optimization.
 *
 * @details If |alpha_trial * max(search)| > restrict_current_, reduce alpha_trial
 * so that the maximum lambda change per step is bounded by restrict_current_.
 *
 * This prevents the optimizer from taking steps that are too large, which
 * could lead to oscillation or divergence.
 *
 * @par Output (when restriction is applied)
 * - "alpha after restrict = X eV/uB^2": The capped step size
 * - "boundary after = X eV/uB": The actual maximum lambda change
 *
 * @param search Current search direction
 * @param alpha_trial Trial step size (modified in place if capped)
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::check_restriction(
    const std::vector<ModuleBase::Vector3<double>>& search,
    double& alpha_trial)
{
    double boundary = std::abs(alpha_trial * maxval_abs_2d(search));

    if (this->restrict_current_ > 0 && boundary > this->restrict_current_)
    {
        alpha_trial = copysign(1.0, alpha_trial) * this->restrict_current_ / maxval_abs_2d(search);
        boundary = std::abs(alpha_trial * maxval_abs_2d(search));
        std::cout << "alpha after restrict = " << alpha_trial * ModuleBase::Ry_to_eV << std::endl;
        std::cout << "boundary after = " << boundary * ModuleBase::Ry_to_eV << std::endl;
    }
}

/**
 * @brief Compute optimal step size via linear interpolation.
 *
 * @par Algorithm
 * Uses the two-point linear interpolation (secant method) to find the
 * step size that would drive Mi to M_target:
 *
 *   alpha_opt = sum_k / sum_k2 * alpha_trial
 *
 * where:
 *   sum_k  = sum((target - spin) . (spin_plus - spin))   over constrained components
 *   sum_k2 = sum(|spin - spin_plus|^2)                   over constrained components
 *
 * This is equivalent to finding the minimum of a quadratic approximation
 * to E(lambda) along the search direction.
 *
 * @par Edge case handling
 * - If |sum_k2| < 1e-30: spin and spin_plus are nearly identical, meaning
 *   the lambda change has no effect on Mi. Return alpha_trial as fallback.
 *   This can happen if the system is already saturated or if lambda is too small.
 *
 * @param spin Mi at current lambda
 * @param spin_plus Mi at trial lambda (current + alpha_trial * search)
 * @param alpha_trial Current trial step size
 * @return Optimal step size alpha_opt
 */
template <>
double spinconstrain::SpinConstrain<std::complex<double>>::cal_alpha_opt(
    std::vector<ModuleBase::Vector3<double>> spin,
    std::vector<ModuleBase::Vector3<double>> spin_plus,
    const double alpha_trial)
{
    int nat = this->get_nat();
    const bool print = false;
    const double zero = 0.0;

    // Mask to only constrained components
    std::vector<ModuleBase::Vector3<double>> spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> target_spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin_plus_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_2(nat, 0.0);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->target_mag_, target_spin_mask);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, spin, spin_mask);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, spin_plus, spin_plus_mask);

    // Compute dot products for linear interpolation
    for (int ia = 0; ia < nat; ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            // sum_k: (target - current) . (trial - current)
            temp_1[ia][ic]
                = (target_spin_mask[ia][ic] - spin_mask[ia][ic]) * (spin_plus_mask[ia][ic] - spin_mask[ia][ic]);
            // sum_k2: |current - trial|^2
            temp_2[ia][ic] = std::pow(spin_mask[ia][ic] - spin_plus_mask[ia][ic], 2);
        }
    }
    double sum_k = sum_2d(temp_1);
    double sum_k2 = sum_2d(temp_2);

    // Debug output (controlled by print flag)
    for(int ia=0; ia<std::min(2,(int)nat); ++ia) {
        if (print) {
        printf("[ALPHA-OPT] nat=%d sum_k=%.6e sum_k2=%.6e alpha_trial=%.6e\n", nat, sum_k, sum_k2, alpha_trial);
        printf("[ALPHA-OPT] spin[%d]=(%.4f,%.4f,%.4f) spin_plus[%d]=(%.4f,%.4f,%.4f)\n",
                ia, spin[ia].x, spin[ia].y, spin[ia].z,
                ia, spin_plus[ia].x, spin_plus[ia].y, spin_plus[ia].z);
        }
    }

    // Guard against division by zero
    if (std::abs(sum_k2) < 1e-30) {
        if (print) {
        printf("[ALPHA-OPT] WARNING: sum_k2 too small, returning alpha_trial\n");
        }
        fflush(stdout);
        return alpha_trial;
    }
    fflush(stdout);
    return sum_k * alpha_trial / sum_k2;
}

/**
 * @brief Check if the magnetic susceptibility gradient dM/dlambda has decayed below threshold.
 *
 * @par Algorithm
 * 1. Compute spin_change = new_spin - spin (change in magnetic moments)
 * 2. Compute nu_change = delta_lambda - dnu_last_step (change in lambda)
 * 3. Compute full gradient matrix: dM[ia][ic]/dlambda[ja][jc] = spin_change[ia][ic] / nu_change[ja][jc]
 * 4. Extract diagonal: dM[ia][ic]/dlambda[ia][ic] (self-susceptibility)
 * 5. Find max diagonal gradient per atom type
 * 6. If max_gradient[itype] < decay_grad[itype], return true (early termination)
 *
 * @par Physical meaning
 * The diagonal gradient dM/dlambda represents how sensitive the magnetic moment
 * is to changes in the Lagrange multiplier. When this gradient becomes very small,
 * further increases in lambda produce diminishing returns in Mi, indicating that
 * the optimization has reached its practical limit.
 *
 * @par Output (when triggered)
 * "Reach limitation of current step ( maximum gradient < X uB^2/eV in atom type Y ), exit."
 *
 * @par Debug output [GRAD-DECAY]
 * - WARNING: nu_change too small: indicates delta_lambda and dnu_last_step are
 *   nearly identical, meaning the optimizer is not making progress. This can happen
 *   if alpha_trial has become very small or if the search direction is nearly zero.
 *   Solution: check that alpha_trial is not vanishing; increase sc_thr if target
 *   is physically unreachable.
 *
 * @param new_spin Mi at current lambda
 * @param spin Mi at previous lambda
 * @param delta_lambda Current lambda change
 * @param dnu_last_step Previous cumulative step
 * @param print Whether to print detailed gradient info
 * @return true if gradient decayed below threshold (should terminate), false otherwise
 */
template <>
bool spinconstrain::SpinConstrain<std::complex<double>>::check_gradient_decay(
    std::vector<ModuleBase::Vector3<double>> new_spin,
    std::vector<ModuleBase::Vector3<double>> spin,
    std::vector<ModuleBase::Vector3<double>> delta_lambda,
    std::vector<ModuleBase::Vector3<double>> dnu_last_step,
    bool print)
{
    const double one = 1.0;
    const double zero = 0.0;
    int nat = this->get_nat();
    int ntype = this->get_ntype();

    // Change in magnetic moments and lambda
    std::vector<ModuleBase::Vector3<double>> spin_change(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> nu_change(nat, 1.0);

    // Full gradient matrix: dM[ia][ic]/dlambda[ja][jc]
    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
        nat,
        std::vector<std::vector<std::vector<double>>>(
            3,
            std::vector<std::vector<double>>(nat, std::vector<double>(3, 0.0))));
    // Diagonal gradient: dM[ia][ic]/dlambda[ia][ic] (self-susceptibility)
    std::vector<ModuleBase::Vector3<double>> spin_nu_gradient_diag(nat, 0.0);
    std::vector<std::pair<int, int>> max_gradient_index(ntype, std::make_pair(0, 0));
    std::vector<double> max_gradient(ntype, 0.0);

    subtract_2d(new_spin, spin, spin_change);
    subtract_2d(delta_lambda, dnu_last_step, nu_change);

    // Mask unconstrained components
    where_fill_scalar_2d(this->constrain_, 0, zero, spin_change);
    where_fill_scalar_2d(this->constrain_, 0, one, nu_change);

    // Calculate full gradient matrix
    for (int ia = 0; ia < nat; ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            for (int ja = 0; ja < nat; ja++)
            {
                for (int jc = 0; jc < 3; jc++)
                {
                    if (std::abs(nu_change[ja][jc]) < 1e-30) {
                        printf("[GRAD-DECAY] WARNING: nu_change[%d][%d] too small! delta_lambda=(%.6e,%.6e,%.6e) dnu_last=(%.6e,%.6e,%.6e)\n",
                               ja, jc, delta_lambda[ja].x, delta_lambda[ja].y, delta_lambda[ja].z,
                               dnu_last_step[ja].x, dnu_last_step[ja].y, dnu_last_step[ja].z);
                        fflush(stdout);
                        nu_change[ja][jc] = 1e-30;
                    }
                    spin_nu_gradient[ia][ic][ja][jc] = spin_change[ia][ic] / nu_change[ja][jc];
                }
            }
        }
    }

    // Extract diagonal gradient and find max per atom type
    for (const auto& sc_elem: this->get_atomCounts())
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        max_gradient[it] = 0.0;
        for (int ia = 0; ia < nat_it; ia++)
        {
            for (int ic = 0; ic < 3; ic++)
            {
                spin_nu_gradient_diag[ia][ic] = spin_nu_gradient[ia][ic][ia][ic];
                if (std::abs(spin_nu_gradient_diag[ia][ic]) > std::abs(max_gradient[it]))
                {
                    max_gradient[it] = spin_nu_gradient_diag[ia][ic];
                    max_gradient_index[it].first = ia;
                    max_gradient_index[it].second = ic;
                }
            }
        }
    }

    if (print)
    {
        print_2d("diagonal gradient: ", spin_nu_gradient_diag, this->nspin_);
        std::cout << "maximum gradient appears at: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            std::cout << "( " << max_gradient_index[it].first << ", " << max_gradient_index[it].second << " )"
                      << std::endl;
        }
        std::cout << "maximum gradient: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            std::cout << max_gradient[it]/ModuleBase::Ry_to_eV << std::endl;
        }
    }

    // Check if any atom type's gradient has decayed below threshold
    for (int it = 0; it < ntype; it++)
    {
        if (this->decay_grad_[it] > 0 && std::abs(max_gradient[it]) < this->decay_grad_[it])
        {
            std::cout << "Reach limitation of current step ( maximum gradient < " << this->decay_grad_[it]/ModuleBase::Ry_to_eV // uB^2/Ry to uB^2/eV
                      << " in atom type " << it << " ), exit." << std::endl;
            return true;
        }
    }
    return false;
}

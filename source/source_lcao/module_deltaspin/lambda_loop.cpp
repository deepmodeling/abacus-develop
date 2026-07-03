#include "spin_constrain.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "basic_funcs.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"

/**
 * @file lambda_loop.cpp
 * @brief Core lambda optimization algorithms for DeltaSpin.
 *
 * @par run_lambda_loop: Conjugate-gradient-like BFGS optimizer
 * Iteratively adjusts Lagrange multipliers (lambda) to drive atomic magnetic
 * moments (Mi) toward target values (M_target).
 *
 * @par Algorithm overview
 * The optimizer follows a modified Polak-Ribiere conjugate gradient scheme:
 *
 *   Step -1 (Initialization):
 *     - Compute initial Mi from current wavefunction
 *     - Save initial lambda (lambda with unconstrained components zeroed)
 *     - Set adaptive convergence threshold: current_sc_thr_ = max(rms_0 * sc_drop_thr_, sc_thr_)
 *
 *   Each inner step (i_step = 0, 1, ..., nsc-1):
 *     1. Update lambda: lambda = initial_lambda + delta_lambda
 *     2. [direction_only] Project out parallel component of lambda
 *     3. cal_mw_from_lambda(): apply lambda -> solve -> compute new Mi
 *     4. Check gradient decay: if dM/dlambda < decay_grad, exit early
 *     5. Compute residual: delta_spin = Mi - M_target
 *     6. Compute RMS error: rms = sqrt(mean(delta_spin^2))
 *     7. Check convergence: if rms < current_sc_thr_, update_psi_charge() and exit
 *        [PW basis] Re-check with cal_mi_pw(), recursively rerun if RMS too large
 *     8. [i_step >= 2] Compute Polak-Ribiere beta = rms^2 / rms_old^2
 *     9. Update search direction: search = delta_spin + beta * search_old
 *     10. Apply restriction: cap alpha_trial so that |alpha_trial * search| < restrict_current_
 *     11. Compute cumulative step: dnu = dnu + alpha_trial * search
 *     12. [direction_only] Project out parallel component of dnu
 *     13. Trial step: compute Mi at dnu, find optimal alpha via linear interpolation
 *     14. Update dnu with optimal alpha
 *     15. Adapt alpha_trial: if |alpha_opt| >> alpha_trial, increase; else decrease
 *
 * @par Key variables
 * - initial_lambda: lambda with unconstrained components set to 0
 * - delta_lambda: current lambda change from initial
 * - dnu: cumulative lambda change (search path integral)
 * - search: current search direction (steepest descent or conjugate)
 * - spin, spin_plus: Mi at current and trial lambda values
 * - alpha_trial: current step size (adaptively adjusted)
 * - alpha_opt: optimal step size from linear interpolation
 *
 * @par Convergence criteria
 * 1. RMS(Mi - M_target) < current_sc_thr_ (adaptive threshold)
 * 2. Maximum gradient dM/dlambda < decay_grad[itype] per atom type
 * 3. Maximum steps reached (nsc)
 *
 * @par Error output and solutions
 * - "RMS error is too large, rerun the loop": The subspace diagonalization
 *   was not accurate enough. The loop is rerun with rerun=false to use the
 *   full PW solver for better precision. If this persists, check:
 *   - PW_DIAG_NMAX and PW_DIAG_THR in DiagoIterAssist
 *   - higher_mag_prec flag for forced high-precision mode
 * - "Reach maximum number of steps": Lambda optimization did not converge
 *   within nsc steps. Check:
 *   - target_mag values are physically reasonable
 *   - alpha_trial is not too small (slow convergence)
 *   - decay_grad thresholds are not too aggressive
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop(int outer_step, bool rerun)
{
    int nat = this->get_nat();
    int ntype = this->get_ntype();

    // =============================================================
    // STATE VECTORS (all sized [nat][3])
    // =============================================================
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat,0.0); ///< Lambda with unconstrained components = 0
    std::vector<ModuleBase::Vector3<double>> delta_lambda(nat,0.0);   ///< Current lambda change from initial
    std::vector<ModuleBase::Vector3<double>> dnu(nat, 0.0), dnu_last_step(nat, 0.0); ///< Cumulative step, previous step
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);        ///< Temporary workspace
    std::vector<ModuleBase::Vector3<double>> spin(nat, 0.0), delta_spin(nat, 0.0);   ///< Current Mi, residual (Mi - M_target)
    std::vector<ModuleBase::Vector3<double>> search(nat, 0.0), search_old(nat, 0.0); ///< Search direction, previous direction
    std::vector<ModuleBase::Vector3<double>> new_spin(nat, 0.0), spin_plus(nat, 0.0); ///< Mi at current and trial lambda

    double alpha_opt, alpha_plus;  ///< Optimal step size, correction to trial
    double beta;                    ///< Polak-Ribiere conjugate gradient parameter
    double mean_error, mean_error_old, rms_error; ///< Mean squared error, RMS error
    double g;                       ///< Adaptation factor for alpha_trial

    double alpha_trial = this->alpha_trial_; ///< Current trial step size (Ry/uB^2)

    const double zero = 0.0;
    const double one = 1.0;

    // Timer initialization (MPI or CPU)
#ifdef __MPI
	auto iterstart = MPI_Wtime();
#else
	auto iterstart = std::chrono::system_clock::now();
#endif

    double inner_loop_duration = 0.0;

    this->print_header();

    // =============================================================
    // MAIN OPTIMIZATION LOOP
    // i_step = -1: initialization (compute initial Mi, save initial lambda)
    // i_step = 0, 1, ..., nsc-1: optimization steps
    // =============================================================
    for (int i_step = -1; i_step < this->nsc_; i_step++)
    {
        double duration = 0.0;
        if (i_step == -1)
        {
            // =============================================================
            // STEP -1: INITIALIZATION
            // Compute initial magnetic moments and save starting state
            // =============================================================
            this->cal_mw_from_lambda(i_step);
            spin = this->Mi_;

            // Save initial lambda: for unconstrained components (constrain==0), set to 0
            where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);

            print_2d("initial lambda (eV/uB): ", initial_lambda, this->nspin_, ModuleBase::Ry_to_eV);
            print_2d("initial spin (uB): ", spin, this->nspin_);
            print_2d("target spin (uB): ", this->target_mag_, this->nspin_);
            i_step++;
        }
        else
        {
            // =============================================================
            // OPTIMIZATION STEP
            // Update lambda, compute new Mi, check convergence
            // =============================================================

            // Mask unconstrained components of delta_lambda to 0
            where_fill_scalar_2d(this->constrain_, 0, zero, delta_lambda);

            // lambda = initial_lambda + delta_lambda
            add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);

            // [direction_only mode] Project out parallel component of lambda
            // This keeps |lambda| -> 0, only constraining spin direction
            if(this->direction_only_)
            for (int ia = 0; ia < nat; ia++)
            {
                const auto& target = this->target_mag_[ia];
                const double norm = std::sqrt(target.x*target.x + target.y*target.y + target.z*target.z);

                if (norm > 1e-8) {
                    const ModuleBase::Vector3<double> dir = target / norm;
                    double parallel = this->lambda_[ia].x*dir.x +
                                    this->lambda_[ia].y*dir.y +
                                    this->lambda_[ia].z*dir.z;
                    this->lambda_[ia].x -= parallel * dir.x;
                    this->lambda_[ia].y -= parallel * dir.y;
                    this->lambda_[ia].z -= parallel * dir.z;
                }
            }

            // Apply lambda and compute new magnetic moments
            this->cal_mw_from_lambda(i_step, delta_lambda.data());
            new_spin = this->Mi_;

            // Check if gradient dM/dlambda has decayed below threshold
            bool GradLessThanBound = this->check_gradient_decay(new_spin, spin, delta_lambda, dnu_last_step);
            if (i_step >= this->nsc_min_ && GradLessThanBound)
            {
                // Gradient has decayed: further optimization yields diminishing returns
                // Apply the last successful step and exit
                add_scalar_multiply_2d(initial_lambda, dnu_last_step, one, this->lambda_);
                this->update_psi_charge(dnu_last_step.data(), true, true);
#ifdef __MPI
		        duration = (double)(MPI_Wtime() - iterstart);
#else
			    duration =
                    (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now()
                    - iterstart)).count() / static_cast<double>(1e6);
#endif
                inner_loop_duration += duration;
                std::cout << "Total TIME(s) = " << inner_loop_duration << std::endl;
                this->print_termination();
                break;
            }
            spin = new_spin;
        }

        // =============================================================
        // COMPUTE RESIDUAL AND RMS ERROR
        // =============================================================
        // delta_spin = spin - target_mag (residual error)
        subtract_2d(spin, this->target_mag_, delta_spin);
        // Mask unconstrained components to 0 (they don't contribute to error)
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);

        // Search direction starts as the residual (steepest descent)
        search = delta_spin;

        // [direction_only mode] Modify residual to exclude parallel component
        // and adjust target direction without mutating target_mag_
        std::vector<ModuleBase::Vector3<double>> target_mag_adj = this->target_mag_;
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++)
        {
            const auto& target = this->target_mag_[ia];
            const double norm = std::sqrt(target.x*target.x + target.y*target.y + target.z*target.z);

            if (norm > 1e-8) {
                const ModuleBase::Vector3<double> dir = target / norm;
                const double parallel = delta_spin[ia].x*dir.x + delta_spin[ia].y*dir.y + delta_spin[ia].z*dir.z;
                // Store perpendicular component squared in temp_1 (for RMS)
                temp_1[ia][0] = std::pow(delta_spin[ia].x,2) + std::pow(delta_spin[ia].y,2) +
                                std::pow(delta_spin[ia].z,2) - std::pow(parallel,2);
                temp_1[ia][1] = 0;
                temp_1[ia][2] = 0;
                // Adjust target to include parallel component (work on copy, don't mutate target_mag_)
                target_mag_adj[ia] += parallel * dir;
            }
            else {
                temp_1[ia][0] = std::pow(delta_spin[ia].x,2) +
                              std::pow(delta_spin[ia].y,2) +
                              std::pow(delta_spin[ia].z,2);
                temp_1[ia][1] = 0;
                temp_1[ia][2] = 0;
            }
        }
        else
        for (int ia = 0; ia < nat; ia++)
        {
            for (int ic = 0; ic < 3; ic++)
            {
                temp_1[ia][ic] = std::pow(delta_spin[ia][ic],2);
            }
        }
        mean_error = sum_2d(temp_1) / nat;
        rms_error = std::sqrt(mean_error);

        // Set adaptive convergence threshold on first step
        if(i_step == 0)
        {
            this->current_sc_thr_ = std::max(rms_error * this->sc_drop_thr_, this->sc_thr_);
        }

        // =============================================================
        // CHECK CONVERGENCE
        // =============================================================
#ifdef __MPI
			duration = (double)(MPI_Wtime() - iterstart);
#else
			duration =
               (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now()
                - iterstart)).count() / static_cast<double>(1e6);
#endif
        inner_loop_duration += duration;
        if (this->check_rms_stop(outer_step, i_step, rms_error, duration, inner_loop_duration))
        {
            // Converged or max steps reached: final update
            this->update_psi_charge(dnu_last_step.data(), rerun, true);

            // [PW basis] Extra verification: re-compute Mi from scratch
            if(PARAM.inp.basis_type == "pw")
            {
                this->cal_mi_pw();
                subtract_2d(this->Mi_, this->target_mag_, delta_spin);
                where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
                search = delta_spin;
                for (int ia = 0; ia < nat; ia++)
                {
                    for (int ic = 0; ic < 3; ic++)
                    {
                        temp_1[ia][ic] = std::pow(delta_spin[ia][ic],2);
                    }
                }
                mean_error = sum_2d(temp_1) / nat;
                rms_error = std::sqrt(mean_error);
                std::cout<<"Current RMS: "<<rms_error<<std::endl;

                // If RMS is still large after full update, recursively rerun
                // with higher precision (full PW solver instead of subspace only)
                if(rms_error > this->current_sc_thr_ * 10 && rerun == true && this->higher_mag_prec == true)
                {
                    std::cout<<"Error: RMS error is too large, rerun the loop"<<std::endl;
                    this->run_lambda_loop(outer_step, false);
                }
            }
            break;
        }

        // Reset timer for next iteration
#ifdef __MPI
		iterstart = MPI_Wtime();
#else
		iterstart = std::chrono::system_clock::now();
#endif

        // =============================================================
        // POLAK-RIBIERE CONJUGATE GRADIENT UPDATE
        // =============================================================
        // For i_step >= 2, compute conjugate direction
        if (i_step >= 2)
        {
            // Polak-Ribiere beta = ||gradient_new||^2 / ||gradient_old||^2
            beta = mean_error / mean_error_old;
            // search = delta_spin + beta * search_old (conjugate direction)
            add_scalar_multiply_2d(search, search_old, beta, search);
        }

        // Cap step size to prevent overshooting
        this->check_restriction(search, alpha_trial);

        // =============================================================
        // CUMULATIVE STEP UPDATE
        // =============================================================
        dnu_last_step = dnu;
        // dnu = dnu + alpha_trial * search
        add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);

        // [direction_only mode] Project out parallel component from dnu
        // Use target_mag_adj (copy with parallel components added) instead of mutating target_mag_
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++) {
            const auto& target = target_mag_adj[ia];
            const double norm = std::sqrt(target.x*target.x + target.y*target.y + target.z*target.z);

            if (norm > 1e-8) {
                const ModuleBase::Vector3<double> dir = target / norm;
                double parallel = dnu[ia].x*dir.x + dnu[ia].y*dir.y + dnu[ia].z*dir.z;
                dnu[ia].x -= parallel * dir.x;
                dnu[ia].y -= parallel * dir.y;
                dnu[ia].z -= parallel * dir.z;
            }
        }
        delta_lambda = dnu;

        // Mask unconstrained components
        where_fill_scalar_else_2d(this->constrain_, 0, zero, delta_lambda, delta_lambda);
        // Update lambda
        add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);

        // =============================================================
        // TRIAL STEP: compute Mi at trial position
        // =============================================================
        this->cal_mw_from_lambda(i_step, delta_lambda.data());
        spin_plus = this->Mi_;

        // Find optimal step size via linear interpolation
        alpha_opt = this->cal_alpha_opt(spin, spin_plus, alpha_trial);
        this->check_restriction(search, alpha_opt);

        // Correct dnu: dnu += (alpha_opt - alpha_trial) * search
        alpha_plus = alpha_opt - alpha_trial;
        scalar_multiply_2d(search, alpha_plus, temp_1);
        add_scalar_multiply_2d(dnu, temp_1, one, dnu);

        // [direction_only] Project out parallel component from corrected dnu
        // Use target_mag_adj (copy) instead of mutating target_mag_
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++) {
            const auto& target = target_mag_adj[ia];
            const double norm = std::sqrt(target.x*target.x + target.y*target.y + target.z*target.z);

            if (norm > 1e-8) {
                const ModuleBase::Vector3<double> dir = target / norm;
                double parallel = dnu[ia].x*dir.x + dnu[ia].y*dir.y + dnu[ia].z*dir.z;
                dnu[ia].x -= parallel * dir.x;
                dnu[ia].y -= parallel * dir.y;
                dnu[ia].z -= parallel * dir.z;
            }
        }
        delta_lambda = dnu;

        // =============================================================
        // ADAPT STEP SIZE FOR NEXT ITERATION
        // =============================================================
        search_old = search;
        mean_error_old = mean_error;

        // Adapt alpha_trial based on ratio of optimal to trial step
        // g = 1.5 * |alpha_opt| / alpha_trial
        // - g > 2.0: alpha_opt was much larger than alpha_trial -> increase alpha_trial
        // - g < 0.5: alpha_opt was much smaller -> decrease alpha_trial
        // - 0.5 <= g <= 2.0: step size is reasonable -> modest adjustment
        g = 1.5 * std::abs(alpha_opt) / alpha_trial;
        if (g > 2.0)
        {
            g = 2;
        }
        else if (g < 0.5)
        {
            g = 0.5;
        }
        alpha_trial = alpha_trial * pow(g, 0.7);
    }

    return;
}

/**
 * @file lambda_loop.cpp (continued)
 * @brief Linear lambda scan mode for energy landscape mapping.
 *
 * @par Purpose
 * Instead of optimizing lambda to match target moments, this function
 * sweeps lambda values from sc_scan_lambda_start to sc_scan_lambda_end
 * in equal steps, computing Mi at each point. Useful for:
 * - Debugging: understanding the Mi vs lambda relationship
 * - Plotting: creating E(lambda) curves for analysis
 * - Validation: checking that Mi responds monotonically to lambda
 *
 * @par Output
 * Results written to lambda_scan_results.dat with columns:
 *   step, lambda_eV_uB, Mi_x_0, Mi_y_0, Mi_z_0, Mi_x_1, ...
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_linear_scan(int outer_step)
{
    int nat = this->get_nat();
    int ntype = this->get_ntype();

    double lambda_start = PARAM.inp.sc_scan_lambda_start;
    double lambda_end = PARAM.inp.sc_scan_lambda_end;
    int nsteps = PARAM.inp.sc_scan_steps;

    if (nsteps <= 0) {
        std::cout << "[DS-DIAG] linear_scan: sc_scan_steps <= 0, skipping" << std::endl;
        return;
    }

    // Convert eV to Ry for internal calculations
    double lambda_start_ry = lambda_start / ModuleBase::Ry_to_eV;
    double lambda_end_ry = lambda_end / ModuleBase::Ry_to_eV;
    double lambda_step = (lambda_end_ry - lambda_start_ry) / (nsteps - 1);

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === LINEAR LAMBDA SCAN START ===" << std::endl;
    std::cout << "[DS-DIAG] Scan range: " << lambda_start << " -> " << lambda_end << " eV/uB" << std::endl;
    std::cout << "[DS-DIAG] Number of steps: " << nsteps << std::endl;
    std::cout << "[DS-DIAG] Lambda step size: " << lambda_step * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;
    std::cout << "[DS-DIAG] nat = " << nat << ", ntype = " << ntype << std::endl;
    std::cout << "[DS-DIAG] nspin_ = " << this->nspin_ << ", npol_ = " << this->npol_ << std::endl;
    std::cout << "[DS-DIAG] p_operator = " << (this->p_operator ? "valid" : "NULL") << std::endl;
    std::cout << "[DS-DIAG] constrain_ size = " << this->constrain_.size() << std::endl;

    // Check if any constraints are defined; if not, set all atoms as constrained
    bool has_constraints = false;
    for (int ia = 0; ia < nat; ia++) {
        if (this->constrain_[ia].x != 0 || this->constrain_[ia].y != 0 || this->constrain_[ia].z != 0) {
            has_constraints = true;
            break;
        }
    }

    if (!has_constraints) {
        std::cout << "[DS-DIAG] No constraints found in STRU, setting all atoms as constrained" << std::endl;
        for (int ia = 0; ia < nat; ia++) {
            if (this->nspin_ == 4) {
                this->constrain_[ia] = ModuleBase::Vector3<int>(1, 1, 1);
            } else {
                this->constrain_[ia] = ModuleBase::Vector3<int>(0, 0, 1);
            }
        }
        this->reset_dspin_operator();
    }

    for (int ia = 0; ia < nat; ia++) {
        std::cout << "[DS-DIAG]   Atom " << ia << " constrain = ("
                  << this->constrain_[ia].x << ", " << this->constrain_[ia].y << ", " << this->constrain_[ia].z << ")"
                  << " target_mag = (" << this->target_mag_[ia].x << ", " << this->target_mag_[ia].y << ", " << this->target_mag_[ia].z << ")" << std::endl;
    }
    std::cout << std::string(80, '=') << "\n" << std::endl;

    // Save initial lambda to restore after scan
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    where_fill_scalar_else_2d(this->constrain_, 0, 0.0, this->lambda_, initial_lambda);

    // Open output file
    std::ofstream ofs_scan;
    if (outer_step == 0) {
        ofs_scan.open("lambda_scan_results.dat");
        ofs_scan << "# Linear Lambda Scan Results" << std::endl;
        ofs_scan << "# lambda_start = " << lambda_start << " eV/uB" << std::endl;
        ofs_scan << "# lambda_end = " << lambda_end << " eV/uB" << std::endl;
        ofs_scan << "# nsteps = " << nsteps << std::endl;
        ofs_scan << "#" << std::endl;
        ofs_scan << "# SCF iteration: " << outer_step << std::endl;
    } else {
        ofs_scan.open("lambda_scan_results.dat", std::ios::app);
        ofs_scan << "#" << std::endl;
        ofs_scan << "# SCF iteration: " << outer_step << std::endl;
    }

    // Write header
    ofs_scan << "# step  lambda_eV_uB";
    for (int ia = 0; ia < nat; ia++) {
        ofs_scan << "  Mi_x_" << ia << "  Mi_y_" << ia << "  Mi_z_" << ia;
    }
    ofs_scan << std::endl;

    double original_sc_thr = this->sc_thr_;

    // Save step 0 Mi for consistency check later
    std::vector<ModuleBase::Vector3<double>> mi_step0;

    // =============================================================
    // SCAN LOOP: sweep lambda from start to end
    // =============================================================
    for (int istep = 0; istep < nsteps; istep++) {
        double lambda_val_ry = lambda_start_ry + istep * lambda_step;
        double lambda_val_ev = lambda_val_ry * ModuleBase::Ry_to_eV;

        // Set lambda for all constrained atoms/components
        for (int ia = 0; ia < nat; ia++) {
            for (int ic = 0; ic < 3; ic++) {
                if (this->constrain_[ia][ic] != 0) {
                    this->lambda_[ia][ic] = lambda_val_ry;
                } else {
                    this->lambda_[ia][ic] = 0.0;
                }
            }
        }

        std::cout << "[DS-DIAG] === Scan step " << istep << "/" << nsteps
                  << " lambda = " << lambda_val_ev << " eV/uB ===" << std::endl;

        // Compute magnetic moments at current lambda
        this->cal_mw_from_lambda(istep);

        // Save step 0 Mi for consistency verification
        if (istep == 0) {
            mi_step0 = this->Mi_;
        }

        // Write results
        ofs_scan << std::scientific << std::setprecision(6);
        ofs_scan << istep << "  " << lambda_val_ev;
        for (int ia = 0; ia < nat; ia++) {
            ofs_scan << "  " << this->Mi_[ia].x
                     << "  " << this->Mi_[ia].y
                     << "  " << this->Mi_[ia].z;
        }
        ofs_scan << std::endl;

        std::cout << "[DS-DIAG]   lambda = " << lambda_val_ev << " eV/uB" << std::endl;
        for (int ia = 0; ia < nat; ia++) {
            std::cout << "[DS-DIAG]   Atom " << ia << " Mi = ("
                      << this->Mi_[ia].x << ", "
                      << this->Mi_[ia].y << ", "
                      << this->Mi_[ia].z << ") uB" << std::endl;
        }
        std::cout << std::endl;
    }

    // =============================================================
    // CONSISTENCY CHECK: restore initial lambda and recompute Mi
    // to verify that the lambda->Mi mapping is numerically stable
    // after multiple lambda updates in the scan loop
    // =============================================================
    std::cout << "[DS-DIAG] === Consistency check: restoring initial lambda ===" << std::endl;
    this->lambda_ = initial_lambda;
    this->cal_mw_from_lambda(nsteps);

    // Write consistency check result
    ofs_scan << std::scientific << std::setprecision(6);
    ofs_scan << "init_recheck  " << lambda_start;
    for (int ia = 0; ia < nat; ia++) {
        ofs_scan << "  " << this->Mi_[ia].x
                 << "  " << this->Mi_[ia].y
                 << "  " << this->Mi_[ia].z;
    }
    ofs_scan << std::endl;

    std::cout << "[DS-DIAG]   lambda = " << lambda_start << " eV/uB (restored)" << std::endl;
    for (int ia = 0; ia < nat; ia++) {
        std::cout << "[DS-DIAG]   Atom " << ia << " Mi = ("
                  << this->Mi_[ia].x << ", "
                  << this->Mi_[ia].y << ", "
                  << this->Mi_[ia].z << ") uB" << std::endl;
    }

    // Compare restored Mi with step 0 Mi to check consistency
    ofs_scan << "# [consistency] step 0 vs init_recheck Mi difference:" << std::endl;
    double max_mi_diff = 0.0;
    for (int ia = 0; ia < nat; ia++) {
        double dx = std::abs(this->Mi_[ia].x - mi_step0[ia].x);
        double dy = std::abs(this->Mi_[ia].y - mi_step0[ia].y);
        double dz = std::abs(this->Mi_[ia].z - mi_step0[ia].z);
        double diff = std::max({dx, dy, dz});
        if (diff > max_mi_diff) max_mi_diff = diff;
        ofs_scan << "#   Atom " << ia << " dM = (" << dx << ", " << dy << ", " << dz << ") uB" << std::endl;
    }
    std::cout << "[DS-DIAG] Max Mi difference between step 0 and init_recheck: " << max_mi_diff << " uB" << std::endl;
    if (max_mi_diff > 1e-8) {
        std::cout << "[DS-DIAG] WARNING: Mi mapping may be inconsistent after multiple lambda updates!" << std::endl;
    } else {
        std::cout << "[DS-DIAG] OK: Mi mapping is consistent." << std::endl;
    }
    ofs_scan << "#   Max Mi difference: " << max_mi_diff << " uB" << std::endl;

    ofs_scan.close();

    // Restore original lambda values (already restored above, but explicit for clarity)
    this->lambda_ = initial_lambda;

    std::cout << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === LINEAR LAMBDA SCAN COMPLETE ===" << std::endl;
    std::cout << "[DS-DIAG] Results written to: lambda_scan_results.dat" << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;

    return;
}

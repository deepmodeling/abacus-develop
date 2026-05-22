#include "spin_constrain.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "basic_funcs.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"

#ifdef __LCAO
#include "source_hsolver/hsolver_lcao.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#endif

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
/**
 * @brief Local diagnostic: compare methods near a reference lambda value.
 *
 * @details Builds subspace at lambda_ref, then scans lambda_ref ± 0.5 eV/uB.
 * Tests hypothesis: small perturbations near converged lambda are well-approximated.
 *
 * Methods compared:
 *   A: First-order eigenvalue response
 *   B: Full subspace diagonalization
 *   C: Full HSolverLCAO diagonalization (ground truth)
 *
 * @param outer_step Current SCF outer iteration number
 * @param lambda_ref_ry Reference lambda in Ry/uB (where subspace is built)
 * @param label Label for output file
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_local_diagnostic(
    int outer_step, double lambda_ref_ry, const std::string& label)
{
#ifdef __LCAO
    if (PARAM.inp.basis_type != "lcao") return;
    if (this->nspin_ != 2) return;
#else
    return;
#endif

    int nat = this->get_nat();
    double scan_half_range_ev = 0.5;
    int nsteps = 11;

    std::cout << "\n" << std::string(80, '-') << std::endl;
    std::cout << "[DS-LOCAL] === Local diagnostic: " << label << " ===" << std::endl;
    std::cout << "[DS-LOCAL] Subspace built at lambda_ref = " << lambda_ref_ry * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;
    std::cout << "[DS-LOCAL] Scan range: lambda_ref ± " << scan_half_range_ev << " eV/uB, " << nsteps << " steps" << std::endl;

    // Set lambda to reference value and build subspace
    this->free_lcao_subspace_cache();
    for (int ia = 0; ia < nat; ia++) {
        for (int ic = 0; ic < 3; ic++) {
            if (this->constrain_[ia][ic] != 0)
                this->lambda_[ia][ic] = lambda_ref_ry;
            else
                this->lambda_[ia][ic] = 0.0;
        }
    }

    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
    hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);

    dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
        ->update_lambda();

    hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
    hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
    elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                 this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                 this->pelec->skip_weights);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    // Cache subspace data at lambda_ref
    const int nk = psi_t->get_nk();
    const int nbands = PARAM.inp.nbands;
    const int nlocal = this->ParaV->nrow;
    const int nn = nbands * nbands;
    this->lcao_nbands_ = nbands;
    this->lcao_nk_ = nk;
    this->lcao_nlocal_ = nlocal;

    delete[] this->lcao_sub_h_save;
    delete[] this->lcao_sub_s_save;
    this->lcao_sub_h_save = new std::complex<double>[nk * nn];
    this->lcao_sub_s_save = new std::complex<double>[nk * nn];
    this->lcao_PI_sub_save_.resize(nk);
    this->lcao_ekb_save_.resize(nk * nbands);

    for (int ik = 0; ik < nk; ik++) {
        psi_t->fix_k(ik);
        this->calculate_lcao_sub_hs(
            this->p_hamilt, psi_t[0], this->ParaV,
            this->lcao_sub_h_save + ik * nn,
            this->lcao_sub_s_save + ik * nn,
            ik, nbands, nlocal);

        auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
            this->p_operator);
        dspin_op->cal_PI_sub(
            this->kv_.kvec_d[ik],
            psi_t->get_pointer(),
            nbands,
            this->lcao_PI_sub_save_[ik]);

        for (int ib = 0; ib < nbands; ib++)
            this->lcao_ekb_save_[ik * nbands + ib] = this->pelec->ekb(ik, ib);
    }
    this->lcao_lambda_in_sub_ = this->lambda_;
    this->lcao_subspace_initialized_ = true;

    // Save original psi wavefunctions at lambda_ref for restoration in each method
    const int nrow = this->ParaV->nrow;
    const int nloc_wfc = this->ParaV->nloc_wfc;
    std::vector<std::vector<std::complex<double>>> psi_save(nk);
    for (int ik = 0; ik < nk; ik++) {
        psi_t->fix_k(ik);
        psi_save[ik].assign(psi_t->get_pointer(), psi_t->get_pointer() + nloc_wfc);
    }

    // Helper lambda: restore psi, optionally rotate with V, then compute Mi via full-space DMR
    auto compute_mi_full_dm = [this, psi_t, nk, nbands, nrow, nloc_wfc, &psi_save](
        const std::vector<std::vector<std::complex<double>>>& vcc_list) -> std::vector<ModuleBase::Vector3<double>>
    {
        // Restore psi to lambda_ref state
        for (int ik = 0; ik < nk; ik++) {
            psi_t->fix_k(ik);
            std::complex<double>* psi_ptr = psi_t->get_pointer();
            std::memcpy(psi_ptr, psi_save[ik].data(), sizeof(std::complex<double>) * nloc_wfc);
        }

        // Rotate psi = C_old * V (no-op if V is identity)
        bool is_identity = true;
        for (int ik = 0; ik < nk && is_identity; ik++) {
            for (int ib = 0; ib < nbands; ib++) {
                if (std::abs(vcc_list[ik][ib + ib * nbands] - std::complex<double>{1.0, 0.0}) > 1e-14) {
                    is_identity = false;
                    break;
                }
                for (int jb = 0; jb < nbands; jb++) {
                    if (ib != jb && std::abs(vcc_list[ik][ib + jb * nbands]) > 1e-14) {
                        is_identity = false;
                        break;
                    }
                }
            }
        }

        if (!is_identity) {
            this->rotate_psi_subspace_lcao(*psi_t, this->ParaV, vcc_list, nbands, nrow, nk);
        }

        // Build full-space DMK from psi and weights, then DMR
        elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
        this->dm_->cal_DMR();

        // Compute magnetic moments via real-space projection (same as cal_mi_lcao)
        this->cal_mi_lcao(0);
        return this->Mi_;
    };

    // Open output file
    std::string fname = "local_diagnostic_" + label + ".dat";
    std::ofstream ofs_diag(fname);
    ofs_diag << "# Local Diagnostic: Methods near lambda_ref = " << lambda_ref_ry * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;
    ofs_diag << "# Methods: (A) First-order response, (B) Subspace diag, (C) Full HSolverLCAO" << std::endl;
    ofs_diag << "# Scan: Logarithmic deltas: 0, +/- 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0 eV/uB" << std::endl;
    ofs_diag << "#" << std::endl;
    ofs_diag << "# step  delta_lambda_eV  ";
    ofs_diag << "ekb_A_ib0  ekb_B_ib0  ekb_C_ib0  |  ";
    ofs_diag << "max_dek_AB  max_dek_AC  max_dek_BC  |  ";
    ofs_diag << "Mi_A_0  Mi_B_0  Mi_C_0  |  ";
    ofs_diag << "Mi_A_1  Mi_B_1  Mi_C_1  |  ";
    ofs_diag << "max_dMi_AB  max_dMi_AC  max_dMi_BC" << std::endl;

    std::vector<double> deltas_ev = {
        0.0,
        -0.00001, 0.00001,
        -0.0001, 0.0001,
        -0.001, 0.001,
        -0.01, 0.01,
        -0.1, 0.1,
        -1.0, 1.0
    };

    int step_count = 0;
    for (double delta_ev : deltas_ev) {
        double delta_lambda_ry = delta_ev / ModuleBase::Ry_to_eV;

        // Set lambda = lambda_ref + delta
        for (int ia = 0; ia < nat; ia++) {
            for (int ic = 0; ic < 3; ic++) {
                if (this->constrain_[ia][ic] != 0)
                    this->lambda_[ia][ic] = lambda_ref_ry + delta_lambda_ry;
                else
                    this->lambda_[ia][ic] = 0.0;
            }
        }

        // Method A: First-order eigenvalue response
        std::vector<std::vector<double>> ekb_A(nk);
        for (int ik = 0; ik < nk; ik++) {
            ekb_A[ik].resize(nbands);
            int spin_sign = this->get_spin_sign(ik);
            for (int ib = 0; ib < nbands; ib++) {
                double delta_epsilon = 0.0;
                for (int iat = 0; iat < nat; iat++) {
                    if (this->lcao_PI_sub_save_[ik][iat].empty()) continue;
                    double p_diag = this->lcao_PI_sub_save_[ik][iat][ib + ib * nbands].real();
                    double dl = this->lambda_[iat].z - this->lcao_lambda_in_sub_[iat].z;
                    delta_epsilon += dl * p_diag;
                }
                ekb_A[ik][ib] = this->lcao_ekb_save_[ik * nbands + ib] + spin_sign * delta_epsilon;
            }
        }
        // Update pelec ekb and recompute weights for method A
        for (int ik = 0; ik < nk; ik++)
            for (int ib = 0; ib < nbands; ib++)
                this->pelec->ekb(ik, ib) = ekb_A[ik][ib];
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);

        // Method A: first-order response does not change wavefunctions (V = identity)
        // Use full-space DMR approach: psi unchanged, cal_dm_psi → cal_DMR → cal_moment
        std::vector<std::vector<std::complex<double>>> vcc_identity(nk);
        for (int ik = 0; ik < nk; ik++) {
            vcc_identity[ik].assign(nbands * nbands, {0.0, 0.0});
            for (int ib = 0; ib < nbands; ib++)
                vcc_identity[ik][ib + ib * nbands] = {1.0, 0.0};
        }
        std::vector<ModuleBase::Vector3<double>> Mi_A = compute_mi_full_dm(vcc_identity);

        // Method B: Full subspace diagonalization
        std::vector<std::vector<double>> ekb_B(nk);
        std::vector<std::vector<std::complex<double>>> vcc_B(nk);
        for (int ik = 0; ik < nk; ik++) {
            std::vector<std::complex<double>> h_tmp(nn), s_tmp(nn);
            std::memcpy(h_tmp.data(), this->lcao_sub_h_save + ik * nn, sizeof(std::complex<double>) * nn);
            std::memcpy(s_tmp.data(), this->lcao_sub_s_save + ik * nn, sizeof(std::complex<double>) * nn);

            this->calculate_delta_hcc_lcao(h_tmp.data(), this->lcao_PI_sub_save_[ik],
                                           this->lambda_.data(), nbands, ik, true);

            std::vector<std::complex<double>> vcc(nn);
            std::vector<double> eigenvalues(nbands, 0.0);
            std::vector<std::complex<double>> s_copy(nn);
            std::memcpy(s_copy.data(), s_tmp.data(), sizeof(std::complex<double>) * nn);

            hsolver::DiagoIterAssist<std::complex<double>>::diag_hegvd(
                nbands, nbands, h_tmp.data(), s_copy.data(), nbands,
                eigenvalues.data(), vcc.data());

            vcc_B[ik].assign(vcc.data(), vcc.data() + nn);
            ekb_B[ik].assign(eigenvalues.data(), eigenvalues.data() + nbands);
        }
        // Update pelec ekb and recompute weights for method B
        for (int ik = 0; ik < nk; ik++)
            for (int ib = 0; ib < nbands; ib++)
                this->pelec->ekb(ik, ib) = ekb_B[ik][ib];
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);

        // Method B: subspace diagonalization gives V, rotate psi = C_old * V
        // Then use full-space DMR: cal_dm_psi → cal_DMR → cal_moment
        std::vector<ModuleBase::Vector3<double>> Mi_B = compute_mi_full_dm(vcc_B);

        // Method C: Full HSolverLCAO diagonalization
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
            ->update_lambda();
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

        std::vector<std::vector<double>> ekb_C(nk);
        for (int ik = 0; ik < nk; ik++) {
            ekb_C[ik].resize(nbands);
            for (int ib = 0; ib < nbands; ib++)
                ekb_C[ik][ib] = this->pelec->ekb(ik, ib);
        }
        this->cal_mi_lcao(0);
        std::vector<ModuleBase::Vector3<double>> Mi_C = this->Mi_;

        // Debug: print Mi comparison at delta=0
        if (std::abs(delta_ev) < 1e-10) {
            std::cout << "[DS-DEBUG] delta=0 Mi comparison:" << std::endl;
            for (int ia = 0; ia < nat; ia++) {
                std::cout << "  Atom " << ia << ": Mi_A=" << Mi_A[ia].z
                          << " Mi_B=" << Mi_B[ia].z << " Mi_C=" << Mi_C[ia].z << std::endl;
            }
            std::cout << "[DS-DEBUG] Pelec Mi (from cal_mi_lcao): ";
            for (int ia = 0; ia < nat; ia++) std::cout << Mi_C[ia].z << " ";
            std::cout << std::endl;
        }

        // Compute differences
        double max_dek_AB = 0.0, max_dek_AC = 0.0, max_dek_BC = 0.0;
        for (int ik = 0; ik < nk; ik++) {
            for (int ib = 0; ib < nbands; ib++) {
                double dAB = std::abs(ekb_A[ik][ib] - ekb_B[ik][ib]);
                double dAC = std::abs(ekb_A[ik][ib] - ekb_C[ik][ib]);
                double dBC = std::abs(ekb_B[ik][ib] - ekb_C[ik][ib]);
                if (dAB > max_dek_AB) max_dek_AB = dAB;
                if (dAC > max_dek_AC) max_dek_AC = dAC;
                if (dBC > max_dek_BC) max_dek_BC = dBC;
            }
        }
        double max_dMi_AB = 0.0, max_dMi_AC = 0.0, max_dMi_BC = 0.0;
        for (int ia = 0; ia < nat; ia++) {
            double dAB = std::abs(Mi_A[ia].z - Mi_B[ia].z);
            double dAC = std::abs(Mi_A[ia].z - Mi_C[ia].z);
            double dBC = std::abs(Mi_B[ia].z - Mi_C[ia].z);
            if (dAB > max_dMi_AB) max_dMi_AB = dAB;
            if (dAC > max_dMi_AC) max_dMi_AC = dAC;
            if (dBC > max_dMi_BC) max_dMi_BC = dBC;
        }

        // Output
        ofs_diag << std::scientific << std::setprecision(8);
        ofs_diag << step_count << "  " << delta_ev << "  ";
        ofs_diag << ekb_A[0][0] << "  " << ekb_B[0][0] << "  " << ekb_C[0][0] << "  |  ";
        ofs_diag << max_dek_AB << "  " << max_dek_AC << "  " << max_dek_BC << "  |  ";
        ofs_diag << Mi_A[0].z << "  " << Mi_B[0].z << "  " << Mi_C[0].z << "  |  ";
        ofs_diag << Mi_A[1].z << "  " << Mi_B[1].z << "  " << Mi_C[1].z << "  |  ";
        ofs_diag << max_dMi_AB << "  " << max_dMi_AC << "  " << max_dMi_BC << std::endl;

        std::cout << "[DS-LOCAL] delta=" << delta_ev << " eV/uB:" << std::endl;
        std::cout << "  max|dek| AB=" << max_dek_AB << " AC=" << max_dek_AC << " BC=" << max_dek_BC << std::endl;
        std::cout << "  max|dMi| AB=" << max_dMi_AB << " AC=" << max_dMi_AC << " BC=" << max_dMi_BC << std::endl;

        step_count++;
    }

    ofs_diag.close();
    std::cout << "[DS-LOCAL] Results written to: " << fname << std::endl;
    std::cout << std::string(80, '-') << "\n" << std::endl;
}

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

            // Mask unconstrained components of delta_lambda to 0.
            // For nspin=2 (collinear), constrain_[ia].x and constrain_[ia].y are
            // forced to 0 in init_sc.cpp, so only delta_lambda.z survives.
            where_fill_scalar_2d(this->constrain_, 0, zero, delta_lambda);

            // lambda = initial_lambda + delta_lambda
            add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);

            // =================================================================
            // [direction_only mode] Project out parallel component of lambda
            // =================================================================
            //
            // Physics motivation:
            //   The constraint energy is E_scon = -sum_i lambda_i . (M_i - M_target_i).
            //   The lambda vector can be decomposed into:
            //     lambda_parallel = (lambda . dir) * dir  (along target direction)
            //     lambda_perp = lambda - lambda_parallel   (perpendicular to target)
            //   - lambda_parallel controls the MAGNITUDE of M_i
            //   - lambda_perp controls the DIRECTION of M_i (rotation)
            //
            //   In direction_only mode, we want to constrain only the spin DIRECTION,
            //   not the magnitude. So we remove lambda_parallel, leaving only
            //   lambda_perp to rotate the spin toward the target direction.
            //
            // CRITICAL: This projection has a devastating effect for nspin=2 (collinear):
            //   - In collinear mode, constrain_[ia].x = constrain_[ia].y = 0 (set in
            //     init_sc.cpp), so only lambda_z is non-zero.
            //   - The target magnetization for collinear mode is also purely along z:
            //     dir = (0, 0, 1).
            //   - The parallel component is: lambda . dir = lambda_z.
            //   - After projection: lambda_z -= lambda_z * 1 = 0.
            //   RESULT: lambda becomes ZERO for all components, and the constraint
            //   is completely disabled!
            //
            // Therefore: direction_only is temporarily disabled during Phase 1
            // (sc_dir_phase1_steps iterations) for collinear calculations.
            // See esolver_ks_lcao.cpp for details.
            // =================================================================
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

        // =================================================================
        // [direction_only mode] Modify residual and adjust target magnitude
        // =================================================================
        //
        // In direction_only mode, we don't care about the MAGNITUDE error
        // (|M_i| - |M_target|), only about the DIRECTION error (angle between
        // M_i and M_target). This block:
        //
        // 1. Computes the perpendicular component of delta_spin:
        //      |delta_spin_perp|^2 = |delta_spin|^2 - (delta_spin . dir)^2
        //    This is the error that direction_only actually tries to minimize.
        //
        // 2. Adjusts target_mag by adding the parallel component:
        //      target_mag_new = target_mag + (delta_spin . dir) * dir
        //    This makes the residual parallel component zero by definition,
        //    so the BFGS optimizer doesn't try to correct the magnitude error.
        //
        // For nspin=2 (collinear):
        //   - target direction = (0, 0, 1), delta_spin = (0, 0, Mz - Mz_target)
        //   - parallel = Mz - Mz_target (the entire residual)
        //   - |delta_spin_perp|^2 = |delta_spin|^2 - parallel^2 = 0
        //   - target_mag is adjusted to include the current Mz
        //   RESULT: RMS error = 0, optimizer thinks it's converged immediately!
            //   This is another reason why direction_only is disabled during
            //   Phase 1 for collinear calculations.
        // =================================================================
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
                // Adjust target to include parallel component
                this->target_mag_[ia] += parallel * dir;
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

        // =============================================================
        // ACCELERATION ACTIVATION CHECK (LCAO nspin=2 only)
        // =============================================================
        //
        // PURPOSE: When the BFGS optimizer is near convergence (RMS error
        // drops below sc_acceleration_rms_thr), lambda changes become small.
        // At this point, subspace methods can approximate the response much
        // faster than full diagonalization while maintaining accuracy.
        //
        // ACTIVATION CONDITIONS (all must be true):
        //   1. PARAM.inp.basis_type == "lcao" (only LCAO basis supported)
        //   2. this->nspin_ == 2 (only collinear spin supported;
        //      nspin=4 requires complex-type operator template instantiation)
        //   3. sc_acceleration_mode != "off" (user must explicitly enable;
        //      valid values: "first_order", "subspace")
        //   4. sc_acceleration_rms_thr > 0 (threshold must be set)
        //   5. rms_error < sc_acceleration_rms_thr (current error is small enough)
        //
        // ONCE-ONLY ACTIVATION:
        //   The check uses `!this->acceleration_active_` to ensure subspace
        //   is built only ONCE per SCF iteration. After activation, all
        //   subsequent cal_mw_from_lambda calls use the accelerated path
        //   until the next SCF iteration resets the flag.
        //
        // WHY build subspace at current lambda?
        //   The subspace approximation is only valid for small perturbations
        //   around the reference lambda. Building it at the current (near-converged)
        //   lambda ensures the reference point is close to the optimal solution.
        //   If we built it at lambda=0, the approximation would be poor for
        //   large lambda values.
        //
        // WHY use i_step=-2?
        //   The main BFGS loop uses i_step = -1, 0, 1, ..., nsc-1.
        //   Using -2 ensures this call is handled by a special branch in
        //   cal_mw_from_lambda that does full diagonalization + cache build
        //   WITHOUT affecting the main loop state (psi, delta_lambda, etc.).
        //   This is a clean separation of concerns.
        //
        // WHY free_lcao_subspace_cache() first?
        //   To prevent memory leaks from previous SCF iterations. Each SCF
        //   iteration should build a fresh subspace at its own converged lambda.
        // =============================================================
        const bool accel_enabled = (PARAM.inp.basis_type == "lcao") &&
                                   (this->nspin_ == 2) &&
                                   (this->sc_acceleration_mode_ != "off") &&
                                   (this->sc_acceleration_rms_thr_ > 0.0) &&
                                   (rms_error < this->sc_acceleration_rms_thr_);

        if (accel_enabled && !this->acceleration_active_)
        {
            // First time crossing threshold: activate and build subspace at current lambda
            this->acceleration_active_ = true;
            this->acceleration_subspace_built_ = false;
            this->lambda_at_acceleration_ = this->lambda_;

            std::cout << "\n[DS-ACCEL] RMS = " << rms_error << " uB < threshold " 
                      << this->sc_acceleration_rms_thr_ << " uB" << std::endl;
            std::cout << "[DS-ACCEL] Acceleration activated: mode = " << this->sc_acceleration_mode_ << std::endl;
            std::cout << "[DS-ACCEL] Building subspace at current lambda...\n" << std::endl;

            // Free any cached subspace data from previous SCF iterations
            this->free_lcao_subspace_cache();

            // Build fresh subspace at current lambda using full diagonalization.
            // i_step=-2 triggers the special cache-build branch in cal_mw_from_lambda.
            // This does NOT modify the BFGS loop state (psi, delta_lambda).
            this->cal_mw_from_lambda(-2, delta_lambda.data());

            this->acceleration_subspace_built_ = true;
        }

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

            // [LCAO nspin=2] Run local diagnostic near converged lambda
            // Trigger when RMS is small (lambda near converged), regardless of drho
            // BFGS changes Hamiltonian each step, so drho may stay large
#ifdef __LCAO
            if (PARAM.inp.basis_type == "lcao" && this->nspin_ == 2
                && rms_error < 5.0 && !this->local_diag_run_)
            {
                // Use the current converged lambda as reference
                double lambda_ref_ry = 0.0;
                for (int ia = 0; ia < nat; ia++) {
                    if (this->constrain_[ia].z != 0) {
                        lambda_ref_ry = this->lambda_[ia].z;
                        break;
                    }
                }
                std::cout << "\n[DS-DIAG] === LOCAL DIAGNOSTIC TRIGGERED (BFGS converged) ===" << std::endl;
                std::cout << "[DS-DIAG] RMS = " << rms_error << " uB (< 0.5 threshold)" << std::endl;
                std::cout << "[DS-DIAG] Converged lambda: " << lambda_ref_ry * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;

                // Test: Build subspace at converged lambda, scan +/-0.5 eV/uB
                this->run_lambda_local_diagnostic(outer_step, lambda_ref_ry, "at_converged");

                // Comparison: Build subspace at lambda=0, scan +/-0.5 eV/uB
                this->run_lambda_local_diagnostic(outer_step, 0.0, "at_zero");

                this->local_diag_run_ = true;
                std::cout << "[DS-DIAG] Local diagnostic complete.\n" << std::endl;
            }
#endif

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
        // dnu is the cumulative search path integral: the total lambda change
        // accumulated over all inner BFGS steps. It starts at 0 and grows as:
        //   dnu += alpha_trial * search  (at each inner step)
        //
        // dnu is what gets applied to lambda at the end of the loop:
        //   lambda = initial_lambda + dnu
        dnu_last_step = dnu;
        // dnu = dnu + alpha_trial * search
        add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);

        // [direction_only] Project out parallel component from dnu.
        // This prevents any accumulation of the parallel (magnitude-controlling)
        // component in the cumulative step. Without this projection, small
        // numerical errors in the parallel component could accumulate over
        // many inner steps, eventually affecting the magnitude.
        //
        // For nspin=2 (collinear): this zeroes dnu.z (the only non-zero component),
        // making dnu = 0. This is why direction_only is disabled during Phase 1 for collinear calculations.
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++) {
            const auto& target = this->target_mag_[ia];
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

        // [direction_only] Project out parallel component from corrected dnu.
        // Same as above: after the optimal step correction, remove any parallel
        // component that may have been introduced.
        // For nspin=2 collinear: this zeroes the only non-zero component (dnu.z).
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++) {
            const auto& target = this->target_mag_[ia];
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

    // Track optimal lambda (minimum RMS to target) for local diagnostic
    double min_rms = std::numeric_limits<double>::max();
    double optimal_lambda_ry = 0.0;

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

        // At step 0 of each SCF iteration, force full diagonalization to refresh
        // subspace cache with current wavefunctions. This avoids using stale
        // subspace data from the previous SCF iteration.
        if (istep == 0) {
            this->free_lcao_subspace_cache();
        }

        // Compute magnetic moments using eigenvalue response method (current path)
        // For istep==0: does full diagonalization and rebuilds subspace cache
        // For istep>0: uses first-order eigenvalue shift with cached subspace data
        this->cal_mw_from_lambda(istep);

        // Save response method Mi
        std::vector<ModuleBase::Vector3<double>> Mi_response = this->Mi_;

        // Restore response Mi as the official result
        this->Mi_ = Mi_response;

        // Compute RMS to target and track optimal lambda
        double rms = 0.0;
        for (int ia = 0; ia < nat; ia++) {
            ModuleBase::Vector3<double> diff = this->Mi_[ia] - this->target_mag_[ia];
            rms += diff * diff;
        }
        rms = std::sqrt(rms / nat);
        if (rms < min_rms) {
            min_rms = rms;
            optimal_lambda_ry = lambda_val_ry;
        }

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

    // Restore original lambda values
    this->lambda_ = initial_lambda;

    // =============================================================
    // LOCAL DIAGNOSTIC: compare methods near optimal lambda
    // Only runs once when charge density is near convergence (drho < 1e-3)
    // =============================================================
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao" && this->nspin_ == 2
        && this->last_drho_ > 0 && this->last_drho_ < 1e-3
        && !this->local_diag_run_)
    {
        std::cout << "\n[DS-DIAG] === LOCAL DIAGNOSTIC TRIGGERED ===" << std::endl;
        std::cout << "[DS-DIAG] drho = " << this->last_drho_ << " (< 1e-3)" << std::endl;
        std::cout << "[DS-DIAG] Optimal lambda from scan: " << optimal_lambda_ry * ModuleBase::Ry_to_eV << " eV/uB (RMS=" << min_rms << " uB)" << std::endl;

        // Test 1: Build subspace at optimal lambda (converged lambda), scan ±0.5 eV/uB
        this->run_lambda_local_diagnostic(outer_step, optimal_lambda_ry, "at_optimal");

        // Test 2: Build subspace at lambda=0, scan ±0.5 eV/uB (for comparison)
        this->run_lambda_local_diagnostic(outer_step, 0.0, "at_zero");

        this->local_diag_run_ = true;
        std::cout << "[DS-DIAG] Local diagnostic complete.\n" << std::endl;
    }
#endif

    // =============================================================
    // DIAGNOSTIC: Compare first-order response vs full subspace diag vs full diag
    // Only runs when charge density is near convergence (drho < 1e-3)
    // =============================================================
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao" && this->nspin_ == 2 && this->last_drho_ > 0 && this->last_drho_ < 1e-3)
    {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "[DS-DIAG] === EIGENVALUE & Mi COMPARISON DIAGNOSTIC ===" << std::endl;
        std::cout << "[DS-DIAG] Methods: (A) First-order response, (B) Full subspace diag, (C) Full HSolverLCAO" << std::endl;
        std::cout << "[DS-DIAG] SCF convergence: drho = " << this->last_drho_ << " (< 1e-3 threshold)" << std::endl;
        std::cout << std::string(80, '=') << "\n" << std::endl;

        // Force full diagonalization to build fresh subspace cache at lambda=0
        this->free_lcao_subspace_cache();
        this->lambda_ = std::vector<ModuleBase::Vector3<double>>(nat, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));

        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);

        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
            ->update_lambda();

        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

        // Cache subspace data
        const int nk = psi_t->get_nk();
        const int nbands = PARAM.inp.nbands;
        const int nlocal = this->ParaV->nrow;
        const int nn = nbands * nbands;
        this->lcao_nbands_ = nbands;
        this->lcao_nk_ = nk;
        this->lcao_nlocal_ = nlocal;

        delete[] this->lcao_sub_h_save;
        delete[] this->lcao_sub_s_save;
        this->lcao_sub_h_save = new std::complex<double>[nk * nn];
        this->lcao_sub_s_save = new std::complex<double>[nk * nn];
        this->lcao_PI_sub_save_.resize(nk);
        this->lcao_ekb_save_.resize(nk * nbands);

        for (int ik = 0; ik < nk; ik++)
        {
            psi_t->fix_k(ik);
            this->calculate_lcao_sub_hs(
                this->p_hamilt, psi_t[0], this->ParaV,
                this->lcao_sub_h_save + ik * nn,
                this->lcao_sub_s_save + ik * nn,
                ik, nbands, nlocal);

            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
                this->p_operator);
            dspin_op->cal_PI_sub(
                this->kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                this->lcao_PI_sub_save_[ik]);

            for (int ib = 0; ib < nbands; ib++)
                this->lcao_ekb_save_[ik * nbands + ib] = this->pelec->ekb(ik, ib);
        }
        this->lcao_lambda_in_sub_ = this->lambda_;
        this->lcao_subspace_initialized_ = true;

        // Save reference ekb at lambda=0
        std::vector<std::vector<double>> ekb_ref(nk);
        for (int ik = 0; ik < nk; ik++) {
            ekb_ref[ik].resize(nbands);
            for (int ib = 0; ib < nbands; ib++)
                ekb_ref[ik][ib] = this->pelec->ekb(ik, ib);
        }

        // Open comparison file
        std::ofstream ofs_diag("eigenvalue_mi_comparison.dat");
        ofs_diag << "# Eigenvalue and Magnetic Moment Comparison" << std::endl;
        ofs_diag << "# Methods: (A) First-order response, (B) Full subspace diag, (C) Full HSolverLCAO" << std::endl;
        ofs_diag << "# lambda range: " << lambda_start << " -> " << lambda_end << " eV/uB" << std::endl;
        ofs_diag << "# nsteps = " << nsteps << std::endl;
        ofs_diag << "#" << std::endl;
        ofs_diag << "# step  lambda_eV  ";
        // Eigenvalue section
        ofs_diag << "ekb_A_ib0  ekb_B_ib0  ekb_C_ib0  |  ekb_A_ib1  ekb_B_ib1  ekb_C_ib1  |  ";
        ofs_diag << "max_dek_AB  max_dek_AC  max_dek_BC  |  ";
        // Mi section
        ofs_diag << "Mi_A_0  Mi_B_0  Mi_C_0  |  Mi_A_1  Mi_B_1  Mi_C_1  |  ";
        ofs_diag << "max_dMi_AB  max_dMi_AC  max_dMi_BC" << std::endl;

        double lambda_start_ry = lambda_start / ModuleBase::Ry_to_eV;
        double lambda_end_ry = lambda_end / ModuleBase::Ry_to_eV;
        double lambda_step_diag = (lambda_end_ry - lambda_start_ry) / (nsteps - 1);

        for (int istep = 0; istep < nsteps; istep++)
        {
            double lambda_val_ry = lambda_start_ry + istep * lambda_step_diag;
            double lambda_val_ev = lambda_val_ry * ModuleBase::Ry_to_eV;

            // Set lambda
            for (int ia = 0; ia < nat; ia++) {
                for (int ic = 0; ic < 3; ic++) {
                    if (this->constrain_[ia][ic] != 0)
                        this->lambda_[ia][ic] = lambda_val_ry;
                    else
                        this->lambda_[ia][ic] = 0.0;
                }
            }

            // ========================================================
            // Method A: First-order eigenvalue response
            // ========================================================
            std::vector<std::vector<double>> ekb_A(nk);
            for (int ik = 0; ik < nk; ik++) {
                ekb_A[ik].resize(nbands);
                int spin_sign = this->get_spin_sign(ik);
                for (int ib = 0; ib < nbands; ib++) {
                    double delta_epsilon = 0.0;
                    for (int iat = 0; iat < this->get_nat(); iat++) {
                        if (this->lcao_PI_sub_save_[ik][iat].empty()) continue;
                        double p_diag = this->lcao_PI_sub_save_[ik][iat][ib + ib * nbands].real();
                        double dl = this->lambda_[iat].z - this->lcao_lambda_in_sub_[iat].z;
                        delta_epsilon += dl * p_diag;
                    }
                    ekb_A[ik][ib] = this->lcao_ekb_save_[ik * nbands + ib] + spin_sign * delta_epsilon;
                }
            }
            // Compute Mi_A using trace formula with identity V
            std::vector<std::vector<std::complex<double>>> vcc_identity(nk);
            for (int ik = 0; ik < nk; ik++) {
                vcc_identity[ik].assign(nbands * nbands, {0.0, 0.0});
                for (int ib = 0; ib < nbands; ib++)
                    vcc_identity[ik][ib + ib * nbands] = {1.0, 0.0};
            }
            this->cal_mi_lcao_subspace(
                this->lcao_PI_sub_save_, vcc_identity, ekb_A, nbands, nk, this->npol_);
            std::vector<ModuleBase::Vector3<double>> Mi_A = this->Mi_;

            // ========================================================
            // Method B: Full subspace diagonalization
            // ========================================================
            std::vector<std::vector<double>> ekb_B(nk);
            std::vector<std::vector<std::complex<double>>> vcc_B(nk);
            for (int ik = 0; ik < nk; ik++) {
                std::vector<std::complex<double>> h_tmp(nn), s_tmp(nn);
                std::memcpy(h_tmp.data(), this->lcao_sub_h_save + ik * nn, sizeof(std::complex<double>) * nn);
                std::memcpy(s_tmp.data(), this->lcao_sub_s_save + ik * nn, sizeof(std::complex<double>) * nn);

                this->calculate_delta_hcc_lcao(h_tmp.data(), this->lcao_PI_sub_save_[ik],
                                               this->lambda_.data(), nbands, ik, true);

                std::vector<std::complex<double>> vcc(nn);
                std::vector<double> eigenvalues(nbands, 0.0);
                std::vector<std::complex<double>> s_copy(nn);
                std::memcpy(s_copy.data(), s_tmp.data(), sizeof(std::complex<double>) * nn);

                hsolver::DiagoIterAssist<std::complex<double>>::diag_hegvd(
                    nbands, nbands, h_tmp.data(), s_copy.data(), nbands,
                    eigenvalues.data(), vcc.data());

                vcc_B[ik].assign(vcc.data(), vcc.data() + nn);
                ekb_B[ik].assign(eigenvalues.data(), eigenvalues.data() + nbands);
            }
            this->cal_mi_lcao_subspace(
                this->lcao_PI_sub_save_, vcc_B, ekb_B, nbands, nk, this->npol_);
            std::vector<ModuleBase::Vector3<double>> Mi_B = this->Mi_;

            // ========================================================
            // Method C: Full HSolverLCAO diagonalization
            // ========================================================
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
            hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
            elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                         this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                         this->pelec->skip_weights);
            elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

            std::vector<std::vector<double>> ekb_C(nk);
            for (int ik = 0; ik < nk; ik++) {
                ekb_C[ik].resize(nbands);
                for (int ib = 0; ib < nbands; ib++)
                    ekb_C[ik][ib] = this->pelec->ekb(ik, ib);
            }
            this->cal_mi_lcao(0);
            std::vector<ModuleBase::Vector3<double>> Mi_C = this->Mi_;

            // ========================================================
            // Compute differences
            // ========================================================
            double max_dek_AB = 0.0, max_dek_AC = 0.0, max_dek_BC = 0.0;
            for (int ik = 0; ik < nk; ik++) {
                for (int ib = 0; ib < nbands; ib++) {
                    double dAB = std::abs(ekb_A[ik][ib] - ekb_B[ik][ib]);
                    double dAC = std::abs(ekb_A[ik][ib] - ekb_C[ik][ib]);
                    double dBC = std::abs(ekb_B[ik][ib] - ekb_C[ik][ib]);
                    if (dAB > max_dek_AB) max_dek_AB = dAB;
                    if (dAC > max_dek_AC) max_dek_AC = dAC;
                    if (dBC > max_dek_BC) max_dek_BC = dBC;
                }
            }
            double max_dMi_AB = 0.0, max_dMi_AC = 0.0, max_dMi_BC = 0.0;
            for (int ia = 0; ia < nat; ia++) {
                double dAB = std::abs(Mi_A[ia].z - Mi_B[ia].z);
                double dAC = std::abs(Mi_A[ia].z - Mi_C[ia].z);
                double dBC = std::abs(Mi_B[ia].z - Mi_C[ia].z);
                if (dAB > max_dMi_AB) max_dMi_AB = dAB;
                if (dAC > max_dMi_AC) max_dMi_AC = dAC;
                if (dBC > max_dMi_BC) max_dMi_BC = dBC;
            }

            // Output
            ofs_diag << std::scientific << std::setprecision(8);
            ofs_diag << istep << "  " << lambda_val_ev << "  ";
            // ekb for first 2 bands (ik=0)
            ofs_diag << ekb_A[0][0] << "  " << ekb_B[0][0] << "  " << ekb_C[0][0] << "  |  ";
            ofs_diag << ekb_A[0][1] << "  " << ekb_B[0][1] << "  " << ekb_C[0][1] << "  |  ";
            ofs_diag << max_dek_AB << "  " << max_dek_AC << "  " << max_dek_BC << "  |  ";
            // Mi
            ofs_diag << Mi_A[0].z << "  " << Mi_B[0].z << "  " << Mi_C[0].z << "  |  ";
            ofs_diag << Mi_A[1].z << "  " << Mi_B[1].z << "  " << Mi_C[1].z << "  |  ";
            ofs_diag << max_dMi_AB << "  " << max_dMi_AC << "  " << max_dMi_BC << std::endl;

            std::cout << "[DS-DIAG] step " << istep << " lambda=" << lambda_val_ev << " eV/uB:" << std::endl;
            std::cout << "  ekb[0] A=" << ekb_A[0][0] << " B=" << ekb_B[0][0] << " C=" << ekb_C[0][0] << std::endl;
            std::cout << "  max|dek| AB=" << max_dek_AB << " AC=" << max_dek_AC << " BC=" << max_dek_BC << std::endl;
            std::cout << "  Mi[0]  A=" << Mi_A[0].z << " B=" << Mi_B[0].z << " C=" << Mi_C[0].z << std::endl;
            std::cout << "  max|dMi| AB=" << max_dMi_AB << " AC=" << max_dMi_AC << " BC=" << max_dMi_BC << std::endl;
        }

        ofs_diag.close();
        std::cout << "\n[DS-DIAG] Comparison written to: eigenvalue_mi_comparison.dat" << std::endl;
        std::cout << std::string(80, '=') << "\n" << std::endl;
    }
#endif

    std::cout << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === LINEAR LAMBDA SCAN COMPLETE ===" << std::endl;
    std::cout << "[DS-DIAG] Results written to: lambda_scan_results.dat" << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;

    return;
}

/**
 * @brief Diagnostic: compare subspace Mi vs full diagonalization Mi at various lambda values.
 *
 * @details For nspin=2 LCAO only. At each lambda point:
 *   1. Compute subspace Mi using trace formula
 *   2. Compute full Mi using full diagonalization
 *   3. Output both and relative error
 *
 * Results written to subspace_vs_full_scan.dat
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_scan_diagnostic(int outer_step)
{
#ifdef __LCAO
    if (PARAM.inp.basis_type != "lcao")
    {
        std::cout << "[DS-DIAG] scan_diagnostic: only supported for LCAO basis" << std::endl;
        return;
    }
    if (this->nspin_ != 2)
    {
        std::cout << "[DS-DIAG] scan_diagnostic: only supported for nspin=2" << std::endl;
        return;
    }
#else
    return;
#endif

    int nat = this->get_nat();

    double lambda_start = PARAM.inp.sc_scan_lambda_start;
    double lambda_end = PARAM.inp.sc_scan_lambda_end;
    int nsteps = PARAM.inp.sc_scan_steps;

    if (nsteps <= 0) {
        std::cout << "[DS-DIAG] scan_diagnostic: sc_scan_steps <= 0, skipping" << std::endl;
        return;
    }

    double lambda_start_ry = lambda_start / ModuleBase::Ry_to_eV;
    double lambda_end_ry = lambda_end / ModuleBase::Ry_to_eV;
    double lambda_step = (lambda_end_ry - lambda_start_ry) / (nsteps - 1);

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === SUBSPACE vs FULL DIAGNOSTIC SCAN START ===" << std::endl;
    std::cout << "[DS-DIAG] Scan range: " << lambda_start << " -> " << lambda_end << " eV/uB" << std::endl;
    std::cout << "[DS-DIAG] Number of steps: " << nsteps << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;

    // Save initial state
    std::vector<ModuleBase::Vector3<double>> initial_lambda = this->lambda_;

    // Open output file
    std::ofstream ofs_scan;
    ofs_scan.open("subspace_vs_full_scan.dat");
    ofs_scan << "# Subspace vs Full Diagonalization Diagnostic Scan" << std::endl;
    ofs_scan << "# lambda range: " << lambda_start << " -> " << lambda_end << " eV/uB" << std::endl;
    ofs_scan << "# nsteps = " << nsteps << std::endl;
    ofs_scan << "#" << std::endl;
    ofs_scan << "# step  lambda_eV  ";
    for (int ia = 0; ia < nat; ia++) {
        ofs_scan << "sub_Mz_" << ia << "  full_Mz_" << ia << "  rel_err_" << ia << "  ";
    }
    ofs_scan << std::endl;

    // Step 0: full diagonalization to initialize subspace cache
    std::cout << "[DS-DIAG] Step 0: Initializing subspace cache with full diagonalization..." << std::endl;
    this->lambda_ = std::vector<ModuleBase::Vector3<double>>(nat, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    
    // Force full diagonalization path by resetting subspace
    this->free_lcao_subspace_cache();

    // Call cal_mw_from_lambda with i_step=-1 to do full diag and cache
    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
    hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);

    if (this->nspin_ == 2)
    {
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
            ->update_lambda();
    }

    hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
    hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
    elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                 this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                 this->pelec->skip_weights);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    // Cache subspace data
    const int nk = psi_t->get_nk();
    const int nbands = PARAM.inp.nbands;
    const int nlocal = this->ParaV->nrow;
    int nn = nbands * nbands;
    this->lcao_nbands_ = nbands;
    this->lcao_nk_ = nk;
    this->lcao_nlocal_ = nlocal;

    delete[] this->lcao_sub_h_save;
    delete[] this->lcao_sub_s_save;
    this->lcao_sub_h_save = new std::complex<double>[nk * nn];
    this->lcao_sub_s_save = new std::complex<double>[nk * nn];
    this->lcao_PI_sub_save_.resize(nk);
    this->lcao_ekb_save_.resize(nk * nbands);

    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);
        this->calculate_lcao_sub_hs(
            this->p_hamilt, psi_t[0], this->ParaV,
            this->lcao_sub_h_save + ik * nn,
            this->lcao_sub_s_save + ik * nn,
            ik, nbands, nlocal);

        this->lcao_PI_sub_save_[ik].resize(this->get_nat());
        if (this->nspin_ == 2)
        {
            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
                this->p_operator);
            dspin_op->cal_PI_sub(
                this->kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                this->lcao_PI_sub_save_[ik]);
        }

        for (int ib = 0; ib < nbands; ib++)
        {
            this->lcao_ekb_save_[ik * nbands + ib] = this->pelec->ekb(ik, ib);
        }
    }

    this->lcao_lambda_in_sub_ = this->lambda_;
    this->lcao_subspace_initialized_ = true;

    // Get full Mi at lambda=0
    this->cal_mi_lcao(0);
    std::vector<ModuleBase::Vector3<double>> Mi_full_init = this->Mi_;

    std::cout << "[DS-DIAG] Full Mi at lambda=0: ";
    for (int ia = 0; ia < nat; ia++)
        std::cout << Mi_full_init[ia].z << " ";
    std::cout << std::endl;

    // =============================================================
    // SCAN LOOP
    // =============================================================
    for (int istep = 0; istep < nsteps; istep++)
    {
        double lambda_val_ry = lambda_start_ry + istep * lambda_step;
        double lambda_val_ev = lambda_val_ry * ModuleBase::Ry_to_eV;

        // Set lambda
        for (int ia = 0; ia < nat; ia++)
        {
            if (this->constrain_[ia].z != 0)
                this->lambda_[ia].z = lambda_val_ry;
            else
                this->lambda_[ia].z = 0.0;
        }

        std::cout << "\n[DS-DIAG] === Step " << istep << "/" << (nsteps-1)
                  << " lambda = " << lambda_val_ev << " eV/uB ===" << std::endl;

        // --- Subspace Mi ---
        // Compute subspace Hamiltonian and diagonalize
        std::vector<std::vector<std::complex<double>>> vcc_all(nk);
        std::vector<std::vector<double>> ekb_all(nk);

        for (int ik = 0; ik < nk; ik++)
        {
            std::vector<std::complex<double>> h_tmp(nn), s_tmp(nn);
            std::memcpy(h_tmp.data(), this->lcao_sub_h_save + ik * nn, sizeof(std::complex<double>) * nn);
            std::memcpy(s_tmp.data(), this->lcao_sub_s_save + ik * nn, sizeof(std::complex<double>) * nn);

            this->calculate_delta_hcc_lcao(h_tmp.data(), this->lcao_PI_sub_save_[ik],
                                           this->lambda_.data(), nbands, ik, true);

            std::vector<std::complex<double>> vcc(nn);
            std::vector<double> eigenvalues(nbands, 0.0);
            std::vector<std::complex<double>> s_copy(nn);
            std::memcpy(s_copy.data(), s_tmp.data(), sizeof(std::complex<double>) * nn);

            hsolver::DiagoIterAssist<std::complex<double>>::diag_hegvd(
                nbands, nbands, h_tmp.data(), s_copy.data(), nbands,
                eigenvalues.data(), vcc.data());

            vcc_all[ik].assign(vcc.data(), vcc.data() + nn);
            ekb_all[ik] = eigenvalues;

            for (int ib = 0; ib < nbands; ib++)
            {
                this->pelec->ekb(ik, ib) = eigenvalues[ib];
            }
        }

        // Calculate weights and subspace Mi
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

        this->cal_mi_lcao_subspace(
            this->lcao_PI_sub_save_, vcc_all, ekb_all, nbands, nk, this->npol_);

        std::vector<ModuleBase::Vector3<double>> Mi_subspace = this->Mi_;

        std::cout << "[DS-DIAG] Subspace Mi: ";
        for (int ia = 0; ia < nat; ia++)
            std::cout << Mi_subspace[ia].z << " ";
        std::cout << std::endl;

        // --- Full Mi ---
        // Full diagonalization
        if (this->nspin_ == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
        }

        hsolver::HSolverLCAO<std::complex<double>> hsolver_full(this->ParaV, PARAM.inp.ks_solver);
        hsolver_full.solve(hamilt_t, psi_t[0], this->pelec, *this->dm_, *this->pelec->charge, this->nspin_, true);
        elecstate::calculate_weights(this->pelec->ekb, this->pelec->wg, this->pelec->klist,
                                     this->pelec->eferm, this->pelec->f_en, this->pelec->nelec_spin,
                                     this->pelec->skip_weights);
        elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

        this->cal_mi_lcao(0);

        std::vector<ModuleBase::Vector3<double>> Mi_full = this->Mi_;

        std::cout << "[DS-DIAG] Full Mi:      ";
        for (int ia = 0; ia < nat; ia++)
            std::cout << Mi_full[ia].z << " ";
        std::cout << std::endl;

        // --- Compare ---
        std::cout << "[DS-DIAG] Rel error(%): ";
        ofs_scan << std::scientific << std::setprecision(6);
        ofs_scan << istep << "  " << lambda_val_ev;

        double max_rel_err = 0;
        for (int ia = 0; ia < nat; ia++)
        {
            double rel_err = 0;
            if (std::abs(Mi_full[ia].z) > 1e-10)
            {
                rel_err = std::abs((Mi_subspace[ia].z - Mi_full[ia].z) / Mi_full[ia].z) * 100.0;
                if (rel_err > max_rel_err) max_rel_err = rel_err;
            }
            std::cout << rel_err << " ";
            ofs_scan << "  " << Mi_subspace[ia].z << "  " << Mi_full[ia].z << "  " << rel_err;
        }
        std::cout << "  (max: " << max_rel_err << "%)" << std::endl;
        ofs_scan << std::endl;
    }

    ofs_scan.close();

    // Restore initial state
    this->lambda_ = initial_lambda;
    this->free_lcao_subspace_cache();

    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === DIAGNOSTIC SCAN COMPLETE ===" << std::endl;
    std::cout << "[DS-DIAG] Results written to: subspace_vs_full_scan.dat" << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;
}


template <>
void spinconstrain::SpinConstrain<double>::run_lambda_scan_diagnostic(int outer_step)
{
    std::cout << "[DS-DIAG] scan_diagnostic: only implemented for complex<double> (nspin=2)" << std::endl;
}

template <>
void spinconstrain::SpinConstrain<double>::run_lambda_local_diagnostic(
    int outer_step, double lambda_ref_ry, const std::string& label)
{
    std::cout << "[DS-LOCAL] local_diagnostic: only implemented for complex<double> (nspin=2)" << std::endl;
}

#include "spin_constrain.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "basic_funcs.h"
#include "source_base/constants.h"
#include "source_io/module_parameter/parameter.h"

/**
 * @file lambda_loop_diag.cpp
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
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_linear_scan(int outer_step, std::ostream& ofs_running)
{
    int nat = this->get_nat();
    int ntype = this->get_ntype();

    double lambda_start = PARAM.inp.sc_scan_lambda_start;
    double lambda_end = PARAM.inp.sc_scan_lambda_end;
    int nsteps = PARAM.inp.sc_scan_steps;

    if (nsteps <= 0) {
        ofs_running << " [DS-DIAG] linear_scan: sc_scan_steps <= 0, skipping" << std::endl;
        return;
    }

    // Convert eV to Ry for internal calculations
    double lambda_start_ry = lambda_start / ModuleBase::Ry_to_eV;
    double lambda_end_ry = lambda_end / ModuleBase::Ry_to_eV;
    double lambda_step = (lambda_end_ry - lambda_start_ry) / (nsteps - 1);

    ofs_running << "\n" << std::string(80, '=') << std::endl;
    ofs_running << " [DS-DIAG] === LINEAR LAMBDA SCAN START ===" << std::endl;
    ofs_running << " [DS-DIAG] Scan range: " << lambda_start << " -> " << lambda_end << " eV/uB" << std::endl;
    ofs_running << " [DS-DIAG] Number of steps: " << nsteps << std::endl;
    ofs_running << " [DS-DIAG] Lambda step size: " << lambda_step * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;
    ofs_running << " [DS-DIAG] nat = " << nat << ", ntype = " << ntype << std::endl;
    ofs_running << " [DS-DIAG] nspin_ = " << this->state_.nspin_ << ", npol_ = " << this->state_.npol_ << std::endl;
    ofs_running << " [DS-DIAG] p_operator = " << (this->p_operator ? "valid" : "NULL") << std::endl;
    ofs_running << " [DS-DIAG] constrain_ size = " << this->state_.constrain_.size() << std::endl;

    // Check if any constraints are defined; if not, set all atoms as constrained
    bool has_constraints = false;
    for (int ia = 0; ia < nat; ia++) {
        if (this->state_.constrain_[ia].x != 0 || this->state_.constrain_[ia].y != 0 || this->state_.constrain_[ia].z != 0) {
            has_constraints = true;
            break;
        }
    }

    if (!has_constraints) {
        ofs_running << " [DS-DIAG] No constraints found in STRU, setting all atoms as constrained" << std::endl;
        for (int ia = 0; ia < nat; ia++) {
            if (this->state_.nspin_ == 4) {
                this->state_.constrain_[ia] = ModuleBase::Vector3<int>(1, 1, 1);
            } else {
                this->state_.constrain_[ia] = ModuleBase::Vector3<int>(0, 0, 1);
            }
        }
        this->reset_dspin_operator();
    }

    for (int ia = 0; ia < nat; ia++) {
        ofs_running << " [DS-DIAG]   Atom " << ia << " constrain = ("
                             << this->state_.constrain_[ia].x << ", " << this->state_.constrain_[ia].y << ", " << this->state_.constrain_[ia].z << ")"
                             << " target_mag = (" << this->state_.target_mag_[ia].x << ", " << this->state_.target_mag_[ia].y << ", " << this->state_.target_mag_[ia].z << ")" << std::endl;
    }
    ofs_running << std::string(80, '=') << "\n" << std::endl;

    // Save initial lambda to restore after scan
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    where_fill_scalar_else_2d(this->state_.constrain_, 0, 0.0, this->state_.get_sc_lambda(), initial_lambda);

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

    double original_sc_thr = this->state_.sc_thr_;

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
                if (this->state_.constrain_[ia][ic] != 0) {
                    this->state_.get_lambda()[ia][ic] = lambda_val_ry;
                } else {
                    this->state_.get_lambda()[ia][ic] = 0.0;
                }
            }
        }

        ofs_running << " [DS-DIAG] === Scan step " << istep << "/" << nsteps
                             << " lambda = " << lambda_val_ev << " eV/uB ===" << std::endl;

        // Compute magnetic moments at current lambda
        this->cal_mw_from_lambda(istep);

        // Save step 0 Mi for consistency verification
        if (istep == 0) {
            mi_step0 = this->state_.get_mi();
        }

        // Write results
        ofs_scan << std::scientific << std::setprecision(6);
        ofs_scan << istep << "  " << lambda_val_ev;
        for (int ia = 0; ia < nat; ia++) {
            ofs_scan << "  " << this->state_.get_mi()[ia].x
                     << "  " << this->state_.get_mi()[ia].y
                     << "  " << this->state_.get_mi()[ia].z;
        }
        ofs_scan << std::endl;

        ofs_running << " [DS-DIAG]   lambda = " << lambda_val_ev << " eV/uB" << std::endl;
        for (int ia = 0; ia < nat; ia++) {
            ofs_running << " [DS-DIAG]   Atom " << ia << " Mi = ("
                                 << this->state_.get_mi()[ia].x << ", "
                                 << this->state_.get_mi()[ia].y << ", "
                                 << this->state_.get_mi()[ia].z << ") uB" << std::endl;
        }
        ofs_running << std::endl;
    }

    // =============================================================
    // CONSISTENCY CHECK: restore initial lambda and recompute Mi
    // to verify that the lambda->Mi mapping is numerically stable
    // after multiple lambda updates in the scan loop
    // =============================================================
    ofs_running << " [DS-DIAG] === Consistency check: restoring initial lambda ===" << std::endl;
    this->state_.get_lambda() = initial_lambda;
    this->cal_mw_from_lambda(nsteps);

    // Write consistency check result
    ofs_scan << std::scientific << std::setprecision(6);
    ofs_scan << "init_recheck  " << lambda_start;
    for (int ia = 0; ia < nat; ia++) {
        ofs_scan << "  " << this->state_.get_mi()[ia].x
                 << "  " << this->state_.get_mi()[ia].y
                 << "  " << this->state_.get_mi()[ia].z;
    }
    ofs_scan << std::endl;

    ofs_running << " [DS-DIAG]   lambda = " << lambda_start << " eV/uB (restored)" << std::endl;
    for (int ia = 0; ia < nat; ia++) {
        ofs_running << " [DS-DIAG]   Atom " << ia << " Mi = ("
                             << this->state_.get_mi()[ia].x << ", "
                             << this->state_.get_mi()[ia].y << ", "
                             << this->state_.get_mi()[ia].z << ") uB" << std::endl;
    }

    // Compare restored Mi with step 0 Mi to check consistency
    ofs_scan << "# [consistency] step 0 vs init_recheck Mi difference:" << std::endl;
    double max_mi_diff = 0.0;
    for (int ia = 0; ia < nat; ia++) {
        double dx = std::abs(this->state_.get_mi()[ia].x - mi_step0[ia].x);
        double dy = std::abs(this->state_.get_mi()[ia].y - mi_step0[ia].y);
        double dz = std::abs(this->state_.get_mi()[ia].z - mi_step0[ia].z);
        double diff = std::max({dx, dy, dz});
        if (diff > max_mi_diff) max_mi_diff = diff;
        ofs_scan << "#   Atom " << ia << " dM = (" << dx << ", " << dy << ", " << dz << ") uB" << std::endl;
    }
    ofs_running << " [DS-DIAG] Max Mi difference between step 0 and init_recheck: " << max_mi_diff << " uB" << std::endl;
    if (max_mi_diff > 1e-8) {
        ofs_running << " [DS-DIAG] WARNING: Mi mapping may be inconsistent after multiple lambda updates!" << std::endl;
    } else {
        ofs_running << " [DS-DIAG] OK: Mi mapping is consistent." << std::endl;
    }
    ofs_scan << "#   Max Mi difference: " << max_mi_diff << " uB" << std::endl;

    ofs_scan.close();

    // Restore original lambda values (already restored above, but explicit for clarity)
    this->state_.get_lambda() = initial_lambda;

    ofs_running << std::string(80, '=') << std::endl;
    ofs_running << " [DS-DIAG] === LINEAR LAMBDA SCAN COMPLETE ===" << std::endl;
    ofs_running << " [DS-DIAG] Results written to: lambda_scan_results.dat" << std::endl;
    ofs_running << std::string(80, '=') << "\n" << std::endl;

    return;
}

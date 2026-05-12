#include "spin_constrain.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "basic_funcs.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop(int outer_step, bool rerun)
{
    int nat = this->get_nat();
    int ntype = this->get_ntype();
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> delta_lambda(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> dnu(nat, 0.0), dnu_last_step(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin(nat, 0.0), delta_spin(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> search(nat, 0.0), search_old(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> new_spin(nat, 0.0), spin_plus(nat, 0.0);

    double alpha_opt, alpha_plus;
    double beta;
    double mean_error, mean_error_old, rms_error;
    double g;

    double alpha_trial = this->alpha_trial_;

    const double zero = 0.0;
    const double one = 1.0;

#ifdef __MPI
	auto iterstart = MPI_Wtime();
#else
	auto iterstart = std::chrono::system_clock::now();
#endif

    double inner_loop_duration = 0.0;

    this->print_header();
    for (int i_step = -1; i_step < this->nsc_; i_step++)
    {
        double duration = 0.0;
        if (i_step == -1)
        {
            this->cal_mw_from_lambda(i_step);
            spin = this->Mi_;
            where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);
            print_2d("initial lambda (eV/uB): ", initial_lambda, this->nspin_, ModuleBase::Ry_to_eV);
            print_2d("initial spin (uB): ", spin, this->nspin_);
            print_2d("target spin (uB): ", this->target_mag_, this->nspin_);
            i_step++;
        }
        else
        {
            where_fill_scalar_else_2d(this->constrain_, 0, zero, delta_lambda, delta_lambda);
            add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);
        
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
            this->cal_mw_from_lambda(i_step, delta_lambda.data()); 
            new_spin = this->Mi_;
            bool GradLessThanBound = this->check_gradient_decay(new_spin, spin, delta_lambda, dnu_last_step);
            if (i_step >= this->nsc_min_ && GradLessThanBound)
            {
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
        subtract_2d(spin, this->target_mag_, delta_spin);
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
        search = delta_spin;
        if(this->direction_only_)
        for (int ia = 0; ia < nat; ia++)
        {
            const auto& target = this->target_mag_[ia];
            const double norm = std::sqrt(target.x*target.x + target.y*target.y + target.z*target.z);
            
            if (norm > 1e-8) {
                const ModuleBase::Vector3<double> dir = target / norm;
                const double parallel = delta_spin[ia].x*dir.x + delta_spin[ia].y*dir.y + delta_spin[ia].z*dir.z;
                temp_1[ia][0] = std::pow(delta_spin[ia].x,2) + std::pow(delta_spin[ia].y,2) + 
                                std::pow(delta_spin[ia].z,2) - std::pow(parallel,2);
                temp_1[ia][1] = 0;
                temp_1[ia][2] = 0;
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
        if(i_step == 0)
        {
            this->current_sc_thr_ = std::max(rms_error * this->sc_drop_thr_, this->sc_thr_);
        }
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
            this->update_psi_charge(dnu_last_step.data(), rerun, true);
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
                if(rms_error > this->current_sc_thr_ * 10 && rerun == true && this->higher_mag_prec == true)
                {
                    std::cout<<"Error: RMS error is too large, rerun the loop"<<std::endl;
                    this->run_lambda_loop(outer_step, false);
                }
            }
            break;
        }
#ifdef __MPI
		iterstart = MPI_Wtime();
#else
		iterstart = std::chrono::system_clock::now();
#endif
        if (i_step >= 2)
        {
            beta = mean_error / mean_error_old;
            add_scalar_multiply_2d(search, search_old, beta, search);
        }
        this->check_restriction(search, alpha_trial);

        dnu_last_step = dnu;
        add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);
        
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

        where_fill_scalar_else_2d(this->constrain_, 0, zero, delta_lambda, delta_lambda);
        add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);
        this->cal_mw_from_lambda(i_step, delta_lambda.data());

        spin_plus = this->Mi_;

        alpha_opt = this->cal_alpha_opt(spin, spin_plus, alpha_trial);
        this->check_restriction(search, alpha_opt);

        alpha_plus = alpha_opt - alpha_trial;
        scalar_multiply_2d(search, alpha_plus, temp_1);
        add_scalar_multiply_2d(dnu, temp_1, one, dnu);
        
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

        search_old = search;
        mean_error_old = mean_error;

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

    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    where_fill_scalar_else_2d(this->constrain_, 0, 0.0, this->lambda_, initial_lambda);

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

    ofs_scan << "# step  lambda_eV_uB";
    for (int ia = 0; ia < nat; ia++) {
        ofs_scan << "  Mi_x_" << ia << "  Mi_y_" << ia << "  Mi_z_" << ia;
    }
    ofs_scan << std::endl;

    double original_sc_thr = this->sc_thr_;

    for (int istep = 0; istep < nsteps; istep++) {
        double lambda_val_ry = lambda_start_ry + istep * lambda_step;
        double lambda_val_ev = lambda_val_ry * ModuleBase::Ry_to_eV;

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

        this->cal_mw_from_lambda(istep);

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

    ofs_scan.close();

    this->lambda_ = initial_lambda;

    std::cout << std::string(80, '=') << std::endl;
    std::cout << "[DS-DIAG] === LINEAR LAMBDA SCAN COMPLETE ===" << std::endl;
    std::cout << "[DS-DIAG] Results written to: lambda_scan_results.dat" << std::endl;
    std::cout << std::string(80, '=') << "\n" << std::endl;

    return;
}

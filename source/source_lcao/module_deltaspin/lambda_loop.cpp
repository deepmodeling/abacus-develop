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
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#endif

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
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_bfgs_v2(int outer_step)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "run_lambda_bfgs_v2");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "run_lambda_bfgs_v2");

    const int nat = this->get_nat();
    const double zero = 0.0;
    const double one = 1.0;

#ifdef __MPI
    auto iterstart = MPI_Wtime();
#else
    auto iterstart = std::chrono::system_clock::now();
#endif
    double inner_loop_duration = 0.0;

    std::cout << "[DS-DIAG] BFGS v2: nspin_ = " << this->nspin_ << ", npol_ = " << this->npol_ << ", nat = " << nat << std::endl;
    std::cout << "[DS-DIAG] BFGS v2: PARAM.inp.nspin = " << PARAM.inp.nspin << std::endl;
    std::cout << "[DS-DIAG] BFGS v2: constrain[0] = (" << this->constrain_[0].x << "," << this->constrain_[0].y << "," << this->constrain_[0].z << ")" << std::endl;
    if (nat > 1) std::cout << "[DS-DIAG] BFGS v2: constrain[1] = (" << this->constrain_[1].x << "," << this->constrain_[1].y << "," << this->constrain_[1].z << ")" << std::endl;
    struct DofInfo { int iat; int ic; };
    std::vector<DofInfo> dof_map;
    for (int ia = 0; ia < nat; ++ia) {
        for (int ic = 0; ic < 3; ++ic) {
            if (this->constrain_[ia][ic] != 0) {
                if (this->nspin_ == 2 && ic != 2) continue;
                dof_map.push_back({ia, ic});
            }
        }
    }
    std::cout << "[DS-DIAG] BFGS v2: dof_map size after initial = " << dof_map.size() << std::endl;
    for (const auto& d : dof_map) {
        std::cout << "[DS-DIAG]   DOF: atom=" << d.iat << " ic=" << d.ic << std::endl;
    }

    int n_dof = dof_map.size();
    if (n_dof == 0) {
        std::cout << "[DS-DIAG] BFGS v2: no constraints found in STRU, auto-setting all atoms as constrained" << std::endl;
        for (int ia = 0; ia < nat; ++ia) {
            if (this->nspin_ == 4) {
                this->constrain_[ia] = ModuleBase::Vector3<int>(1, 1, 1);
            } else {
                this->constrain_[ia] = ModuleBase::Vector3<int>(0, 0, 1);
            }
        }
        if (this->p_operator != nullptr) {
            if (this->nspin_ == 4) {
                auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator);
                if (dspin) dspin->reset_initialized();
            } else {
                auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator);
                if (dspin) dspin->reset_initialized();
            }
        }
        dof_map.clear();
        for (int ia = 0; ia < nat; ++ia) {
            for (int ic = 0; ic < 3; ++ic) {
                if (this->constrain_[ia][ic] != 0) {
                    if (this->nspin_ == 2 && ic != 2) continue;
                    dof_map.push_back({ia, ic});
                }
            }
        }
        n_dof = dof_map.size();
        if (n_dof == 0) {
            std::cout << "[DS-DIAG] BFGS v2: still no constrained DOFs after auto-setting, skipping" << std::endl;
            ModuleBase::timer::end("spinconstrain::SpinConstrain", "run_lambda_bfgs_v2");
            return;
        }
    }

    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);

    auto to_flat = [&](const std::vector<ModuleBase::Vector3<double>>& v) {
        std::vector<double> f(n_dof);
        for (int i = 0; i < n_dof; ++i)
            f[i] = v[dof_map[i].iat][dof_map[i].ic];
        return f;
    };

    auto apply_dl = [&](const std::vector<double>& dl_flat) {
        this->lambda_ = initial_lambda;
        for (int i = 0; i < n_dof; ++i)
            this->lambda_[dof_map[i].iat][dof_map[i].ic] = initial_lambda[dof_map[i].iat][dof_map[i].ic] + dl_flat[i];
    };

    this->cal_mw_from_lambda(-1);

    std::vector<ModuleBase::Vector3<double>> delta_spin(nat, 0.0);
    subtract_2d(this->Mi_, this->target_mag_, delta_spin);
    where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);

    std::vector<double> r = to_flat(delta_spin);
    double F = 0.0;
    for (int i = 0; i < n_dof; ++i) F += r[i] * r[i];
    F *= 0.5;

    double rms = std::sqrt(F * 2.0 / n_dof);

    if (outer_step == 0) {
        this->current_sc_thr_ = std::max(rms * this->sc_drop_thr_, this->sc_thr_);
    }

    std::vector<double> dl_flat(n_dof);
    auto dl_init = to_flat(initial_lambda);
    auto lambda_init_flat = to_flat(this->lambda_);
    for (int i = 0; i < n_dof; ++i) dl_flat[i] = lambda_init_flat[i] - dl_init[i];

    double alpha_init = this->alpha_trial_;
    std::vector<double> H_inv(n_dof * n_dof, 0.0);
    for (int i = 0; i < n_dof; ++i) H_inv[i * n_dof + i] = alpha_init;

    double restrict_cur = this->restrict_current_;

    this->print_header();
    std::cout << "[DS-DIAG] n_dof = " << n_dof << ", nsc = " << this->nsc_ << std::endl;
    std::cout << "[DS-DIAG] alpha_init = " << alpha_init * ModuleBase::Ry_to_eV << " eV/uB^2" << std::endl;
    std::cout << "[DS-DIAG] restrict_current = " << restrict_cur * ModuleBase::Ry_to_eV << " eV/uB" << std::endl;
    std::cout << "[DS-DIAG] current_sc_thr = " << this->current_sc_thr_ << std::endl;
    print_2d("initial lambda (eV/uB): ", initial_lambda, this->nspin_, ModuleBase::Ry_to_eV);
    print_2d("initial spin (uB): ", this->Mi_, this->nspin_);
    print_2d("target spin (uB): ", this->target_mag_, this->nspin_);
    std::cout << "[DS-DIAG] Initial F = " << F << ", RMS = " << rms << std::endl;

    std::ofstream ofs_bfgs;
    if (outer_step == 0) {
        ofs_bfgs.open("bfgs2_convergence.dat");
        ofs_bfgs << "# BFGS v2 Lambda Optimization Convergence" << std::endl;
        ofs_bfgs << "# n_dof = " << n_dof << ", alpha_init = " << alpha_init * ModuleBase::Ry_to_eV << " eV/uB^2" << std::endl;
    } else {
        ofs_bfgs.open("bfgs2_convergence.dat", std::ios::app);
    }
    ofs_bfgs << "#" << std::endl;
    ofs_bfgs << "# SCF outer step: " << outer_step << std::endl;
    ofs_bfgs << "# step  F  rms  |p_max|  alpha_eff  sy  reset" << std::endl;

    const int max_iter = this->nsc_;
    int reset_count = 0;
    bool converged = false;
    double g_norm = 0.0;
    for (int i = 0; i < n_dof; ++i) g_norm += r[i] * r[i];
    g_norm = std::sqrt(g_norm);

    for (int k = 0; k < max_iter; ++k)
    {
        std::vector<double> p(n_dof, 0.0);
        for (int i = 0; i < n_dof; ++i)
            for (int j = 0; j < n_dof; ++j)
                p[i] += H_inv[i * n_dof + j] * r[j];

        double rp = 0.0;
        for (int i = 0; i < n_dof; ++i) rp += r[i] * p[i];

        if (rp <= 0.0) {
            std::cout << "[DS-DIAG] BFGS v2: non-descent direction at step " << k
                      << " (rp=" << rp << "), resetting H" << std::endl;
            std::fill(H_inv.begin(), H_inv.end(), 0.0);
            for (int i = 0; i < n_dof; ++i) H_inv[i * n_dof + i] = alpha_init;
            for (int i = 0; i < n_dof; ++i) p[i] = alpha_init * r[i];
            rp = 0.0;
            for (int i = 0; i < n_dof; ++i) rp += r[i] * p[i];
            reset_count++;
        }

        double p_max = 0.0;
        for (int i = 0; i < n_dof; ++i) p_max = std::max(p_max, std::abs(p[i]));

        double alpha_eff = 1.0;
        if (p_max > restrict_cur) {
            alpha_eff = restrict_cur / p_max;
        }

        std::vector<double> dl_new(n_dof);
        for (int i = 0; i < n_dof; ++i) dl_new[i] = dl_flat[i] + alpha_eff * p[i];

        apply_dl(dl_new);

        std::vector<ModuleBase::Vector3<double>> dl_full(nat, 0.0);
        for (int i = 0; i < n_dof; ++i)
            dl_full[dof_map[i].iat][dof_map[i].ic] = dl_new[i];

        this->cal_mw_from_lambda(k, dl_full.data());
        std::vector<ModuleBase::Vector3<double>> Mi_new = this->Mi_;

        subtract_2d(Mi_new, this->target_mag_, delta_spin);
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
        std::vector<double> r_new = to_flat(delta_spin);

        std::vector<double> s(n_dof);
        for (int i = 0; i < n_dof; ++i) s[i] = alpha_eff * p[i];

        std::vector<double> y(n_dof);
        for (int i = 0; i < n_dof; ++i) y[i] = r_new[i] - r[i];

        double sy = 0.0;
        for (int i = 0; i < n_dof; ++i) sy += s[i] * y[i];

        if (sy > 1e-14) {
            double rho = 1.0 / sy;
            std::vector<double> Hy(n_dof, 0.0);
            for (int i = 0; i < n_dof; ++i)
                for (int j = 0; j < n_dof; ++j)
                    Hy[i] += H_inv[i * n_dof + j] * y[j];

            double yHy = 0.0;
            for (int i = 0; i < n_dof; ++i) yHy += y[i] * Hy[i];

            std::vector<double> H_new(n_dof * n_dof);
            for (int i = 0; i < n_dof; ++i)
                for (int j = 0; j < n_dof; ++j)
                    H_new[i * n_dof + j] = H_inv[i * n_dof + j]
                        - rho * (s[i] * Hy[j] + Hy[i] * s[j])
                        + rho * (1.0 + rho * yHy) * s[i] * s[j];

            H_inv = H_new;
        } else {
            std::cout << "[DS-DIAG] BFGS v2: curvature violated at step "
                      << k << " (sy=" << sy << "), skipping H update" << std::endl;
        }

        dl_flat = dl_new;
        r = r_new;
        F = 0.0;
        for (int i = 0; i < n_dof; ++i) F += r[i] * r[i];
        F *= 0.5;

        g_norm = 0.0;
        for (int i = 0; i < n_dof; ++i) g_norm += r[i] * r[i];
        g_norm = std::sqrt(g_norm);

        rms = std::sqrt(F * 2.0 / n_dof);
        for (int i = 0; i < n_dof; ++i) H_inv[i * n_dof + i] = alpha_init;

#ifdef __MPI
        double duration = (double)(MPI_Wtime() - iterstart);
#else
        double duration =
            (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now()
            - iterstart)).count() / static_cast<double>(1e6);
#endif
        inner_loop_duration += duration;

        std::cout << "[DS-DIAG] BFGS v2 step " << k
                  << " F=" << F << " rms=" << rms
                  << " |g|=" << g_norm << " |p_max|=" << p_max * ModuleBase::Ry_to_eV
                  << " alpha_eff=" << alpha_eff << " alpha_init=" << alpha_init * ModuleBase::Ry_to_eV << std::endl;

        ofs_bfgs << std::scientific << std::setprecision(6);
        ofs_bfgs << k << "  " << F << "  " << rms << "  " << p_max * ModuleBase::Ry_to_eV
                 << "  " << alpha_eff
                 << "  " << sy
                 << "  " << reset_count << std::endl;

        if (rms < this->current_sc_thr_) {
            std::cout << "[DS-DIAG] BFGS v2 converged at step " << k
                      << " (rms=" << rms << " < thr=" << this->current_sc_thr_ << ")" << std::endl;
            converged = true;
            break;
        }

#ifdef __MPI
        iterstart = MPI_Wtime();
#else
        iterstart = std::chrono::system_clock::now();
#endif
    }

    apply_dl(dl_flat);
    std::vector<ModuleBase::Vector3<double>> dl_final(nat, 0.0);
    for (int i = 0; i < n_dof; ++i)
        dl_final[dof_map[i].iat][dof_map[i].ic] = dl_flat[i];

    if (converged) {
        std::cout << "Meet convergence criterion ( < " << this->current_sc_thr_ << " ), exit.";
    } else {
        std::cout << "Reach maximum number of steps ( " << this->nsc_ << " ), exit.";
    }
    std::cout << "       Total TIME(s) = " << inner_loop_duration << std::endl;

    this->update_psi_charge(dl_final.data(), true, true);

    if (PARAM.inp.basis_type == "pw") {
        this->cal_mi_pw();
        subtract_2d(this->Mi_, this->target_mag_, delta_spin);
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
        double final_rms = 0.0;
        for (int ia = 0; ia < nat; ++ia)
            for (int ic = 0; ic < 3; ++ic)
                final_rms += delta_spin[ia][ic] * delta_spin[ia][ic];
        final_rms = std::sqrt(final_rms / nat);
        std::cout << "[DS-DIAG] BFGS v2 final RMS (after update_psi_charge): " << final_rms << std::endl;
    }

    ofs_bfgs.close();
    this->print_termination();

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "run_lambda_bfgs_v2");
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
    std::cout << "[DS-DIAG] dm_ = " << (this->dm_ ? "valid" : "NULL") << std::endl;
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
        if (this->p_operator != nullptr) {
            if (this->nspin_ == 4) {
                auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator);
                if (dspin) {
                    dspin->reset_initialized();
                }
            } else {
                auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator);
                if (dspin) {
                    dspin->reset_initialized();
                }
            }
        }
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

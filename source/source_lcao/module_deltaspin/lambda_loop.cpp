#include "spin_constrain.h"

#include <iostream>
#include <cmath>
#include <chrono>

#include "basic_funcs.h"
#include "source_io/module_parameter/parameter.h"

// lambda = initial_lambda + delta_lambda/(spin2 - spin1) * (target_spin - spin1)
/*inline void next_lambda(std::vector<ModuleBase::Vector3<double>>& initial_lambda,
                        std::vector<ModuleBase::Vector3<double>>& delta_lambda,
                        std::vector<ModuleBase::Vector3<double>>& lambda,
                        std::vector<ModuleBase::Vector3<double>>& spin1,
                        std::vector<ModuleBase::Vector3<double>>& spin2,
                        std::vector<ModuleBase::Vector3<double>>& target_spin)
{
    for (int ia = 0; ia < lambda.size(); ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            lambda[ia][ic] = initial_lambda[ia][ic] + delta_lambda[ia][ic] / (spin2[ia][ic] - spin1[ia][ic]) * (target_spin[ia][ic] - spin1[ia][ic]);
        }
    }
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop(int outer_step)
{
    // init parameters
    int nat = this->get_nat();
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> delta_lambda(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin2(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> delta_spin(nat, 0.0);
    // current lambda is this->lambda_
    // current spin is this->Mi_
    // target spin is this->target_mag_
    // loop to optimize lambda to get target spin
    int step = -1;
    do
    {
        // set initial lambda
        where_fill_scalar_else_2d(this->constrain_, 0, 0.0, this->lambda_, initial_lambda);
        // save current spin to spin1 if step > 0
        if (step > 0)
        {
            spin1 = this->Mi_;
        }
        // calculate current spin
        this->cal_mw_from_lambda(step);
        // save current spin to spin2
        spin2 = this->Mi_;
        // calculate delta_spin = target_spin - spin
        subtract_2d(this->target_mag_, spin2, delta_spin);
        // check RMS error and stop if needed
        // calculate RMS error
        double sum = 0.0;
        for (int ia = 0; ia < nat; ia++)
        {
            for (int ic = 0; ic < 3; ic++)
            {
                sum += std::pow(delta_spin[ia][ic],2);
            }
        }
        double rms_error = std::sqrt(sum/nat);
        std::cout << "RMS error = " << rms_error <<" in step:" <<step << std::endl;
        // check RMS error and stop if needed
        if(rms_error < 1e-5)
        {
            std::cout<<"success"<<std::endl;
            break;
        }
        // calculate delta_lambda
        if(1)//step == 0)
        {
            for(int ia = 0; ia < nat; ia++)
            {
                for(int ic = 2; ic < 3; ic++)
                {
                    delta_lambda[ia][ic] = 0.01;//- delta_spin[ia][ic] / 10.0;
                    this->lambda_[ia][ic] = initial_lambda[ia][ic] + delta_lambda[ia][ic];
                    std::cout<<__LINE__<<"lambda["<<ia<<"] = "<<this->lambda_[ia][ic]<<std::endl;
                }
            }
        }
        else
        {
            //calculate next lambda
            next_lambda(initial_lambda, delta_lambda, this->lambda_, spin1, spin2, this->target_mag_);
            // calculate delta_lambda = this->lambda - initial_lambda
            subtract_2d(this->lambda_, initial_lambda, delta_lambda);
        }
        step++;
    } while (step < this->nsc_);
    
}*/


template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop(
        int outer_step,
		bool rerun)
{
    // init controlling parameters
    int nat = this->get_nat();
    int ntype = this->get_ntype();
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> delta_lambda(nat,0.0);
    // set nu, dnu and dnu_last_step
    std::vector<ModuleBase::Vector3<double>> dnu(nat, 0.0), dnu_last_step(nat, 0.0);
    // two controlling temp variables
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin(nat, 0.0), delta_spin(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> search(nat, 0.0), search_old(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> new_spin(nat, 0.0), spin_plus(nat, 0.0);

    double alpha_opt, alpha_plus;
    double beta = 0.0, g = 0.0, mean_error = 0.0, mean_error_old = 0.0, rms_error = 0.0;

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
    // lambda loop
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
        
            // set the lambda component along the target magnetic moment direction to zero
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
                this->update_psi_charge(dnu_last_step.data());
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
        // continue the lambda loop
        subtract_2d(spin, this->target_mag_, delta_spin);
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
        search = delta_spin;
        // calculate the residual perpendicular to the target magnetic moment direction
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
            // set current_sc_thr_ to max(rms_error * sc_drop_thr, this->sc_thr_)
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
            //add_scalar_multiply_2d(initial_lambda, dnu_last_step, 1.0, this->lambda_);
            this->update_psi_charge(dnu_last_step.data(), rerun);
            if(PARAM.inp.basis_type == "pw")
            {
                //double check Atomic spin moment
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
        /// check if restriction is needed
        this->check_restriction(search, alpha_trial);

        dnu_last_step = dnu;
        add_scalar_multiply_2d(dnu, search, alpha_trial, dnu);
        
        // project delta_lambda to the target direction to ensure the increment update also meets the constraints
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

        // Cap delta_lambda to prevent explosion
        for(int ia=0; ia<nat; ++ia) {
            for(int ic=0; ic<3; ++ic) {
                if(std::abs(delta_lambda[ia][ic]) > 10.0) {
                    delta_lambda[ia][ic] = 10.0 * (delta_lambda[ia][ic] > 0 ? 1.0 : -1.0);
                }
            }
        }

        where_fill_scalar_else_2d(this->constrain_, 0, zero, delta_lambda, delta_lambda);
        add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);

        this->cal_mw_from_lambda(i_step, delta_lambda.data());

        spin_plus = this->Mi_;

        alpha_opt = this->cal_alpha_opt(spin, spin_plus, alpha_trial);
        /// check if restriction is needed
        this->check_restriction(search, alpha_opt);

        alpha_plus = alpha_opt - alpha_trial;
        scalar_multiply_2d(search, alpha_plus, temp_1);
        add_scalar_multiply_2d(dnu, temp_1, one, dnu);
        
        // project delta_lambda to ensure the increment update also meets the constraints
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

#ifdef __LCAO
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/elecstate_tools.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop_lcao(int outer_step)
{
    const int nat = this->get_nat();
    const int nks = this->kv_.get_nks();   // total k-points (spin-up + spin-down for nspin=2)
    const int nk = nks / 2;                // k-points per spin channel
    psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
    const int nbands = psi_t->get_nbands();
    const double alpha_damp = 0.8;
    const int max_inner_iter = 2;

    this->print_header();

    // ── Phase 1: Full diagonalization to get C_k, e_k, Mi ──
    this->cal_mw_from_lambda(-1);
    std::vector<ModuleBase::Vector3<double>> spin(nat);
    spin = this->Mi_;

    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat, 0.0);
    const double zero = 0.0;
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);

    print_2d("initial lambda (eV/uB): ", initial_lambda, this->nspin_, ModuleBase::Ry_to_eV);
    print_2d("initial spin (uB): ", spin, this->nspin_);
    print_2d("target spin (uB): ", this->target_mag_, this->nspin_);

    // Check initial convergence
    std::vector<ModuleBase::Vector3<double>> delta_spin(nat, 0.0);
    subtract_2d(spin, this->target_mag_, delta_spin);
    where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
    double rms_error = 0.0;
    {
        double sum = 0.0;
        for (int ia = 0; ia < nat; ia++)
            for (int ic = 0; ic < 3; ic++)
                sum += std::pow(delta_spin[ia][ic], 2);
        rms_error = std::sqrt(sum / nat);
    }
    this->current_sc_thr_ = std::max(rms_error * this->sc_drop_thr_, this->sc_thr_);

    if (rms_error < this->current_sc_thr_)
    {
        std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- 0"
                  << "       RMS = " << rms_error << std::endl;
        std::cout << "Meet convergence criterion ( < " << this->current_sc_thr_ << " ), exit." << std::endl;
        this->print_termination();
        // Update charge from current psi
        this->pelec->psiToRho(*psi_t);
        return;
    }

    // ── Phase 2: Compute P_I_sub for all k-points ──
    auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator);

    // PI_sub[ik][iat] = nbands × nbands Hermitian matrix
    std::vector<std::vector<std::vector<std::complex<double>>>> PI_sub(nks);
    for (int ik = 0; ik < nks; ik++)
    {
        psi_t->fix_k(ik);
        dspin_op->cal_PI_sub(this->kv_.kvec_d[ik], psi_t->get_pointer(), nbands, PI_sub[ik]);
    }

    // ── Phase 3: Analytical Jacobian ──
    // chi_I = dM_I^z / dlambda_I
    // For nspin=2: M_I = sum_k [sum_n f_n_up * P_I_nn_up - sum_n f_n_down * P_I_nn_down]
    // dM/dlambda uses perturbation theory with both spin channels
    std::vector<double> chi(nat, 0.0);
    for (int iat = 0; iat < nat; iat++)
    {
        if (this->constrain_[iat].z == 0) { continue;
        }
        double chi_val = 0.0;
        for (int ik = 0; ik < nks; ik++)
        {
            if (PI_sub[ik][iat].empty()) { continue;
            }
            // sign: +1 for spin-up (ik < nk), -1 for spin-down (ik >= nk)
            // dH_up/dlambda = +P_I, dH_down/dlambda = -P_I
            // dM/dlambda = d(M_up - M_down)/dlambda
            // For spin-up channel: dM_up/dlambda = sum_{n,m} 2*(f_n-f_m)*|P_nm|^2/(e_n-e_m) * (+1)
            // For spin-down channel: dM_down/dlambda = sum_{n,m} 2*(f_n-f_m)*|P_nm|^2/(e_n-e_m) * (-1)
            // dM/dlambda = dM_up/dlambda - dM_down/dlambda
            // Both channels contribute with same sign to chi
            const double sign = static_cast<double>(this->get_spin_sign(ik));
            const auto& P = PI_sub[ik][iat];
            for (int n = 0; n < nbands; n++)
            {
                const double fn = this->pelec->wg(ik, n);
                for (int m = n + 1; m < nbands; m++)
                {
                    const double fm = this->pelec->wg(ik, m);
                    const double de = this->pelec->ekb(ik, n) - this->pelec->ekb(ik, m);
                    if (std::abs(de) < 1e-10) { continue;
                    }
                    const double P_nm_sq = std::norm(P[n * nbands + m]);
                    // sign * sign = 1 always, so both channels add
                    chi_val += 2.0 * (fn - fm) * P_nm_sq / de;
                }
            }
        }
        chi[iat] = chi_val;
    }

    // ── Phase 4: Newton update + subspace verification ──
    // Storage for subspace diag results
    ModuleBase::matrix ekb_new(nks, nbands);
    ModuleBase::matrix wg_new(nks, nbands);
    std::vector<std::vector<std::complex<double>>> V_save(nks);

    for (int inner = 0; inner < max_inner_iter; inner++)
    {
        // Newton step: delta_lambda_I = alpha_damp * (target - current) / chi_I
        for (int iat = 0; iat < nat; iat++)
        {
            if (this->constrain_[iat].z == 0) { continue;
            }
            if (std::abs(chi[iat]) < 1e-15) { continue;
            }
            const double delta_lambda_z = alpha_damp * (this->target_mag_[iat].z - spin[iat].z) / chi[iat];
            this->lambda_[iat].z = initial_lambda[iat].z + delta_lambda_z;
        }

        // Subspace diag for each k-point
        for (int ik = 0; ik < nks; ik++)
        {
            const double sign = static_cast<double>(this->get_spin_sign(ik));

            // Build H_sub = diag(e_k) + sign * sum_I delta_lambda_I * P_I_sub(k)
            std::vector<std::complex<double>> H_sub(nbands * nbands, {0.0, 0.0});
            for (int n = 0; n < nbands; n++)
            {
                H_sub[n * nbands + n] = {this->pelec->ekb(ik, n), 0.0};
            }
            for (int iat = 0; iat < nat; iat++)
            {
                if (PI_sub[ik][iat].empty()) { continue;
                }
                const double dlambda = sign * (this->lambda_[iat].z - initial_lambda[iat].z);
                for (int i = 0; i < nbands * nbands; i++)
                {
                    H_sub[i] += dlambda * PI_sub[ik][iat][i];
                }
            }

            // Diag with LAPACK zheev
            std::vector<double> e_new(nbands);
            V_save[ik] = H_sub; // zheev overwrites with eigenvectors
            int lwork = 2 * nbands;
            std::vector<std::complex<double>> work(lwork);
            std::vector<double> rwork(3 * nbands);
            int info = 0;
            zheev_("V", "U", &nbands, V_save[ik].data(), &nbands,
                   e_new.data(), work.data(), &lwork, rwork.data(), &info);
            if (info != 0)
            {
                std::cout << "WARNING: zheev failed with info=" << info << " at ik=" << ik << std::endl;
            }
            for (int n = 0; n < nbands; n++)
            {
                ekb_new(ik, n) = e_new[n];
            }
        }

        // Recompute weights from new eigenvalues
        elecstate::calculate_weights(ekb_new,
                                     wg_new,
                                     this->pelec->klist,
                                     this->pelec->eferm,
                                     this->pelec->f_en,
                                     this->pelec->nelec_spin,
                                     this->pelec->skip_weights);

        // Compute Mi_new from subspace rotation
        std::vector<ModuleBase::Vector3<double>> Mi_new(nat, 0.0);
        for (int iat = 0; iat < nat; iat++)
        {
            if (this->constrain_[iat].z == 0) { continue;
            }
            double mi_z = 0.0;
            for (int ik = 0; ik < nks; ik++)
            {
                if (PI_sub[ik][iat].empty()) { continue;
                }
                const double sign = static_cast<double>(this->get_spin_sign(ik));
                const auto& V = V_save[ik];
                const auto& P = PI_sub[ik][iat];

                // P_rotated = V^dag P V, we only need diagonal elements
                // P_rotated[n,n] = sum_{a,b} conj(V[a,n]) * P[a,b] * V[b,n]
                for (int n = 0; n < nbands; n++)
                {
                    std::complex<double> pnn = {0.0, 0.0};
                    for (int a = 0; a < nbands; a++)
                    {
                        std::complex<double> tmp = {0.0, 0.0};
                        for (int b = 0; b < nbands; b++)
                        {
                            tmp += P[a * nbands + b] * V[b * nbands + n];
                        }
                        pnn += std::conj(V[a * nbands + n]) * tmp;
                    }
                    mi_z += sign * wg_new(ik, n) * pnn.real();
                }
            }
            Mi_new[iat].z = mi_z;
        }

        // Check convergence
        subtract_2d(Mi_new, this->target_mag_, delta_spin);
        where_fill_scalar_2d(this->constrain_, 0, zero, delta_spin);
        {
            double sum = 0.0;
            for (int ia = 0; ia < nat; ia++)
                for (int ic = 0; ic < 3; ic++)
                    sum += std::pow(delta_spin[ia][ic], 2);
            rms_error = std::sqrt(sum / nat);
        }

        std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << std::left << std::setw(5) << inner + 1
                  << "       RMS = " << rms_error << " (subspace)" << std::endl;

        if (rms_error < this->current_sc_thr_)
        {
            std::cout << "Meet convergence criterion ( < " << this->current_sc_thr_ << " ), exit." << std::endl;
            break;
        }

        // Update spin for next iteration
        spin = Mi_new;
    }

    this->print_termination();

    // ── Phase 5: Finalize — rotate wavefunctions and update DM/charge ──
    // C_new_k = C_k * V_k via pzgemm (2D-block distributed)
    // V_k is nbands × nbands (small, replicated on all procs)
    // C_k is nlocal × nbands (2D-block distributed)
    for (int ik = 0; ik < nks; ik++)
    {
        psi_t->fix_k(ik);
        const int nlocal = this->ParaV->get_row_size();
        const int ncol_local = this->ParaV->ncol_bands;

        // Temporary storage for rotated wavefunction
        std::vector<std::complex<double>> psi_new(nlocal * ncol_local, {0.0, 0.0});

        // C_new[irow, jcol_local] = sum_m C[irow, m_local] * V[m_global, jcol_global]
        // Since V is replicated, we can do this locally per process
        const std::complex<double>* psi_old = psi_t->get_pointer();
        for (int jcol_local = 0; jcol_local < ncol_local; jcol_local++)
        {
            const int jcol_global = this->ParaV->local2global_col(jcol_local);
            for (int mcol_local = 0; mcol_local < ncol_local; mcol_local++)
            {
                const int mcol_global = this->ParaV->local2global_col(mcol_local);
                // V[mcol_global, jcol_global] — V is column-major from zheev
                const std::complex<double> v_mj = V_save[ik][mcol_global * nbands + jcol_global];
                // psi_new[:, jcol_local] += psi_old[:, mcol_local] * v_mj
                for (int irow = 0; irow < nlocal; irow++)
                {
                    psi_new[irow + jcol_local * nlocal] += psi_old[irow + mcol_local * nlocal] * v_mj;
                }
            }
        }

        // Copy back
        std::complex<double>* psi_ptr = const_cast<std::complex<double>*>(psi_t->get_pointer());
        std::copy(psi_new.begin(), psi_new.end(), psi_ptr);

        // Update eigenvalues
        for (int n = 0; n < nbands; n++)
        {
            this->pelec->ekb(ik, n) = ekb_new(ik, n);
        }
    }

    // Update weights, DM, and charge
    elecstate::calculate_weights(this->pelec->ekb,
                                 this->pelec->wg,
                                 this->pelec->klist,
                                 this->pelec->eferm,
                                 this->pelec->f_en,
                                 this->pelec->nelec_spin,
                                 this->pelec->skip_weights);
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    elecstate::cal_dm_psi(this->ParaV, this->pelec->wg, *psi_t, *this->dm_);
    this->dm_->cal_DMR();
    this->pelec->psiToRho(*psi_t);
}
#endif // __LCAO

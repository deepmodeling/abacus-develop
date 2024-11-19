#include "spin_constrain.h"

#include <iostream>
#include <cmath>
#include <chrono>

#include "basic_funcs.h"
#include "module_parameter/parameter.h"

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
void spinconstrain::SpinConstrain<std::complex<double>>::run_lambda_loop(int outer_step, bool rerun)
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
            /*if(PARAM.inp.basis_type == "pw")
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
            }*/
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
        delta_lambda = dnu;

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

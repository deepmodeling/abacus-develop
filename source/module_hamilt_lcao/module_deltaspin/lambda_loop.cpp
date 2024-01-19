#include "spin_constrain.h"

#include <iostream>
#include <cmath>

#include "basic_funcs.h"

template<>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::run_lambda_loop(int outer_step)
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

    this->print_header();
    // lambda loop
    for (int i_step = 0; i_step < this->nsc_; i_step++)
    {
        if (i_step == 0)
        {
            spin = this->Mi_;
            where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);
            print_2d("initial lambda: ", initial_lambda, this->nspin_);
            print_2d("initial spin: ", spin, this->nspin_);
            print_2d("target spin: ", this->target_mag_, this->nspin_);
        }
        else
        {
            where_fill_scalar_else_2d(this->constrain_, 0, zero, delta_lambda, delta_lambda);
            add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);
            this->check_lock();
            this->cal_mw_from_lambda(i_step);
            new_spin = this->Mi_;
            bool GradLessThanBound = this->check_gradient_decay(new_spin, spin, delta_lambda, dnu_last_step);
            if (i_step >= this->nsc_min_ && GradLessThanBound)
            {
                add_scalar_multiply_2d(initial_lambda, dnu_last_step, one, this->lambda_);
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
        if (this->check_rms_stop(outer_step, i_step, rms_error))
        {
            add_scalar_multiply_2d(initial_lambda, dnu_last_step, 1.0, this->lambda_);
            break;
        }
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
        this->check_lock();
        this->cal_mw_from_lambda(i_step);

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
}

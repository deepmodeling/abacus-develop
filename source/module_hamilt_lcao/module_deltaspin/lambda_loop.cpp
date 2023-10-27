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
    std::vector<ModuleBase::Vector3<double>> dnu(nat,0.0), dnu_last_step(nat,0.0), nu_change(nat,1.0);
    // two controlling temp variables
    std::vector<ModuleBase::Vector3<double>> temp_1(nat,0.0), temp_2(nat,0.0);
    // MW during loop
    // new_spin replace MW in the original code
    // spin gradient
    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
        nat, std::vector<std::vector<std::vector<double>>>(
            3, std::vector<std::vector<double>>(
                nat, std::vector<double>(
                    3,0.0))));
    std::vector<ModuleBase::Vector3<double>> spin_nu_gradient_diag(nat,0.0);
    std::pair<int, int> maxloc(std::make_pair(0,0));
    std::vector<std::pair<int,int>> max_gradient_index(ntype, std::make_pair(0,0));
    std::vector<double> max_gradient(ntype,0.0);
    // question: bound_gradient should be read from INPUT?
    std::vector<double> bound_gradient(ntype,0.0);

    std::vector<ModuleBase::Vector3<double>> spin(nat,0.0), delta_spin(nat,0.0), delta_spin_old(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> search(nat,0.0), search_old(nat,0.0);

    std::vector<ModuleBase::Vector3<double>> spin_mask(nat,0.0), target_spin_mask(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> new_spin(nat,0.0), spin_change(nat,0.0), spin_plus(nat,0.0);
    std::vector<ModuleBase::Vector3<double>> spin_plus_mask(nat,0.0);

    double alpha_opt, alpha_plus;
    double beta;
    double mean_error, mean_error_old, rms_error;
    double boundary;
    double sum_k, sum_k2;
    double g;

    // calculate number of components to be constrained
    double num_component = sum_2d(this->constrain_);

    double alpha_trial = this->alpha_trial_;

    const double zero = 0.0;
    const double one = 1.0;

    for (int it = 0; it < ntype; it++)
    {
        bound_gradient[it] = this->decay_grad_[it];
        //std::cout << "bound_gradient " << bound_gradient[it] << std::endl;
    }

    this->print_header();
    // lambda loop
    for (int i_step = 0; i_step < this->nsc_; i_step++)
    {
        if (i_step == 0)
        {
            spin = this->Mi_;
            where_fill_scalar_else_2d(this->constrain_, 0, zero, this->lambda_, initial_lambda);
            print_2d("initial lambda: ", initial_lambda);
            print_2d("initial spin: ", spin);
            print_2d("target spin: ", this->target_mag_);
        }
        else
        {
            add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);
            this->cal_mw_from_lambda(i_step);
            new_spin = this->Mi_;
            subtract_2d(new_spin, spin, spin_change);
            subtract_2d(delta_lambda, dnu_last_step, nu_change);
            where_fill_scalar_2d(this->constrain_, 0, zero, spin_change);
            where_fill_scalar_2d(this->constrain_, 0, one, nu_change);
            // calculate spin_nu_gradient
            for (int ia = 0; ia < nat; ia++)
            {
                for (int ic = 0; ic < 3; ic++)
                {
                    for (int ja = 0; ja < nat; ja++)
                    {
                        for(int jc = 0; jc < 3; jc++)
                        {
                            spin_nu_gradient[ia][ic][ja][jc] = spin_change[ia][ic] / nu_change[ja][jc];
                        }
                    }
                }
            }
            for (const auto& sc_elem : this->get_atomCounts())
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
            //print_2d("diagonal gradient: ", spin_nu_gradient_diag);
            //std::cout << "maximum gradient appears at: " << std::endl;
            //for (int it = 0; it < ntype; it++)
            //{
            //    std::cout << "( " << max_gradient_index[it].first << ", " << max_gradient_index[it].second << " )" << std::endl;
            //}
            //std::cout << "maximum gradient: " << std::endl;
            //for (int it = 0; it < ntype; it++)
            //{
            //    std::cout << max_gradient[it] << std::endl;
            //}
            for (int it = 0; it < ntype; it++)
            {
                if (i_step >= this->nsc_min_ && bound_gradient[it] > 0 && std::abs(max_gradient[it]) < bound_gradient[it])
                {
                    std::cout << "Reach limitation of current step ( maximum gradient < " 
                        << bound_gradient[it] << " in atom type " << it << " ), exit." << std::endl;
                    // roll back to the last step
                    add_scalar_multiply_2d(initial_lambda, dnu_last_step, one, this->lambda_);
                    this->print_termination();
                    break;
                }
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
        mean_error = sum_2d(temp_1) / num_component;
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

        add_scalar_multiply_2d(initial_lambda, delta_lambda, one, this->lambda_);
        this->cal_mw_from_lambda(i_step);

        spin_plus = this->Mi_;

        where_fill_scalar_else_2d(this->constrain_, 0, zero, this->target_mag_, target_spin_mask);
        where_fill_scalar_else_2d(this->constrain_, 0, zero, spin, spin_mask);
        where_fill_scalar_else_2d(this->constrain_, 0, zero, spin_plus, spin_plus_mask);

        for (int ia = 0; ia < nat; ia++)
        {
            for (int ic = 0; ic < 3; ic++)
            {
                temp_1[ia][ic] = (target_spin_mask[ia][ic] - spin_mask[ia][ic]) * (spin_plus_mask[ia][ic] - spin_mask[ia][ic]);
                temp_2[ia][ic] = std::pow(spin_mask[ia][ic] - spin_plus_mask[ia][ic], 2);
            }
        }
        sum_k = sum_2d(temp_1);
        sum_k2 = sum_2d(temp_2);
        alpha_opt = sum_k * alpha_trial / sum_k2;
        /// check if restriction is needed
        this->check_restriction(search, alpha_opt);

        alpha_plus = alpha_opt - alpha_trial;
        scalar_multiply_2d(search, alpha_plus, temp_1);
        add_scalar_multiply_2d(dnu, temp_1, one, dnu);
        delta_lambda = dnu;

        search_old = search;
        delta_spin_old = delta_spin;
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
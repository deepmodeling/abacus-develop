#include "spin_constrain.h"

#include <iostream>
#include <cmath>

#include "basic_funcs.h"

//void LambdaLoop::init_input_parameters()
//{
    /// todo
    /// init input parameters from reading INPUT file
    //this->spin = GlobalV::MW;
    //this->out_lambda = GlobalV::OUT_LAMBDA;
    // question: how to out_lambda?
    //this->target_spin = GlobalV::M_CONSTR;
    //this->constrain = GlobalV::CONSTRL;
    //this->alpha_trial = GlobalV::INISC;
    //this->sc_thr = GlobalV::SCDIFF;
    //this->bound_gradient = GlobalV::SCCONV_GRAD;
    //this->nsc_ = GlobalV::NSC;
    //this->nsc_min_ = GlobalV::NSCMIN;
    //this->restrict_current = GlobalV::SCCUT;
//}

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::run_lambda_loop(int outer_step)
{
    std::cout << "outer_step = " << outer_step << std::endl;
    std::cout << "sc_thr " << this->sc_thr_ << std::endl;
    std::cout << "nsc " << this->nsc_ << std::endl;
    std::cout << "nsc_min " << this->nsc_min_ << std::endl;
//    // init controlling parameters
//    int nat = this->get_nat();
//    int ntype = this->get_ntype();
//    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat,0.0), delta_lambda(nat,0.0);
//    // question: how to set initial_lambda?
//    // question: is delta_lambda initially zero?
//    // set nu, dnu and dnu_last_step
//    std::vector<ModuleBase::Vector3<double>> nu(nat,0.0), dnu(nat,0.0), dnu_last_step(nat,0.0), nu_change(nat,0.0);
//    // two controlling temp variables
//    std::vector<ModuleBase::Vector3<double>> temp_1(nat,0.0), temp_2(nat,0.0);
//    // MW during loop
//    // new_spin replace MW in the original code
//    // spin gradient
//    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
//        nat, std::vector<std::vector<std::vector<double>>>(
//            3, std::vector<std::vector<double>>(
//                nat, std::vector<double>(
//                    3,0.0))));
//    std::vector<ModuleBase::Vector3<double>> spin_nu_gradient_diag(nat,0.0);
//    std::pair<int, int> maxloc(std::make_pair(0,0));
//    std::vector<std::pair<int,int>> max_gradient_index(ntype, std::make_pair(0,0));
//    std::vector<double> max_gradient(ntype,-1.0);
//    // question: bound_gradient should be read from INPUT?
//    std::vector<double> bound_gradient(ntype,0.0);
//    // temp variables
//
//    // calculate number of components to be constrained
//    int num_component = sum_2d(this->constrain);
//    // delta spin
//    std::vector<ModuleBase::Vector3<double>> delta_spin(nat,0.0), delta_spin_old(nat,0.0);
//    std::vector<ModuleBase::Vector3<double>> search(nat,0.0), search_old(nat,0.0);
//
//    std::vector<ModuleBase::Vector3<double>> spin_mask(nat,0.0), target_spin_mask(nat,0.0);
//    std::vector<ModuleBase::Vector3<double>> new_spin(nat,0.0), spin_change(nat,0.0), spin_plus(nat,0.0);
//    std::vector<ModuleBase::Vector3<double>> spin_plus_mask(nat,0.0);
//
//    double alpha_trial, alpha_opt, alpha_plus;
//    double beta;
//    double mean_error, mean_error_old, rms_error;
//    double restrict_current;
//    double boundary;
//    double sum_k, sum_k2;
//    double g;
//
//    std::cout << "===============================================================================" << std::endl;
//    std::cout << "Inner optimization for lambda begins ..." << std::endl;
//    std::cout << "Covergence criterion for the iteration: " << this->sc_thr << std::endl;
//    // lambda loop
//    for (int i_step = 0; i_step < this->nsc_; i_step++)
//    {
//        if (i_step == 0)
//        {
//            nu = this->out_lambda;
//            where_fill_scalar_else_2d(this->constrain, 0, 0.0, this->out_lambda, initial_lambda);
//            print_2d("initial lambda: ", initial_lambda);
//            print_2d("initial spin: ", this->spin);
//            print_2d("target spin: ", this->sc_mag_);
//        }
//        else
//        {
//            std::cout << "optimal delta lambda: " << std::endl;
//            for (int i=0; i< nat; i++)
//            {
//                delta_lambda[i].print();
//            }
//            add_scalar_multiply_2d(initial_lambda, delta_lambda, 1.0, temp_1);
//            /**
//             * TODO, also in-place change density CHTOT and orbital W, const 3 regular 3
//             * basically, use CHTOTL_RESERVE and temp_1(LAMBDA) recalculate V, then calculate H, then
//             * diagonalize H (in the subspace spanned by W_RESERVE), then use obtained orbitals W 
//             * calculate density CHTOT and orbital mag MW.
//             * Note that using CHTOTL instead of CHTOT is to recreate the H that has W_RESERVE as 
//             * eigenvectors
//            */
//            /// calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, new_spin, CHTOT, W);
//            subtract_2d(new_spin, spin, spin_change);
//            subtract_2d(delta_lambda, dnu_last_step, nu_change);
//            where_fill_scalar_2d(constrain, 0, 0.0, spin_change);
//            where_fill_scalar_2d(constrain, 0, 0.0, nu_change);
//            // calculate spin_nu_gradient
//            for (int ia = 0; ia < nat; ia++)
//            {
//                for (int ic = 0; ic < 3; ic++)
//                {
//                    for (int ja = 0; ja < nat; ja++)
//                    {
//                        for(int jc = 0; jc < 3; jc++)
//                        {
//                            spin_nu_gradient[ia][ic][ja][jc] = spin_change[ia][ic] / nu_change[ja][jc];
//                        }
//                    }
//                }
//            }
//            for (const auto& sc_elem : this->get_atomCounts())
//            {
//                int it = sc_elem.first;
//                int nat_it = sc_elem.second;
//                for (int ia = 0; ia < nat_it; ia++)
//                {
//                    for (int ic = 0; ic < 3; ic++)
//                    {
//                        spin_nu_gradient_diag[ia][ic] = spin_nu_gradient[ia][ic][ia][ic];
//                        if (spin_nu_gradient_diag[ia][ic] > max_gradient[it])
//                        {
//                            max_gradient[it] = spin_nu_gradient_diag[ia][ic];
//                            max_gradient_index[it].first = ia;
//                            max_gradient_index[it].second = ic;
//                        }
//                    }
//                }
//            }
//            print_2d("diagonal gradient: ", spin_nu_gradient_diag);
//            std::cout << "maximum gradient appears at: " << std::endl;
//            for (int it = 0; it < ntype; it++)
//            {
//                std::cout << "( " << max_gradient_index[it].first << ", " << max_gradient_index[it].second << " )" << std::endl;
//            }
//            std::cout << "maximum gradient: " << std::endl;
//            for (int it = 0; it < ntype; it++)
//            {
//                std::cout << max_gradient[it] << std::endl;
//            }
//            for (int it = 0; it < ntype; it++)
//            {
//                if (i_step >= this->nsc_min_ && bound_gradient[it] > 0 && max_gradient[it] < bound_gradient[it])
//                {
//                    std::cout << "Reach limitation of current step ( maximum gradient < " 
//                        << bound_gradient[it] << " in atom type " << it << " ), exit." << std::endl;
//                    // roll back to the last step
//                    // TODO
//                    // CHTOT = CHTOT_last_step;
//                    add_scalar_multiply_2d(initial_lambda, dnu_last_step, 1.0, out_lambda);
//                    goto CG_STOP;
//                }
//            }
//            spin = new_spin;
//            print_2d("new spin: ", spin);
//            print_2d("target spin: ", this->sc_mag_);
//        }
//        // continue the lambda loop
//        subtract_2d(spin, this->sc_mag_, delta_spin);
//        where_fill_scalar_2d(constrain, 0, 0.0, delta_spin);
//        search = delta_spin;
//        for (int ia = 0; ia < nat; ia++)
//        {
//            for (int ic = 0; ic < 3; ic++)
//            {
//                temp_1[ia][ic] = std::pow(delta_spin[ia][ic],2);
//            }
//        }
//        mean_error = sum_2d(temp_1) / num_component;
//        rms_error = std::sqrt(mean_error);
//        std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << i_step + 1 << "       RMS =" << rms_error << std::endl;
//
//        if (rms_error < this->sc_thr || i_step == this->nsc_ - 1)
//        {
//            if (rms_error < this->sc_thr)
//            {
//                std::cout << "Meet convergence criterion ( < " << sc_thr << " ), exit." << std::endl;
//            }
//            else if (i_step == this->nsc_ - 1)
//            {
//                std::cout << "Reach maximum number of steps ( " << this->nsc_ << " ), exit." << std::endl;
//            }
//            add_scalar_multiply_2d(initial_lambda, delta_lambda, 1.0, out_lambda);
//            goto CG_STOP;
//        }
//
//        if(i_step>=1)
//        {
//            beta = mean_error / mean_error_old;
//            temp_1 = search;
//            add_scalar_multiply_2d(temp_1, search_old, beta, search);
//        }
//
//        boundary = abs(alpha_trial * maxval_abs_2d(search));
//        std::cout << "restriction of this step = " << restrict_current << std::endl;
//        std::cout << "alpha_trial before restrict = " << alpha_trial << std::endl;
//        std::cout << "boundary before = " << boundary << std::endl;
//        std::cout << "trial need restriction: false" << std::endl;
//        scalar_multiply_2d(search, alpha_trial, temp_1);
//        print_2d("delta delta lambda: ", temp_1);
//
//        // CHTOT_last_step = CHTOT;
//        dnu_last_step = dnu;
//        temp_1 = dnu;
//        add_scalar_multiply_2d(temp_1, search, alpha_trial, dnu);
//        delta_lambda = dnu;
//
//        print_2d("trial delta lambda:", delta_lambda);
//
//        if (debug)
//        {
//            print_2d("(Debug) before-trial-step spin:", spin);
//            print_2d("(Debug) target spin:", this->sc_mag_);
//        }
//
//        add_scalar_multiply_2d(initial_lambda, delta_lambda, 1.0, temp_1);
//        // TODO
//        // calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, new_spin, CHTOT, W);
//
//        spin_plus = new_spin;
//        print_2d("current spin(trial):", spin_plus);
//        //
//        where_fill_scalar_else_2d(constrain, 0, 0.0, this->sc_mag_, target_spin_mask);
//        where_fill_scalar_else_2d(constrain, 0, 0.0, spin, spin_mask);
//        where_fill_scalar_else_2d(constrain, 0, 0.0, spin_plus, spin_plus_mask);
//
//        for (int ia = 0; ia < nat; ia++)
//        {
//            for (int ic = 0; ic < 3; ic++)
//            {
//                temp_1[ia][ic] = (target_spin_mask[ia][ic] - spin_mask[ia][ic]) * (spin_plus_mask[ia][ic] - spin_mask[ia][ic]);
//                temp_2[ia][ic] = pow(spin_mask[ia][ic] - spin_plus_mask[ia][ic], 2);
//            }
//        }
//        sum_k = sum_2d(temp_1);
//        sum_k2 = sum_2d(temp_2);
//        alpha_opt = sum_k * alpha_trial / sum_k2;
//
//        boundary = abs(alpha_trial * maxval_abs_2d(search));
//        std::cout << "alpha_opt before restrict = " << alpha_opt << std::endl;
//        std::cout << "boundary before = " << boundary << std::endl;
//
//        if (restrict_current > 0 && boundary > restrict_current)
//        {
//            alpha_opt = copysign(1.0, alpha_opt) * restrict_current / maxval_abs_2d(search);
//            boundary = abs(alpha_opt * maxval_abs_2d(search));
//            std::cout << "restriction needed: true" << std::endl;
//            std::cout << "alpha_opt after restrict = " << alpha_opt << std::endl;
//            std::cout << "boundary after = " << boundary << std::endl;
//        }
//        else
//        {
//            std::cout << "restriction needed: false" << std::endl;
//        }
//
//        alpha_plus = alpha_opt - alpha_trial;
//        scalar_multiply_2d(search, alpha_plus, temp_1);
//        print_2d("delta delta lambda:", temp_1);
//
//        temp_2 = dnu;
//        add_scalar_multiply_2d(temp_2, temp_1, 1.0, dnu);
//        delta_lambda = dnu;
//
//        search_old = search;
//        delta_spin_old = delta_spin;
//        mean_error_old = mean_error;
//
//        g = 1.5 * abs(alpha_opt) / alpha_trial;
//        if (g > 2.0)
//        {
//            g = 2;
//        }
//        else if (g < 0.5)
//        {
//            g = 0.5;
//        }
//        else
//        {
//            std::cout << "Warning: g is not in the range of [0.5, 2.0], g = " << g << std::endl;
//        }
//        alpha_trial = alpha_trial * pow(g, 0.7);
//    }
//    
//CG_STOP:
//
//    // TODO
//    // CHTOTL = CHTOTL_RESERVE;
//    if (debug)
//    {
//        print_2d("(Debug) after-optimization spin: (print in the inner loop): ", new_spin);
//        print_2d("target spin: ", this->sc_mag_);
//    }
//    std::cout << "Inner optimization for lambda ends." << std::endl;
//    std::cout << "===============================================================================" << std::endl;
//
}

template class SpinConstrain<double, psi::DEVICE_CPU>;
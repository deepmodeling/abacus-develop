#include "lambda_loop.h"

#include <iostream>

#include "basic_funcs.h"

LambdaLoop::LambdaLoop()
{}

LambdaLoop::~LambdaLoop()
{}

void LambdaLoop::init_input_parameters()
{
    /// todo
    /// init input parameters from reading INPUT file
    //this->spin = GlobalV::MW;
    //this->out_lambda = GlobalV::OUT_LAMBDA;
    //this->target_spin = GlobalV::M_CONSTR;
    //this->constrain = GlobalV::CONSTRL;
    //this->alpha_trial = GlobalV::INISC;
    //this->epsilon = GlobalV::SCDIFF;
    //this->bound_gradient = GlobalV::SCCONV_GRAD;
    //this->num_step = GlobalV::NSC;
    //this->num_min = GlobalV::NSCMIN;
    //this->restrict_current = GlobalV::SCCUT;
}

void LambdaLoop::run_lambda_loop(int outer_step)
{
    // init controlling parameters
    std::vector<ModuleBase::Vector3<double>> initial_lambda(nat,0.0), delta_lambda(nat,0.0);
    // question: is delta_lambda initially zero?
    // set nu, dnu and dnu_last_step
    std::vector<ModuleBase::Vector3<double>> nu(nat,0.0), dnu(nat,0.0), dnu_last_step(nat,0.0), nu_change(nat,0.0);
    // two controlling temp variables
    std::vector<ModuleBase::Vector3<double>> temp_1(nat,0.0), temp_2(nat,0.0);
    // MW during loop
    std::vector<ModuleBase::Vector3<double>> new_spin(nat,0.0), spin_change(nat,0.0);
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
    std::vector<double> max_gradient(ntype,-1.0);
    std::vector<double> bound_gradient(ntype,0.0);
    // temp variables

    // calculate number of components to be constrained
    int num_component = sum_2d(this->constrain);

    std::cout << "===============================================================================" << std::endl;
    std::cout << "Inner optimization for lambda begins ..." << std::endl;
    std::cout << "Covergence criterion for the iteration: " << this->epsilon << std::endl;
    // lambda loop
    for (int i_step = 0; i_step < num_step; i_step++)
    {
        if (i_step == 0)
        {
            nu = this->out_lambda;
            where_fill_scalar_else_2d(constrain, 0, 0.0, out_lambda, initial_lambda);
            print_2d("initial lambda: ", initial_lambda);
            print_2d("initial spin: ", spin);
            print_2d("target spin: ", target_spin);
        }
        else
        {
            std::cout << "optimal delta lambda: " << std::endl;
            for (int i=0; i< this->nat; i++)
            {
                delta_lambda[i].print();
            }
            add_scalar_multiply_2d(initial_lambda, delta_lambda, 1.0, temp_1);
            /**
             * TODO, also in-place change density CHTOT and orbital W, const 3 regular 3
             * basically, use CHTOTL_RESERVE and temp_1(LAMBDA) recalculate V, then calculate H, then
             * diagonalize H (in the subspace spanned by W_RESERVE), then use obtained orbitals W 
             * calculate density CHTOT and orbital mag MW.
             * Note that using CHTOTL instead of CHTOT is to recreate the H that has W_RESERVE as 
             * eigenvectors
            */
            /// calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, new_spin, CHTOT, W);
            subtract_2d(new_spin, spin, spin_change);
            subtract_2d(delta_lambda, dnu_last_step, nu_change);
            where_fill_scalar_2d(constrain, 0, 0.0, spin_change);
            where_fill_scalar_2d(constrain, 0, 0.0, nu_change);
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
            for (int ia = 0; ia < nat; ia++)
            {
                int it = iat2it[ia];
                for (int ic = 0; ic < 3; ic++)
                {
                    spin_nu_gradient_diag[ia][ic] = spin_nu_gradient[ia][ic][ia][ic];
                    if (spin_nu_gradient_diag[ia][ic] > max_gradient[it])
                    {
                        max_gradient[it] = spin_nu_gradient_diag[ia][ic];
                        max_gradient_index[it].first = ia;
                        max_gradient_index[it].second = ic;
                    }
                }
            }
            print_2d("diagonal gradient: ", spin_nu_gradient_diag);
            //new_spin[ia] = spin_nu_gradient_diag[ia];
            //int it = iat2it[ia];
            //maxloc = maxloc_abs_2d(new_spin);
        }
    }
}
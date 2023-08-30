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
    std::vector<ModuleBase::Vector3<double>> initial_lambda, delta_lambda;
    initial_lambda.resize(nat);
    // question: is delta_lambda initially zero?
    delta_lambda.resize(nat);
    // set nu, dnu and dnu_last_step
    std::vector<ModuleBase::Vector3<double>> nu, dnu, dnu_last_step;
    nu.resize(nat);
    dnu.resize(nat);
    dnu_last_step.resize(nat);
    // two controlling temp variables
    std::vector<ModuleBase::Vector3<double>> temp_1, temp_2;
    temp_1.resize(nat);
    temp_2.resize(nat);

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
            std::cout << "initial lambda:" << std::endl;
            print_2d(initial_lambda);
            std::cout << "initial spin: " << std::endl;
            print_2d(spin);
            std::cout << "target spin: " << std::endl;
            print_2d(target_spin);
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
            /// calculate_MW_from_lambda(temp_1, CHTOTL_RESERVE, W_RESERVE, MW, CHTOT, W);
        }
    }
}
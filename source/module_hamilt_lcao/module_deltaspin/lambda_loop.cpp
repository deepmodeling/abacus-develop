#include "lambda_loop.h"

#include <iostream>

#include "basic_funcs.h"

LambdaLoop::LambdaLoop()
{}

LambdaLoop::~LambdaLoop()
{}

void LambdaLoop::init_input_parameters()
{
    // todo
    // init input parameters from reading INPUT file
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

int LambdaLoop::cal_num_component()
{
    // Vector3 is zero automatically
    ModuleBase::Vector3<int> sum_control;
    for (const auto& v : this->constrain) {
        sum_control += v;
    }
    int num_component = sum_control.x + sum_control.y + sum_control.z;
    return num_component;
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

    // calculate number of components to be constrained
    int num_component = this->cal_num_component();

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
            for (int i=0; i< this->nat; i++)
            {
                initial_lambda[i].print();
            }
            std::cout << "initial spin: " << std::endl;
            for (int i=0; i< this->nat; i++)
            {
                this->spin[i].print();
            }
            std::cout << "target spin: " << std::endl;
            for (int i=0; i< this->nat; i++)
            {
                this->target_spin[i].print();
            }
        }
        else
        {
            //nu = nu + dnu;
        }
    }
}
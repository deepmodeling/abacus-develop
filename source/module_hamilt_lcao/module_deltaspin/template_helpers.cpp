#include "spin_constrain.h"

template <>
void spinconstrain::SpinConstrain<double>::cal_mw_from_lambda(int i_step, const ModuleBase::Vector3<double>* delta_lambda)
{
}

template <>
void spinconstrain::SpinConstrain<double>::cal_mi_lcao(const int& step, bool print)
{
}

template <>
void spinconstrain::SpinConstrain<double>::run_lambda_loop(int outer_step, bool rerun)
{
}

template <>
bool spinconstrain::SpinConstrain<double>::check_rms_stop(int outer_step,
                                                                    int i_step,
                                                                    double rms_error,
                                                                    double duration,
                                                                    double total_duration)
{
    return false;
}

template <>
void spinconstrain::SpinConstrain<double>::check_restriction(
    const std::vector<ModuleBase::Vector3<double>>& search,
    double& alpha_trial)
{
}

/// calculate alpha_opt
template <>
double spinconstrain::SpinConstrain<double>::cal_alpha_opt(std::vector<ModuleBase::Vector3<double>> spin,
                                                                     std::vector<ModuleBase::Vector3<double>> spin_plus,
                                                                     const double alpha_trial)
{
    return 0.0;
}

template <>
void spinconstrain::SpinConstrain<double>::print_termination()
{
}

template <>
void spinconstrain::SpinConstrain<double>::print_header()
{
}

template <>
bool spinconstrain::SpinConstrain<double>::check_gradient_decay(
    std::vector<ModuleBase::Vector3<double>> new_spin,
    std::vector<ModuleBase::Vector3<double>> old_spin,
    std::vector<ModuleBase::Vector3<double>> new_delta_lambda,
    std::vector<ModuleBase::Vector3<double>> old_delta_lambda,
    bool print)
{
    return false;
}
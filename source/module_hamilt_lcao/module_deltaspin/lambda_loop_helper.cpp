#include "basic_funcs.h"
#include "spin_constrain.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::print_termination()
{
    print_2d("after-optimization spin: (print in the inner loop): ", this->Mi_);
    std::cout << "Inner optimization for lambda ends." << std::endl;
    std::cout << "===============================================================================" << std::endl;
}

template <>
bool SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_rms_stop(int outer_step, int i_step, double rms_error)
{
    std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << std::setw(5) << i_step + 1
              << "       RMS =" << rms_error << std::endl;
    if (rms_error < this->sc_thr_ || i_step == this->nsc_ - 1)
    {
        if (rms_error < this->sc_thr_)
        {
            std::cout << "Meet convergence criterion ( < " << this->sc_thr_ << " ), exit." << std::endl;
        }
        else if (i_step == this->nsc_ - 1)
        {
            std::cout << "Reach maximum number of steps ( " << this->nsc_ << " ), exit." << std::endl;
        }
        this->print_termination();
        return true;
    }
    return false;
}

/// print header
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::print_header()
{
    std::cout << "===============================================================================" << std::endl;
    std::cout << "Inner optimization for lambda begins ..." << std::endl;
    std::cout << "Covergence criterion for the iteration: " << this->sc_thr_ << std::endl;
}

/// check restriction
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_restriction(
    const std::vector<ModuleBase::Vector3<double>>& search,
    double& alpha_trial)
{
    double boundary = std::abs(alpha_trial * maxval_abs_2d(search));

    if (this->restrict_current_ > 0 && boundary > this->restrict_current_)
    {
        alpha_trial = copysign(1.0, alpha_trial) * this->restrict_current_ / maxval_abs_2d(search);
        boundary = std::abs(alpha_trial * maxval_abs_2d(search));
        std::cout << "alpha after restrict = " << alpha_trial * ModuleBase::Ry_to_eV << std::endl;
        std::cout << "boundary after = " << boundary * ModuleBase::Ry_to_eV << std::endl;
    }
}

/// calculate alpha_opt
template <>
double SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_alpha_opt(
    std::vector<ModuleBase::Vector3<double>> spin,
    std::vector<ModuleBase::Vector3<double>> spin_plus,
    const double alpha_trial)
{
    int nat = this->get_nat();
    const double zero = 0.0;
    std::vector<ModuleBase::Vector3<double>> spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> target_spin_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> spin_plus_mask(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_1(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> temp_2(nat, 0.0);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, this->target_mag_, target_spin_mask);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, spin, spin_mask);
    where_fill_scalar_else_2d(this->constrain_, 0, zero, spin_plus, spin_plus_mask);

    for (int ia = 0; ia < nat; ia++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            temp_1[ia][ic]
                = (target_spin_mask[ia][ic] - spin_mask[ia][ic]) * (spin_plus_mask[ia][ic] - spin_mask[ia][ic]);
            temp_2[ia][ic] = std::pow(spin_mask[ia][ic] - spin_plus_mask[ia][ic], 2);
        }
    }
    double sum_k = sum_2d(temp_1);
    double sum_k2 = sum_2d(temp_2);
    return sum_k * alpha_trial / sum_k2;
}

/// check gradient decay
template <>
bool SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_gradient_decay(
    std::vector<ModuleBase::Vector3<double>> new_spin,
    std::vector<ModuleBase::Vector3<double>> old_spin,
    std::vector<ModuleBase::Vector3<double>> new_delta_lambda,
    std::vector<ModuleBase::Vector3<double>> old_delta_lambda)
{
    return false;
}
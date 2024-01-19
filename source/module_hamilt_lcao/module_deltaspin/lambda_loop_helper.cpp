#include "basic_funcs.h"
#include "spin_constrain.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::print_termination()
{
    print_2d("after-optimization spin: (print in the inner loop): ", this->Mi_, this->nspin_);
    print_2d("after-optimization lambda: (print in the inner loop): ", this->lambda_, this->nspin_);
    std::cout << "Inner optimization for lambda ends." << std::endl;
    std::cout << "===============================================================================" << std::endl;
}

template <>
bool SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_rms_stop(int outer_step, int i_step, double rms_error)
{
    std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << std::left << std::setw(5) << i_step + 1
              << "       RMS = " << rms_error << std::endl;
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
    std::vector<ModuleBase::Vector3<double>> spin,
    std::vector<ModuleBase::Vector3<double>> delta_lambda,
    std::vector<ModuleBase::Vector3<double>> dnu_last_step,
    bool print)
{
    const double one = 1.0;
    const double zero = 0.0;
    int nat = this->get_nat();
    int ntype = this->get_ntype();
    std::vector<ModuleBase::Vector3<double>> spin_change(nat, 0.0);
    std::vector<ModuleBase::Vector3<double>> nu_change(nat, 1.0);
    std::vector<std::vector<std::vector<std::vector<double>>>> spin_nu_gradient(
        nat,
        std::vector<std::vector<std::vector<double>>>(
            3,
            std::vector<std::vector<double>>(nat, std::vector<double>(3, 0.0))));
    std::vector<ModuleBase::Vector3<double>> spin_nu_gradient_diag(nat, 0.0);
    std::vector<std::pair<int, int>> max_gradient_index(ntype, std::make_pair(0, 0));
    std::vector<double> max_gradient(ntype, 0.0);
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
                for (int jc = 0; jc < 3; jc++)
                {
                    spin_nu_gradient[ia][ic][ja][jc] = spin_change[ia][ic] / nu_change[ja][jc];
                }
            }
        }
    }
    for (const auto& sc_elem: this->get_atomCounts())
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
    if (print)
    {
        print_2d("diagonal gradient: ", spin_nu_gradient_diag, this->nspin_);
        std::cout << "maximum gradient appears at: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            std::cout << "( " << max_gradient_index[it].first << ", " << max_gradient_index[it].second << " )"
                      << std::endl;
        }
        std::cout << "maximum gradient: " << std::endl;
        for (int it = 0; it < ntype; it++)
        {
            std::cout << max_gradient[it] << std::endl;
        }
    }
    for (int it = 0; it < ntype; it++)
    {
        if (this->decay_grad_[it] > 0 && std::abs(max_gradient[it]) < this->decay_grad_[it])
        {
            std::cout << "Reach limitation of current step ( maximum gradient < " << this->decay_grad_[it]
                      << " in atom type " << it << " ), exit." << std::endl;
            return true;
        }
    }
    return false;
}

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_lock()
{
    int nat = this->get_nat();
    for (int iat = 0; iat < nat; iat++)
    {
        for (int jat = 0; jat < iat; jat++)
        {
            for (int idir = 0; idir <3; idir++)
            {
                if (this->constrain_[iat][idir] == 1 && this->constrain_[jat][idir] == 1)
                {
                    if (std::abs(this->target_mag_[iat][idir]) < 1e-9 &&
                        std::abs(this->target_mag_[jat][idir]) < 1e-9 &&
                        std::abs(this->lambda_[iat][idir]) > 1e-9 &&
                        std::abs(this->lambda_[jat][idir]) > 1e-9 &&
                        std::abs(this->lambda_[iat][idir] + this->lambda_[jat][idir]) < 1e-9)
                    {
                        this->lambda_[iat][idir] = 0.0;
                        this->lambda_[jat][idir] = 0.0;
                    }
                }
            }
        }
    }
}
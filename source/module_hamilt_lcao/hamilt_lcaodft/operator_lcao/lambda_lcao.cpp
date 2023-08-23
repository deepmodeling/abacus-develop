#include "lambda_lcao.h"

namespace hamilt
{

template class OperatorLambda<OperatorLCAO<double>>;

template class OperatorLambda<OperatorLCAO<std::complex<double>>>;

// destructor
template <>
OperatorLambda<OperatorLCAO<double>>::~OperatorLambda()
{
}

template <>
OperatorLambda<OperatorLCAO<std::complex<double>>>::~OperatorLambda()
{
}

// contribute to HR is not needed.
template <>
void OperatorLambda<OperatorLCAO<double>>::contributeHR()
{
}

// contribute to HR is not needed.
template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::contributeHR()
{
}

template <>
void OperatorLambda<OperatorLCAO<double>>::set_nat(int nat_in)
{
    this->nat_ = nat_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_nat(int nat_in)
{
    this->nat_ = nat_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::set_lambda(const std::vector<ModuleBase::Vector3<double>>& lambda_in)
{
    if (lambda_in.size() != this->nat_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_loc_lambda", "lambda_in size mismatch with nat");
    }
    this->loc_lambda = lambda_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_lambda(
    const std::vector<ModuleBase::Vector3<double>>& lambda_in)
{
    if (lambda_in.size() != this->nat_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_loc_lambda", "lambda_in size mismatch with nat");
    }
    this->loc_lambda = lambda_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::cal_h_lambda(int ik, double* h_lambda)
{
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::cal_h_lambda(int ik, std::complex<double>* h_lambda)
{
    ModuleBase::TITLE("OperatorLambda", "cal_h_lambda");
    ModuleBase::timer::tick("OperatorLambda", "cal_h_lambda");
    // Pauli matrix is here
    std::vector<std::vector<std::complex<double>>> sigma_x;
    std::vector<std::vector<std::complex<double>>> sigma_y;
    std::vector<std::vector<std::complex<double>>> sigma_z;
    sigma_x = {
        {std::complex<double>(0, 0), std::complex<double>(1, 0)},
        {std::complex<double>(1, 0), std::complex<double>(0, 0)}
    };

    sigma_y = {
        {std::complex<double>(0, 0), std::complex<double>(0, -1)},
        {std::complex<double>(0, 1), std::complex<double>(0, 0) }
    };

    sigma_z = {
        {std::complex<double>(1, 0), std::complex<double>(0,  0)},
        {std::complex<double>(0, 0), std::complex<double>(-1, 0)}
    };
    // lambda is transferred from outside in the constructor of this class
    // this->loc_lambda;
    // W_mu_nu is temporarily just the overlap matrix
    // this->LM->Sloc2;
    // h_lambda = W_i * (lambda_x *sigma_x + lambda_y * sigma_y + lambda_z * sigma_z)
    ModuleBase::timer::tick("OperatorLambda", "cal_h_lambda");
}

// contribute to Hk.
template <>
void OperatorLambda<OperatorLCAO<double>>::contributeHk(int ik)
{
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorLambda", "contributeHk");
    ModuleBase::timer::tick("OperatorLambda", "contributeHk");
    std::vector<std::complex<double>> h_lambda(this->LM->ParaV->nloc);
    this->cal_h_lambda(ik, &h_lambda[0]);

    for (int irc = 0; irc < this->LM->ParaV->nloc; irc++)
    {
        this->LM->Hloc2[irc] += h_lambda[irc];
    }
    ModuleBase::timer::tick("OperatorLambda", "contributeHk");
}

} // namespace hamilt
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
void OperatorLambda<OperatorLCAO<double>>::set_nloc(int nloc_in)
{
    this->nloc_ = nloc_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_nloc(int nloc_in)
{
    this->nloc_ = nloc_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::set_npol(int npol_in)
{
    this->npol_ = npol_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_npol(int npol_in)
{
    this->npol_ = npol_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::set_iwt2iat(const std::vector<int>& iwt2iat_in)
{
    if (iwt2iat_in.size() != this->nloc_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_iwt2iat", "iwt2iat_in size mismatch with nloc");
    }
    this->iwt2iat_ = iwt2iat_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_iwt2iat(const std::vector<int>& iwt2iat_in)
{
    if (iwt2iat_in.size() != this->nloc_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_iwt2iat", "iwt2iat_in size mismatch with nloc");
    }
    this->iwt2iat_ = iwt2iat_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::set_lambda(const std::vector<ModuleBase::Vector3<double>>& lambda_in)
{
    if (lambda_in.size() != this->nat_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_lambda", "lambda_in size mismatch with nat");
    }
    this->lambda_ = lambda_in;
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::set_lambda(
    const std::vector<ModuleBase::Vector3<double>>& lambda_in)
{
    if (lambda_in.size() != this->nat_)
    {
        ModuleBase::WARNING_QUIT("OperatorLambda::set_lambda", "lambda_in size mismatch with nat");
    }
    this->lambda_ = lambda_in;
}

template <>
void OperatorLambda<OperatorLCAO<double>>::cal_weight_func(const std::vector<double>& Sloc2)
{
    this->W_i_.reserve(this->nloc_ * this->npol_);
    for (int i = 0; i < nloc_ * npol_; i++)
    {
        int iat = iwt2iat_[i];
        for (int j = i; j < nloc_ * npol_; j++)
        {
            int jat = iwt2iat_[j];
            this->W_i_[i * nloc_ * npol_ + j]
                = (iat == jat) ? Sloc2[i * nloc_ * npol_ + j] : Sloc2[i * nloc_ * npol_ + j] * 0.5;
        }
    }
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::cal_weight_func(const std::vector<std::complex<double>>& Sloc2)
{
    this->W_i_.reserve(this->nloc_ * this->npol_);
    for (int i = 0; i < nloc_ * npol_; i++)
    {
        int iat = iwt2iat_[i];
        for (int j = i; j < nloc_ * npol_; j++)
        {
            int jat = iwt2iat_[j];
            this->W_i_[i * nloc_ * npol_ + j]
                = (iat == jat) ? Sloc2[i * nloc_ * npol_ + j] : Sloc2[i * nloc_ * npol_ + j] * 0.5;
        }
    }
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
    // h_lambda = W_i * (lambda_x *sigma_x + lambda_y * sigma_y + lambda_z * sigma_z)
    // Pauli matrix is used implicitly
    for (int i = 0; i < nloc_ * npol_; i++)
    {
        int iat = iwt2iat_[i];
        for (int j = i; j < nloc_ * npol_; j++)
        {
            int index = i * nloc_ * npol_ + j;
            if (i % 2 == 0)
            {
                if (j % 2 == 0)
                {
                    // H11
                    h_lambda[index] = W_i_[index] * lambda_[iat][2];
                }
                else
                {
                    // H12
                    h_lambda[index] = W_i_[index] * (lambda_[iat][0] + lambda_[iat][1] * std::complex<double>(0, -1));
                }
            }
            else
            {
                if (j % 2 == 0)
                {
                    // H21
                    h_lambda[index] = W_i_[index] * (lambda_[iat][0] + lambda_[iat][1] * std::complex<double>(0, 1));
                }
                else
                {
                    // H22
                    h_lambda[index] = W_i_[index] * (-lambda_[iat][2]);
                }
            }
        }
    }
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

    for (int irc = 0; irc < nloc_; irc++)
    {
        this->LM->Hloc2[irc] += h_lambda[irc];
    }
    ModuleBase::timer::tick("OperatorLambda", "contributeHk");
}

} // namespace hamilt
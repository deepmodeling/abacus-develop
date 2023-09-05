#include "sc_lambda_lcao.h"
#include "module_hamilt_lcao/module_sc/spin_constrain.h"

namespace hamilt
{

// contribute to HR is not needed.
template<typename T>
void OperatorScLambda<OperatorLCAO<T>>::contributeHR()
{
    return;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::set_nat(int nat_in)
{
    this->nat_ = nat_in;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::set_nloc(int nloc_in)
{
    this->nloc_ = nloc_in;
}

template<>
int OperatorScLambda<OperatorLCAO<std::complex<double>>>::get_nat()
{
    return this->nat_;
}

template<>
int OperatorScLambda<OperatorLCAO<std::complex<double>>>::get_nloc()
{
    return this->nloc_;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::set_iwt2iat(const int* iwt2iat_in)
{
    this->iwt2iat_.resize(this->nloc_);
    for (int i = 0; i < this->nloc_; i++)
    {
        this->iwt2iat_[i] = iwt2iat_in[i];
    }
}

template<>
const std::vector<int>& OperatorScLambda<OperatorLCAO<std::complex<double>>>::get_iwt2iat() const
{
    return this->iwt2iat_;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::set_lambda(const std::vector<ModuleBase::Vector3<double>>& lambda_in)
{
    if (lambda_in.size() != this->nat_)
    {
        ModuleBase::WARNING_QUIT("OperatorScLambda::set_lambda", "lambda_in size mismatch with nat");
    }
    this->lambda_ = lambda_in;
}

// get lambda
template<>
const std::vector<ModuleBase::Vector3<double>>& OperatorScLambda<OperatorLCAO<std::complex<double>>>::get_lambda() const
{
    return this->lambda_;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::cal_weight_func(const std::vector<std::complex<double>>& Sloc2)
{
    if (Sloc2.size() != this->nloc_ * this->nloc_)
    {
        ModuleBase::WARNING_QUIT("OperatorScLambda::cal_weight_func", "Sloc2 size mismatch with nloc * npol");
    }
    this->Wi_.resize(this->nloc_ * this->nloc_);
    for (int i = 0; i < this->nloc_; i++)
    {
        int iat = iwt2iat_[i];
        for (int j = 0; j < this->nloc_; j++)
        {
            int jat = iwt2iat_[j];
            this->Wi_[i * this->nloc_ + j]
                = (iat == jat) ? Sloc2[i * this->nloc_ + j] : Sloc2[i * this->nloc_ + j] * 0.5;
        }
    }
}

template<typename T>
void OperatorScLambda<OperatorLCAO<T>>::set_Wi(const std::vector<T>& Wi_in)
{
    if (Wi_in.size() != this->nloc_*this->nloc_)
    {
        ModuleBase::WARNING_QUIT("OperatorScLambda::set_Wi", "Wi size mismatch with nloc * nloc");
    }
    this->Wi_ = Wi_in;
}

// get_Wi
template<typename T>
const std::vector<T>& OperatorScLambda<OperatorLCAO<T>>::get_Wi() const
{
    return this->Wi_;
}

template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::cal_h_lambda(int ik, std::complex<double>* h_lambda)
{
    ModuleBase::TITLE("OperatorScLambda", "cal_h_lambda");
    ModuleBase::timer::tick("OperatorScLambda", "cal_h_lambda");
    // h_lambda = Wi * (lambda_x *sigma_x + lambda_y * sigma_y + lambda_z * sigma_z)
    // Pauli matrix is used implicitly
    // Pauli matrices
    // sigma_x = {{0, 1}, {1, 0}}
    // sigma_y = {{0, -i}, {i, 0}}
    // sigma_z = {{1, 0}, {0, -1}}
    // lambda_x * sigma_x + lambda_y * sigma_y + lambda_z * sigma_z
    // = {{lambda_z, lambda_x - i * lambda_y}, {lambda_x + i * lambda_y, -lambda_z}}
    for (int i = 0; i < nloc_; i++)
    {
        int iat = iwt2iat_[i];
        for (int j = 0; j < nloc_; j++)
        {
            int index = i * nloc_ + j;
            if (i % 2 == 0)
            {
                if (j % 2 == 0)
                {
                    // H11
                    h_lambda[index] = Wi_[index] * lambda_[iat][2];
                }
                else
                {
                    // H12
                    h_lambda[index] = Wi_[index] * (lambda_[iat][0] + lambda_[iat][1] * std::complex<double>(0, -1));
                }
            }
            else
            {
                if (j % 2 == 0)
                {
                    // H21
                    h_lambda[index] = Wi_[index] * (lambda_[iat][0] + lambda_[iat][1] * std::complex<double>(0, 1));
                }
                else
                {
                    // H22
                    h_lambda[index] = Wi_[index] * (-lambda_[iat][2]);
                }
            }
        }
    }
    ModuleBase::timer::tick("OperatorScLambda", "cal_h_lambda");
}

// contribute to Hk
template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorScLambda", "contributeHk");
    ModuleBase::timer::tick("OperatorScLambda", "contributeHk");
    //std::vector<std::complex<double>> h_lambda(this->nloc_ * this->nloc_);
    //this->cal_h_lambda(ik, &h_lambda[0]);

    //for (int irc = 0; irc < this->nloc_; irc++)
    //{
    //    this->LM->Hloc2[irc] += h_lambda[irc];
    //}
    SpinConstrain& sc = SpinConstrain::getInstance();
    std::cout << "nw in operator = " << sc.get_nw() << std::endl;
    std::cout << "OperatorScLambda contributeHk" << std::endl;
    ModuleBase::timer::tick("OperatorScLambda", "contributeHk");
}

template class OperatorScLambda<OperatorLCAO<std::complex<double>>>;

} // namespace hamilt
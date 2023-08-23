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
void OperatorLambda<OperatorLCAO<double>>::cal_h_lambda(int ik, double* h_lambda)
{
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::cal_h_lambda(int ik, std::complex<double>* h_lambda)
{
    ModuleBase::TITLE("OperatorLambda", "cal_h_lambda");
    ModuleBase::timer::tick("OperatorLambda", "cal_h_lambda");
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
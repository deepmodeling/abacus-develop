#include "sc_lambda_lcao.h"
#include "module_hamilt_lcao/module_sc/spin_constrain.h"
#include <algorithm>

namespace hamilt
{

// contribute to HR is not needed.
template<typename T>
void OperatorScLambda<OperatorLCAO<T>>::contributeHR()
{
    return;
}

// contribute to Hk
template<>
void OperatorScLambda<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorScLambda", "contributeHk");
    ModuleBase::timer::tick("OperatorScLambda", "contributeHk");
    SpinConstrain<double, psi::DEVICE_CPU>& sc = SpinConstrain<double, psi::DEVICE_CPU>::getInstance();
    std::vector<std::complex<double>> h_lambda(this->LM->ParaV->nloc);
    std::fill(h_lambda.begin(), h_lambda.end(), std::complex<double>(0, 0));
    sc.cal_weight_func(this->LM->Sloc2);
    sc.cal_h_lambda(&h_lambda[0]);
    for (int irc = 0; irc < this->LM->ParaV->nloc; irc++)
    {
        this->LM->Hloc2[irc] += h_lambda[irc];
    }
    //std::cout << "OperatorScLambda contributeHk" << std::endl;
    ModuleBase::timer::tick("OperatorScLambda", "contributeHk");
}

template class OperatorScLambda<OperatorLCAO<std::complex<double>>>;

} // namespace hamilt
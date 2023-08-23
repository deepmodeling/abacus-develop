#include "lambda_lcao.h"

namespace hamilt
{

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

// contribute to Hk.
template <>
void OperatorLambda<OperatorLCAO<double>>::contributeHk(int ik)
{
}

template <>
void OperatorLambda<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
}

template class OperatorLambda<OperatorLCAO<double>>;

template class OperatorLambda<OperatorLCAO<std::complex<double>>>;

} // namespace hamilt
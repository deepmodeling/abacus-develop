#include "meta_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class Meta<OperatorLCAO<double, double>>;

template class Meta<OperatorLCAO<std::complex<double>, double>>;

template class Meta<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<>
Meta<OperatorLCAO<double, double>>::~Meta()
{
}

template<>
Meta<OperatorLCAO<std::complex<double>, double>>::~Meta()
{
}

template<>
Meta<OperatorLCAO<std::complex<double>, std::complex<double>>>::~Meta()
{
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<double, double>>::contributeHR()
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<std::complex<double>, double>>::contributeHR()
{
    return;
}
//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHR()
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<double, double>>::contributeHk(int ik)
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    return;
}

template<>
void Meta<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    return;
}

}
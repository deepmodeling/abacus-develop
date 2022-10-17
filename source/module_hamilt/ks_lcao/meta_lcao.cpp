#include "meta_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "src_pw/global.h"

namespace hamilt
{

template class Meta<OperatorLCAO<double>>;

template class Meta<OperatorLCAO<std::complex<double>>>;

template<>
Meta<OperatorLCAO<double>>::~Meta()
{
}

template<>
Meta<OperatorLCAO<std::complex<double>>>::~Meta()
{
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<double>>::contributeHR()
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<std::complex<double>>>::contributeHR()
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<double>>::contributeHk(int ik)
{
    return;
}

//nothing to do in LCAO base for meta operator
template<>
void Meta<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    return;
}

}
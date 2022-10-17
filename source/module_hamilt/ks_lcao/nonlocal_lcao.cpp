#include "nonlocal_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

namespace hamilt
{

template class Nonlocal<OperatorLCAO<double>>;

template class Nonlocal<OperatorLCAO<std::complex<double>>>;

template<typename T>
void Nonlocal<OperatorLCAO<T>>::contributeHR()
{
    ModuleBase::TITLE("Nonlocal<OperatorLCAO>", "contributeHR");
    if(!this->HR_fixed_done)
    {
        ModuleBase::timer::tick("Nonlocal<OperatorLCAO>", "contributeHR");
        this->genH->calculate_NL_no(this->HR_pointer->data());
        ModuleBase::timer::tick("Nonlocal<OperatorLCAO>", "contributeHR");
        this->HR_fixed_done = true;
    }
    return;
}

}
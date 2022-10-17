#include "overlap_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

namespace hamilt
{

template class Overlap<OperatorLCAO<double>>;

template class Overlap<OperatorLCAO<std::complex<double>>>;

template<typename T>
void Overlap<OperatorLCAO<T>>::contributeHR()
{
    ModuleBase::TITLE("Overlap", "contributeHR");
    if(!this->SR_fixed_done)
    {
        ModuleBase::timer::tick("Overlap", "contributeHR");
        this->genH->calculate_S_no(this->SR_pointer->data());
        ModuleBase::timer::tick("Overlap", "contributeHR");
        this->SR_fixed_done = true;
    }
    return;
}

}
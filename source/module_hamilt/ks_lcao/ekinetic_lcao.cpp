#include "ekinetic_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

namespace hamilt
{

template class Ekinetic<OperatorLCAO<double>>;

template class Ekinetic<OperatorLCAO<std::complex<double>>>;

template<typename T>
void Ekinetic<OperatorLCAO<T>>::contributeHR()
{
    ModuleBase::TITLE("Ekinetic<OperatorLCAO>", "contributeHR");
    if(!this->HR_fixed_done)
    {
        ModuleBase::timer::tick("Ekinetic<OperatorLCAO>", "contributeHR");
        this->genH->calculate_T_no(this->HR_pointer->data());
        ModuleBase::timer::tick("Ekinetic<OperatorLCAO>", "contributeHR");
        this->HR_fixed_done = true;
    }
    return;
}

}
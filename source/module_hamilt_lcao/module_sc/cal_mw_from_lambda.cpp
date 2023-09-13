#include "spin_constrain.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::cal_mw_from_lambda(
    const std::vector<ModuleBase::Vector3<double>> delta_lambda,
    std::vector<ModuleBase::Vector3<double>> new_spin)
{
    ModuleBase::TITLE("SpinConstrain","cal_mw_from_lambda");
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
    this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, this->KS_SOLVER);
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
}

template class SpinConstrain<double, psi::DEVICE_CPU>;
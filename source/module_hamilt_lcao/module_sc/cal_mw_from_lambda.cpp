#include "spin_constrain.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/elecstate_lcao.h"

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::cal_mw_from_lambda(int i_step)
{
    ModuleBase::TITLE("SpinConstrain","cal_mw_from_lambda");
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
    this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, this->KS_SOLVER);
    elecstate::ElecStateLCAO* pelec_lcao = dynamic_cast<elecstate::ElecStateLCAO*>(this->pelec);
    this->cal_MW(i_step, *(this->LM), pelec_lcao->get_loc()->dm_k, GlobalC::ucell);
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
}

template class SpinConstrain<double, psi::DEVICE_CPU>;
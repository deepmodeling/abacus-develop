#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "spin_constrain.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_mw_from_lambda(int i_step)
{
    ModuleBase::TITLE("SpinConstrain","cal_mw_from_lambda");
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
    // diagonalization without update charge
    this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, this->KS_SOLVER, true);
    elecstate::ElecStateLCAO<std::complex<double>>* pelec_lcao
        = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
    this->pelec->calculate_weights();
    this->pelec->calEBand();
    if (this->KS_SOLVER == "genelpa" || this->KS_SOLVER == "scalapack_gvx" || this->KS_SOLVER == "lapack")
    {
        elecstate::cal_dm_psi(this->ParaV, pelec_lcao->wg, *(this->psi), *(pelec_lcao->get_DM()));
    }
    this->cal_MW(i_step, this->LM, GlobalC::ucell);
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
}
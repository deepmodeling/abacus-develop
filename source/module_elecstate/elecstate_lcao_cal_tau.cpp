#include "elecstate_lcao.h"
#include "module_hamilt_lcao/module_gint/new_grid_tech/gint_interface.h"
#include "elecstate_lcao_cal_tau.h"
#include "module_base/timer.h"

namespace elecstate
{

// calculate the kinetic energy density tau, multi-k case
void lcao_cal_tau_k(Gint_k* gint_k, 
                    Charge* charge)
{
    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
    }
    ModuleGint::cal_gint_tau(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);

    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");
    return;
}

// calculate the kinetic energy density tau, gamma-only case
void lcao_cal_tau_gamma(Gint_Gamma* gint_gamma,
                        Charge* charge)
{
    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(charge->kin_r[is], charge->nrxx);
    }
    ModuleGint::cal_gint_tau(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);

    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");
    return;
}
template <> 
void lcao_cal_tau<double>(Gint_Gamma* gint_gamma, 
                  Gint_k* gint_k, 
                  Charge* charge)
{
    lcao_cal_tau_gamma(gint_gamma, charge);
}
template <> 
void lcao_cal_tau<complex<double>>(Gint_Gamma* gint_gamma, 
                    Gint_k* gint_k, 
                    Charge* charge)
{
    lcao_cal_tau_k(gint_k, charge);
}

} // namespace elecstate
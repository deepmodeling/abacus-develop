#include "elecstate_lcao.h"
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
    Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau, PARAM.inp.nspin);
    gint_k->cal_gint(&inout1);

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
    Gint_inout inout1(charge->kin_r, Gint_Tools::job_type::tau, PARAM.inp.nspin);
    gint_gamma->cal_gint(&inout1);

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
#include "elecstate_lcao.h"

#include "module_base/timer.h"

namespace elecstate
{

// calculate the kinetic energy density tau, multi-k case
template <>
void ElecStateLCAO<std::complex<double>>::cal_tau(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
    }
    Gint_inout inout1(this->charge->kin_r, Gint_Tools::job_type::tau, PARAM.inp.nspin);
    this->gint_k->cal_gint(&inout1);

    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");
    return;
}

// calculate the kinetic energy density tau, gamma-only case
template <>
void ElecStateLCAO<double>::cal_tau(const psi::Psi<double>& psi)
{
    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
    }
    Gint_inout inout1(this->charge->kin_r, Gint_Tools::job_type::tau, PARAM.inp.nspin);
    this->gint_gamma->cal_gint(&inout1);

    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");
    return;
}
}
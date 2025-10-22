#include "elecstate_lcao.h"
#include "source_lcao/module_gint/gint_interface.h"

#include "source_base/timer.h"

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
    ModuleGint::cal_gint_tau(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);
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
    ModuleGint::cal_gint_tau(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);

    ModuleBase::timer::tick("ElecStateLCAO", "cal_tau");
    return;
}
}
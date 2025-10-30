#include "elecstate_lcao.h"

#include "cal_dm.h"
#include "source_base/timer.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/module_parameter/parameter.h"

#include "source_lcao/module_gint/gint_interface.h"

#include <vector>

namespace elecstate
{


template <>
double ElecStateLCAO<double>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<double>& sc = spinconstrain::SpinConstrain<double>::getScInstance();
    return sc.cal_escon();
}

template <>
double ElecStateLCAO<std::complex<double>>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

#ifdef __PEXSI
template <>
void ElecStateLCAO<double>::dm2Rho(std::vector<double*> pexsi_DM, std::vector<double*> pexsi_EDM)
{
    ModuleBase::timer::tick("ElecStateLCAO", "dm2Rho");

    int nspin = PARAM.inp.nspin;
    if (PARAM.inp.nspin == 4)
    {
        nspin = 1;
    }

    this->get_DM()->pexsi_EDM = pexsi_EDM;

    for (int is = 0; is < nspin; is++)
    {
        this->DM->set_DMK_pointer(is, pexsi_DM[is]);
    }
    DM->cal_DMR();

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is],
                                      this->charge->nrxx); // mohan 2009-11-10
    }

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    ModuleGint::cal_gint_rho(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->rho);
    if (XC_Functional::get_ked_flag())
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        ModuleGint::cal_gint_tau(this->DM->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "dm2Rho");
    return;
}

template <>
void ElecStateLCAO<std::complex<double>>::dm2rho(std::vector<std::complex<double>*> pexsi_DM,
                                                  std::vector<std::complex<double>*> pexsi_EDM)
{
    ModuleBase::WARNING_QUIT("ElecStateLCAO", "pexsi is not completed for multi-k case");
}

#endif

template class ElecStateLCAO<double>;               // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case

} // namespace elecstate

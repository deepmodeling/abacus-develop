#include "elecstate.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"

namespace elecstate
{

double ElecState::get_hartree_energy()
{
    return H_Hartree_pw::hartree_energy;
}

double ElecState::get_etot_efield()
{
    return Efield::etotefield;
}

double ElecState::get_etot_gatefield()
{
    return Gatefield::etotgatefield;
}

double ElecState::get_solvent_model_Ael()
{
    return GlobalC::solvent_model.cal_Ael(GlobalC::ucell, this->charge->nrxx, this->charge->nxyz);
}

double ElecState::get_solvent_model_Acav()
{
    return GlobalC::solvent_model.cal_Acav(GlobalC::ucell, this->charge->nxyz);
}

#ifdef __LCAO
double ElecState::get_spin_constrain_energy()
{
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
    return sc.cal_escon();
}

double ElecState::get_dftu_energy()
{
    return GlobalC::dftu.get_energy();
}
#endif

#ifdef __DEEPKS
double ElecState::get_deepks_E_delta()
{
    return GlobalC::ld.E_delta;
}
double ElecState::get_deepks_E_delta_band()
{
    return GlobalC::ld.e_delta_band;
}
#endif

} // namespace elecstate
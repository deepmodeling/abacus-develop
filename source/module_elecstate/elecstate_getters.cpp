#include "module_elecstate/elecstate_getters.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef __LCAO
#include "module_hamilt_lcao/module_dftu/dftu.h" //Quxin adds for DFT+U on 20201029
#endif
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_hamilt_general/module_surchem/surchem.h"

namespace elecstate
{

const int get_rhopw_nrxx()
{
    return GlobalC::rhopw->nrxx;
}
const int get_rhopw_nxyz()
{
    return GlobalC::rhopw->nxyz;
}
const double get_ucell_omega()
{
    return GlobalC::ucell.omega;
}
#ifdef __LCAO
const double get_dftu_energy()
{
    return GlobalC::dftu.get_energy();
}
#endif
#ifdef __DEEPKS
const double get_lcao_deepks_E_delta()
{
    return GlobalC::ld.E_delta;
}
const double get_lcao_deepks_e_delta_band()
{
    return GlobalC::ld.e_delta_band;
}
#endif
const double get_solvent_model_Ael()
{
    return GlobalC::solvent_model.cal_Ael(GlobalC::ucell, GlobalC::rhopw);
}
const double get_solvent_model_Acav()
{
    return GlobalC::solvent_model.cal_Acav(GlobalC::ucell, GlobalC::rhopw);
}
const double get_tot_magnetization()
{
    return GlobalC::ucell.magnet.tot_magnetization;
}
const double get_abs_magnetization()
{
    return GlobalC::ucell.magnet.abs_magnetization;
}
const double get_tot_magnetization_nc_x()
{
    return GlobalC::ucell.magnet.tot_magnetization_nc[0];
}
const double get_tot_magnetization_nc_y()
{
    return GlobalC::ucell.magnet.tot_magnetization_nc[1];
}
const double get_tot_magnetization_nc_z()
{
    return GlobalC::ucell.magnet.tot_magnetization_nc[2];
}

} // namespace elecstate

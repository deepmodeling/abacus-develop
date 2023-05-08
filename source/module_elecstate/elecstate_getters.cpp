#include "module_elecstate/elecstate_getters.h"

#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h"
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
#ifdef __EXX
#ifdef __LCAO
const double get_exx_lip_exx_energy()
{
    return GlobalC::exx_lip.get_exx_energy();
}
const bool get_exx_info_ri_real_number()
{
    return GlobalC::exx_info.info_ri.real_number;
}
const double get_exx_lri_double_Eexx()
{
    return GlobalC::exx_lri_double.Eexx;
}
const std::complex<double> get_exx_lri_complex_Eexx()
{
    return GlobalC::exx_lri_complex.Eexx;
}
const bool get_exx_info_global_cal_exx()
{
    return GlobalC::exx_info.info_global.cal_exx;
}
const double get_exx_info_global_hybrid_alpha()
{
    return GlobalC::exx_info.info_global.hybrid_alpha;
}
#endif // __LCAO
#endif // __EXX

} // namespace elecstate

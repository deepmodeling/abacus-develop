#include "pot_xc.h"

#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

void PotXC::cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff)
{
    ModuleBase::TITLE("PotXC", "cal_v_eff");
    ModuleBase::timer::tick("PotXC", "cal_v_eff");
    const int nrxx_current = chg->nrxx;
    
    //----------------------------------------------------------
    //  calculate the exchange-correlation potential
    //----------------------------------------------------------

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc_meta(nrxx_current, ucell->omega, ucell->tpiba, chg);
        *(this->etxc_) = std::get<0>(etxc_vtxc_v);
        *(this->vtxc_) = std::get<1>(etxc_vtxc_v);
        v_eff += std::get<2>(etxc_vtxc_v);
        *(this->vofk) = std::get<3>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("v_of_rho", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc(nrxx_current, chg, ucell);
        *(this->etxc_) = std::get<0>(etxc_vtxc_v);
        *(this->vtxc_) = std::get<1>(etxc_vtxc_v);
        v_eff += std::get<2>(etxc_vtxc_v);
    }
    ModuleBase::timer::tick("PotXC", "cal_v_eff");
}

} // namespace elecstate
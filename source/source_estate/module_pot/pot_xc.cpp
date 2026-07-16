#include "pot_xc.h"

#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

#ifdef USE_LIBXC
#include "source_hamilt/module_xc/libxc_abacus.h"
#endif

namespace elecstate
{

void PotXC::cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff)
{
    ModuleBase::TITLE("PotXC", "cal_veff");
    ModuleBase::timer::start("PotXC", "cal_veff");
    const int nrxx_current = chg->nrxx;
    
    //----------------------------------------------------------
    //  calculate the exchange-correlation potential
    //----------------------------------------------------------

    if (XC_Functional::get_ked_flag())
    {
#ifdef USE_LIBXC
        const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
#ifdef __EXX
        const double hse_omega = XC_Functional::get_hse_omega();
#else
        const double hse_omega = 0.0;
#endif
        const std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional_Libxc::v_xc_meta(XC_Functional::get_func_id(), nrxx_current, ucell->omega, ucell->tpiba, chg,
                                             PARAM.inp.nspin, hybrid_alpha, hse_omega);
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
        const double hybrid_alpha = XC_Functional::get_hybrid_alpha();
#ifdef __EXX
        const double hse_omega = XC_Functional::get_hse_omega();
#else
        const double hse_omega = 0.0;
#endif
        const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc(nrxx_current, chg, ucell,
                                  PARAM.inp.nspin,
                                  PARAM.globalv.domag,
                                  PARAM.globalv.domag_z,
                                  hybrid_alpha,
                                  hse_omega);
        *(this->etxc_) = std::get<0>(etxc_vtxc_v);
        *(this->vtxc_) = std::get<1>(etxc_vtxc_v);
        v_eff += std::get<2>(etxc_vtxc_v);
    }
    ModuleBase::timer::end("PotXC", "cal_veff");
}

} // namespace elecstate

#ifdef USE_LIBXC

#include "libxc_abacus.h"
#ifdef __EXX
#include "source_hamilt/module_xc/exx_info.h"
#endif

void XC_Functional_Libxc::xc_spin_libxc(
        const std::vector<int> &func_id,
        const double &rhoup, const double &rhodw,
        double &exc, double &vxcup, double &vxcdw)
{
    const std::vector<double> rho_ud = {rhoup, rhodw};
    exc = vxcup = vxcdw = 0.0;

    double hybrid_alpha = 0.0;
    double hse_omega = 0.0;
#ifdef __EXX
    hybrid_alpha = GlobalC::exx_info.info_global.hybrid_alpha;
    hse_omega = GlobalC::exx_info.info_global.hse_omega;
#endif
    std::vector<xc_func_type> funcs = XC_Functional_Libxc::init_func(
        /* func_id = */ func_id, 
        /* xc_polarized = */ XC_POLARIZED,
        /* hybrid_alpha = */ hybrid_alpha,
        /* hse_omega = */ hse_omega);

    for(xc_func_type &func : funcs)
    {
        double e = 0.0;
        std::vector<double> vxc_ud(2);
        if( func.info->family == XC_FAMILY_LDA)
        {
            // call Libxc function: xc_lda_exc_vxc
            xc_lda_exc_vxc( &func, 1, rho_ud.data(), &e, vxc_ud.data());
        }
        exc += e;
        vxcup += vxc_ud[0];
        vxcdw += vxc_ud[1];
    }

    XC_Functional_Libxc::finish_func(funcs);
}

#endif
#include "chg_precond.h"

#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace module_charge
{

void kerker_screen_recip(const MixingConfig& cfg,
                         ModulePW::PW_Basis* rhopw,
                         double tpiba,
                         std::complex<double>* drhog)
{
    ModuleBase::TITLE("module_charge", "kerker_screen_recip");

    if (cfg.mixing_gg0 <= 0.0 || cfg.mixing_beta <= 0.1)
    {
        return;
    }

    ModuleBase::timer::start("module_charge", "kerker_screen_recip");

    const int nspin = cfg.nspin;

    double fac = 0.0;
    double gg0 = 0.0;
    double amin = 0.0;

    /// consider a resize for mixing_angle
    int resize_tmp = 1;
    if (nspin == 4 && cfg.mixing_angle > 0)
    {
        resize_tmp = 2;
    }

    /// implement Kerker for density and magnetization separately
    for (int is = 0; is < nspin / resize_tmp; ++is)
    {
        const int is_idx = is * rhopw->npw;
        /// new mixing method only support nspin=2 not nspin=4
        if (is >= 1)
        {
            if (cfg.mixing_gg0_mag <= 0.0001 || cfg.mixing_beta_mag <= 0.1)
            {
#ifdef __DEBUG
                assert(is == 1); // make sure break works
#endif
                double is_mag = nspin - 1;
                //for (int ig = 0; ig < rhopw->npw * is_mag; ig++)
                //{
                //    drhog[is_idx + ig] *= 1;
                //}
                break;
            }
            fac = cfg.mixing_gg0_mag;
            amin = cfg.mixing_beta_mag;
        }
        else
        {
            fac = cfg.mixing_gg0;
            amin = cfg.mixing_beta;
        }

        gg0 = std::pow(fac * ModuleBase::BOHR_TO_A / tpiba, 2);

        const double gg0_amin = cfg.mixing_gg0_min / amin;

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig = 0; ig < rhopw->npw; ++ig)
        {
            double gg = rhopw->gg[ig];
            double filter_g = std::max(gg / (gg + gg0), gg0_amin);
            drhog[is_idx + ig] *= filter_g;
        }
    }

    ModuleBase::timer::end("module_charge", "kerker_screen_recip");
    return;
}

void kerker_screen_real(const MixingConfig& cfg,
                        ModulePW::PW_Basis* rhopw,
                        double tpiba,
                        double* drhor)
{
    ModuleBase::TITLE("module_charge", "kerker_screen_real");

    if (cfg.mixing_gg0 <= 0.0001 || cfg.mixing_beta <= 0.1)
    {
        return;
    }

    ModuleBase::timer::start("module_charge", "kerker_screen_real");

    const int nspin = cfg.nspin;
    assert(nspin == 1 || nspin == 2 || nspin == 4);

    /// consider a resize for mixing_angle
    const int spin_channels = (nspin == 4 && cfg.mixing_angle > 0) ? 2 : nspin;
    const bool screen_magnetization = cfg.mixing_gg0_mag > 0.0001 && cfg.mixing_beta_mag > 0.1;
    const int screened_channels = screen_magnetization ? spin_channels : 1;

    std::vector<std::complex<double>> drhog(rhopw->npw);
    std::vector<double> drhor_filter(rhopw->nrxx);
    for (int spin = 0; spin < screened_channels; ++spin)
    {
        double* residual = drhor + spin * rhopw->nrxx;
        const double screening = spin == 0 ? cfg.mixing_gg0 : cfg.mixing_gg0_mag;
        const double beta = spin == 0 ? cfg.mixing_beta : cfg.mixing_beta_mag;
        const double gg0 = std::pow(screening * ModuleBase::BOHR_TO_A / tpiba, 2);
        const double gg0_amin = cfg.mixing_gg0_min / beta;
        // Note after this process some G which is higher than Gmax will be filtered.
        // Thus we cannot use kerker_screen_recip(drhog.data()) directly after it.
        rhopw->real2recip(residual, drhog.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ig = 0; ig < rhopw->npw; ++ig)
        {
            const double gg = rhopw->gg[ig];
            const double filter = std::max(gg / (gg + gg0), gg0_amin);
            drhog[ig] *= (1 - filter);
        }
        // Subtracting the low-frequency correction preserves residual components above the basis cutoff.
        // A disabled magnetic filter is the identity and requires no transform.
        rhopw->recip2real(drhog.data(), drhor_filter.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int ir = 0; ir < rhopw->nrxx; ++ir)
        {
            residual[ir] -= drhor_filter[ir];
        }
    }

    ModuleBase::timer::end("module_charge", "kerker_screen_real");
    return;
}

} // namespace module_charge

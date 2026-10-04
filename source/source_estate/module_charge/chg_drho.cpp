#include "charge.h"
#include "chg_drho.h"
#include "chg_drho_detail.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"

#include <cassert>
#include <functional>

namespace module_charge
{

// Charge residual between chr->rho and chr->rho_save, normalized per electron.
double cal_drho(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega,
                const double tpiba)
{
    assert(chr != nullptr);
    ModuleBase::TITLE("module_charge", "cal_drho");
    ModuleBase::timer::start("module_charge", "cal_drho");
    const int nspin = cfg.nspin;
    assert(nspin==1 || nspin==2 || nspin==4);
    double drho = 0.0;

    if (cfg.scf_thr_type == 1)
    {
        for (int is = 0; is < nspin; ++is)
        {
            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho(r) to obtain rho(G).");
            chr->rhopw->real2recip(chr->rho[is], chr->rhog[is]);

            ModuleBase::GlobalFunc::NOTE("Perform FFT on rho_save(r) to obtain rho_save(G).");
            chr->rhopw->real2recip(chr->rho_save[is], chr->rhog_save[is]);
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the charge difference between rho(G) and rho_save(G)");
        std::vector<std::complex<double>> drhog(nspin * rhopw.npw);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 512)
#endif
        for (int is = 0; is < nspin; ++is)
        {
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                drhog[is * rhopw.npw + ig] = chr->rhog[is][ig] - chr->rhog_save[is][ig];
            }
        }

        ModuleBase::GlobalFunc::NOTE("Calculate the norm of the Residual std::vector: < R[rho] | R[rho_save] >");
        drho = module_charge::detail::inner_product_recip_rho(
            drhog.data(), drhog.data(), rhopw, cfg, omega, tpiba);
    }
    else
    {
        // Note: Maybe it is wrong.
        //       The inner_product_real function (L1-norm) is different from that (L2-norm) in mixing.
        for (int is = 0; is < nspin; is++)
        {
            if (is != 0 && is != 3 && cfg.domag_z)
            {
                continue;
            }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : drho)
#endif
            for (int ir = 0; ir < rhopw.nrxx; ir++)
            {
                drho += std::abs(chr->rho[is][ir] - chr->rho_save[is][ir]);
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(drho);
#endif
        assert(nelec != 0);
        assert(omega > 0);
        assert(rhopw.nxyz > 0);
        drho *= omega / static_cast<double>(rhopw.nxyz);
        drho /= nelec;
    }

    ModuleBase::timer::end("module_charge", "cal_drho");
    return drho;
}

// Kinetic-energy-density residual between chr->kin_r and chr->kin_r_save.
double cal_dkin(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega)
{
    assert(chr != nullptr);
    if (!(chr->meta_gga))
    {
        return 0.0;
    };
    ModuleBase::TITLE("module_charge", "cal_dkin");
    ModuleBase::timer::start("module_charge", "cal_dkin");
    double dkin = 0.0;

    // Get dkin from kin_r and kin_r_save for PW and LCAO both, which is different from drho.
    for (int is = 0; is < cfg.nspin; is++)
    {
        if (is != 0 && is != 3 && cfg.domag_z)
        {
            continue;
        }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : dkin)
#endif
        for (int ir = 0; ir < rhopw.nrxx; ir++)
        {
            dkin += std::abs(chr->kin_r[is][ir] - chr->kin_r_save[is][ir]);
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(dkin);
#endif
    assert(nelec != 0);
    assert(omega > 0);
    assert(rhopw.nxyz > 0);
    dkin *= omega / static_cast<double>(rhopw.nxyz);
    dkin /= nelec;

    ModuleBase::timer::end("module_charge", "cal_dkin");
    return dkin;
}

double inner_product_real(const double* rho1,
                          const double* rho2,
                          const ModulePW::PW_Basis& rhopw,
                          const MixingConfig& cfg)
{
    assert(rho1 != nullptr);
    assert(rho2 != nullptr);
    double rnorm = 0.0;
    // consider a resize for mixing_angle
    int resize_tmp = 1;
    if (cfg.nspin == 4 && cfg.mixing_angle > 0)
    {
        resize_tmp = 2;
    }

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : rnorm)
#endif
    for (int ir = 0; ir < rhopw.nrxx * cfg.nspin / resize_tmp; ++ir)
    {
        rnorm += rho1[ir] * rho2[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(rnorm);
#endif
    return rnorm;
}

} // namespace module_charge

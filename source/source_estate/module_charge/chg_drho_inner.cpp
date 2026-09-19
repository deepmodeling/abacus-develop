#include "chg_drho_detail.h"
#include "chg_mix_cfg.h"

#include <cassert>
#include <functional>
#include <vector>

#include "source_basis/module_pw/pw_basis.h"
#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

namespace module_charge
{
namespace detail
{

namespace
{

/// Coulomb-metric sum over G!=0 for a single spin channel
double coulomb_sum_single(const std::complex<double>* g1,
                         const std::complex<double>* g2,
                         const ModulePW::PW_Basis& rhopw,
                         const double fac)
{
    const int ig0 = rhopw.ig_gge0;
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < rhopw.npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(g1[ig]) * g2[ig]).real() / rhopw.gg[ig];
    }
    return sum * fac;
}

/// Non-magnetic case (nspin==1 or nspin==4 without domag)
double recip_rho_nspin1(const std::complex<double>* rho1,
                        const std::complex<double>* rho2,
                        const ModulePW::PW_Basis& rhopw,
                        const double fac)
{
    return coulomb_sum_single(rho1, rho2, rhopw, fac);
}

/// Collinear magnetic case (nspin==2)
double recip_rho_nspin2(const std::complex<double>* rho1,
                        const std::complex<double>* rho2,
                        const ModulePW::PW_Basis& rhopw,
                        const MixingConfig& cfg,
                        const double fac,
                        const double fac2)
{
    const int npw = rhopw.npw;
    const int ig0 = rhopw.ig_gge0;

    // (1) density part: |rho_up + rho_dn|^2 / G^2
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(rho1[ig] + rho1[ig + npw])
                * (rho2[ig] + rho2[ig + npw])).real() / rhopw.gg[ig];
    }
    sum *= fac;
    if (cfg.gamma_only_pw)
    {
        sum *= 2.0;
    }

    // (2) magnetization part: |rho_up - rho_dn|^2 (G=0 included)
    double sum2 = 0.0;
    if (ig0 >= 0)
    {
        sum2 += fac2 * (conj(rho1[ig0] - rho1[ig0 + npw])
                        * (rho2[ig0] - rho2[ig0 + npw])).real();
    }
    double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        mag += (conj(rho1[ig] - rho1[ig + npw])
                * (rho2[ig] - rho2[ig + npw])).real();
    }
    mag *= fac2;
    if (cfg.gamma_only_pw)
    {
        mag *= 2.0;
    }
    sum2 += mag;
    return sum + sum2;
}

/// Non-collinear magnetic case (nspin==4) with magnetization
double recip_rho_nspin4_mag(const std::complex<double>* rho1,
                            const std::complex<double>* rho2,
                            const ModulePW::PW_Basis& rhopw,
                            const MixingConfig& cfg,
                            const double fac,
                            const double fac2)
{
    const int npw = rhopw.npw;
    const int ig0 = rhopw.ig_gge0;

    // charge part
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(rho1[ig]) * rho2[ig]).real() / rhopw.gg[ig];
    }
    sum *= fac;

    // G=0 magnetization term
    if (ig0 > 0)
    {
        sum += fac2
               * ((conj(rho1[ig0 + npw]) * rho2[ig0 + npw]).real()
                  + (conj(rho1[ig0 + 2 * npw]) * rho2[ig0 + 2 * npw]).real()
                  + (conj(rho1[ig0 + 3 * npw]) * rho2[ig0 + 3 * npw]).real());
    }

    // G!=0 magnetization term
    double fac3 = fac2;
    if (cfg.gamma_only_pw)
    {
        fac3 *= 2.0;
    }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += fac3
               * ((conj(rho1[ig + npw]) * rho2[ig + npw]).real()
                  + (conj(rho1[ig + 2 * npw]) * rho2[ig + 2 * npw]).real()
                  + (conj(rho1[ig + 3 * npw]) * rho2[ig + 3 * npw]).real());
    }
    return sum;
}

/// Non-collinear with domag, traditional mixing (nspin==4, mixing_angle<=0)
double recip_hartree_nspin4_trad(const std::complex<double>* rhog1,
                                 const std::complex<double>* rhog2,
                                 const ModulePW::PW_Basis& rhopw,
                                 const MixingConfig& cfg,
                                 const double fac,
                                 const double fac2)
{
    const int npw = rhopw.npw;
    const int ig0 = rhopw.ig_gge0;

    // charge part
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
    }
    sum *= fac;

    // G=0 magnetization
    if (ig0 > 0)
    {
        sum += fac2
               * ((conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real()
                  + (conj(rhog1[ig0 + 2 * npw]) * rhog2[ig0 + 2 * npw]).real()
                  + (conj(rhog1[ig0 + 3 * npw]) * rhog2[ig0 + 3 * npw]).real());
    }

    // G!=0 magnetization
    double fac3 = fac2;
    if (cfg.gamma_only_pw)
    {
        fac3 *= 2.0;
    }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += fac3
               * ((conj(rhog1[ig + npw]) * rhog2[ig + npw]).real()
                  + (conj(rhog1[ig + 2 * npw]) * rhog2[ig + 2 * npw]).real()
                  + (conj(rhog1[ig + 3 * npw]) * rhog2[ig + 3 * npw]).real());
    }
    return sum;
}

/// Non-collinear with angle mixing (nspin==4, mixing_angle>0)
double recip_hartree_nspin4_angle(const std::complex<double>* rhog1,
                                  const std::complex<double>* rhog2,
                                  const ModulePW::PW_Basis& rhopw,
                                  const MixingConfig& cfg,
                                  const double fac,
                                  const double fac2)
{
    const int npw = rhopw.npw;
    const int ig0 = rhopw.ig_gge0;

    // charge part (same as rho only)
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
    }
    sum *= fac;

    // G=0 |m| term
    if (ig0 > 0)
    {
        sum += fac2 * (conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real();
    }

    // G!=0 |m| term
    double fac3 = fac2;
    if (cfg.gamma_only_pw)
    {
        fac3 *= 2.0;
    }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += fac3 * (conj(rhog1[ig + npw]) * rhog2[ig + npw]).real();
    }
    return sum;
}

/// Collinear magnetic case for hartree metric (nspin==2)
double recip_hartree_nspin2(const std::complex<double>* rhog1,
                            const std::complex<double>* rhog2,
                            const ModulePW::PW_Basis& rhopw,
                            const MixingConfig& cfg,
                            const double fac,
                            const double fac2)
{
    const int npw = rhopw.npw;
    const int ig0 = rhopw.ig_gge0;

    // charge density part
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        sum += (conj(rhog1[ig]) * rhog2[ig]).real() / rhopw.gg[ig];
    }
    sum *= fac;
    if (cfg.gamma_only_pw)
    {
        sum *= 2.0;
    }

    // magnetization part (G=0 included)
    double sum2 = 0.0;
    if (ig0 >= 0)
    {
        sum2 += fac2 * (conj(rhog1[ig0 + npw]) * rhog2[ig0 + npw]).real();
    }
    double mag = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : mag)
#endif
    for (int ig = 0; ig < npw; ++ig)
    {
        if (ig == ig0)
        {
            continue;
        }
        mag += (conj(rhog1[ig + npw]) * rhog2[ig + npw]).real();
    }
    mag *= fac2;
    if (cfg.gamma_only_pw)
    {
        mag *= 2.0;
    }
    sum2 += mag;
    return sum + sum2;
}

} // anonymous namespace

double inner_product_recip_rho(const std::complex<double>* rho1,
                               const std::complex<double>* rho2,
                               const ModulePW::PW_Basis& rhopw,
                               const MixingConfig& cfg,
                               const double omega,
                               const double tpiba)
{
    assert(rho1 != nullptr);
    assert(rho2 != nullptr);
    assert(cfg.nspin == 1 || cfg.nspin == 2 || cfg.nspin == 4);
    ModuleBase::TITLE("Charge_Mixing", "recip_rho");
    ModuleBase::timer::start("Charge_Mixing", "recip_rho");

    const int nspin = cfg.nspin;
    std::vector<const std::complex<double>*> rhog1(nspin);
    std::vector<const std::complex<double>*> rhog2(nspin);
    for (int is = 0; is < nspin; is++)
    {
        rhog1[is] = rho1 + is * rhopw.npw;
        rhog2[is] = rho2 + is * rhopw.npw;
    }

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba * tpiba);
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;

    switch (nspin)
    {
    case 1:
        sum += recip_rho_nspin1(rhog1[0], rhog2[0], rhopw, fac);
        break;
    case 2:
        sum += recip_rho_nspin2(rhog1[0], rhog2[0], rhopw, cfg, fac, fac2);
        break;
    case 4:
        if (!cfg.domag && !cfg.domag_z)
        {
            sum += recip_rho_nspin1(rhog1[0], rhog2[0], rhopw, fac);
        }
        else
        {
            sum += recip_rho_nspin4_mag(rhog1[0], rhog2[0], rhopw, cfg, fac, fac2);
        }
        break;
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif
    sum *= omega * 0.5;

    ModuleBase::timer::end("Charge_Mixing", "recip_rho");
    return sum;
}

} // namespace detail

double inner_product_recip_hartree(const std::complex<double>* rhog1,
                                   const std::complex<double>* rhog2,
                                   const ModulePW::PW_Basis& rhopw,
                                   const MixingConfig& cfg,
                                   const double omega,
                                   const double tpiba)
{
    assert(rhog1 != nullptr);
    assert(rhog2 != nullptr);
    assert(cfg.nspin == 1 || cfg.nspin == 2 || cfg.nspin == 4);
    ModuleBase::TITLE("Charge_Mixing", "recip_hartree");
    ModuleBase::timer::start("Charge_Mixing", "recip_hartree");

    static const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (tpiba * tpiba);
    static const double fac2 = ModuleBase::e2 * ModuleBase::FOUR_PI / (ModuleBase::TWO_PI * ModuleBase::TWO_PI);

    double sum = 0.0;

    if (cfg.nspin == 1)
    {
        sum += detail::coulomb_sum_single(rhog1, rhog2, rhopw, fac);
    }
    else if (cfg.nspin == 2)
    {
        sum += detail::recip_hartree_nspin2(rhog1, rhog2, rhopw, cfg, fac, fac2);
    }
    else if (cfg.nspin == 4)
    {
        if (!cfg.domag && !cfg.domag_z)
        {
            sum += detail::coulomb_sum_single(rhog1, rhog2, rhopw, fac);
        }
        else if (cfg.mixing_angle <= 0)
        {
            sum += detail::recip_hartree_nspin4_trad(rhog1, rhog2, rhopw, cfg, fac, fac2);
        }
        else
        {
            sum += detail::recip_hartree_nspin4_angle(rhog1, rhog2, rhopw, cfg, fac, fac2);
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum);
#endif
    sum *= omega * 0.5;

    ModuleBase::timer::end("Charge_Mixing", "recip_hartree");
    return sum;
}

} // namespace module_charge

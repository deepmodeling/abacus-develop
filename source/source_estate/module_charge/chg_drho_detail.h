#ifndef CHG_DRHO_DETAIL_H
#define CHG_DRHO_DETAIL_H

// Internal reciprocal inner product for the charge residual (cal_drho).
// Not part of the public module_charge API: only chg_drho.cpp and the
// charge mixing unit test are expected to include this header.

#include <complex>

#include "chg_mix_cfg.h"

namespace ModulePW
{
class PW_Basis;
}

namespace module_charge
{
namespace detail
{

/**
 * @brief Coulomb-metric reciprocal inner product of the charge residual.
 *
 * @param rho1 first reciprocal-space vector
 * @param rho2 second reciprocal-space vector
 * @param rhopw plane-wave basis supplying npw/gg and the G=0 index
 * @param cfg mixing config (spin channels, gamma-only and magnetism flags)
 * @param omega cell volume
 * @param tpiba 2*pi/lattice constant
 * @return pooled Coulomb-metric inner product
 */
double inner_product_recip_rho(const std::complex<double>* rho1,
                               const std::complex<double>* rho2,
                               const ModulePW::PW_Basis& rhopw,
                               const MixingConfig& cfg,
                               const double omega,
                               const double tpiba);

} // namespace detail
} // namespace module_charge

#endif // CHG_DRHO_DETAIL_H

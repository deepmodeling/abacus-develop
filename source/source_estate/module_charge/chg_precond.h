#ifndef CHG_PRECOND_H
#define CHG_PRECOND_H

// Stateless Kerker preconditioning kernels extracted from Charge_Mixing.
// Every input (grid, geometry, mixing config) is passed explicitly; the
// functions do not read Charge_Mixing members or PARAM/GlobalV.

#include <complex>

#include "chg_mix_cfg.h"

namespace ModulePW
{
class PW_Basis;
}

namespace module_charge
{

/**
 * @brief Apply Kerker screening in reciprocal space.
 *
 * Multiplies drhog[is*npw + ig] by max(gg/(gg+gg0), gg0_min/amin) per spin
 * channel, where gg0 is derived from cfg.mixing_gg0 (density) or
 * cfg.mixing_gg0_mag (magnetization). Early return if Kerker is disabled.
 *
 * @param cfg mixing config (spin count, betas, gg0s, angle, gg0_min)
 * @param rhopw plane-wave basis supplying npw and gg[]
 * @param tpiba 2*pi/lattice constant used to convert gg0 to atomic units
 * @param drhog[in,out] reciprocal-space density residual, length nspin*npw
 */
void kerker_screen_recip(const MixingConfig& cfg,
                         ModulePW::PW_Basis* rhopw,
                         double tpiba,
                         std::complex<double>* drhog);

/**
 * @brief Apply Kerker screening in real space via FFT.
 *
 * Forward-transforms drhor to drhog, applies (1 - filter_g) in reciprocal
 * space, backward-transforms the filtered residual, and subtracts it from
 * drhor in place. Early return if Kerker is disabled.
 *
 * @param cfg mixing config (spin count, betas, gg0s, angle, gg0_min)
 * @param rhopw plane-wave basis supplying npw, nrxx, gg[], real2recip/recip2real
 * @param tpiba 2*pi/lattice constant used to convert gg0 to atomic units
 * @param drhor[in,out] real-space density residual, length nspin*nrxx
 */
void kerker_screen_real(const MixingConfig& cfg,
                        ModulePW::PW_Basis* rhopw,
                        double tpiba,
                        double* drhor);

} // namespace module_charge

#endif // CHG_PRECOND_H

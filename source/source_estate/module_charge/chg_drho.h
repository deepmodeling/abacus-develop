#ifndef CHG_DRHO_H
#define CHG_DRHO_H

// Stateless residual kernels extracted from Charge_Mixing. Every input
// (grid, geometry, mixing config) is passed explicitly; the functions do
// not read Charge_Mixing members or PARAM/GlobalV.

#include <complex>

#include "chg_mix_cfg.h"

namespace ModulePW
{
class PW_Basis;
}

class Charge;

namespace module_charge
{

/**
 * @brief Charge residual between chr->rho and chr->rho_save, normalized per electron.
 *
 * @param chr charge object supplying rho/rho_save (and rhog buffers for the reciprocal case)
 * @param nelec number of electrons, used to normalize the real-space residual
 * @param rhopw plane-wave basis supplying the real/reciprocal grid sizes
 * @param cfg mixing config (nspin, scf_thr_type and magnetization flags select the loops)
 * @param omega cell volume, used to normalize the real-space residual
 * @param tpiba 2*pi/lattice constant, used by the reciprocal metric
 * @return pooled residual value
 */
double cal_drho(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega,
                const double tpiba);

/**
 * @brief Kinetic-energy-density residual between chr->kin_r and chr->kin_r_save.
 *
 * @param chr charge object supplying kin_r/kin_r_save
 * @param nelec number of electrons, used to normalize the residual
 * @param rhopw plane-wave basis supplying the real-space grid size
 * @param cfg mixing config (nspin and magnetization flags select the loops)
 * @param omega cell volume, used to normalize the residual
 * @return pooled residual value
 */
double cal_dkin(Charge* chr,
                const double nelec,
                const ModulePW::PW_Basis& rhopw,
                const MixingConfig& cfg,
                const double omega);

/**
 * @brief Inner product of two real-space vectors used in real-space mixing.
 *
 * @param rho1 first real-space vector
 * @param rho2 second real-space vector
 * @param rhopw plane-wave basis supplying the real-space grid size
 * @param cfg mixing config (nspin and mixing_angle select the loop bound)
 * @return pooled inner product
 */
double inner_product_real(const double* rho1,
                          const double* rho2,
                          const ModulePW::PW_Basis& rhopw,
                          const MixingConfig& cfg);

/**
 * @brief Hartree-like reciprocal inner product used in charge mixing.
 *
 * @param rhog1 first reciprocal-space vector
 * @param rhog2 second reciprocal-space vector
 * @param rhopw plane-wave basis supplying npw/gg and the G=0 index
 * @param cfg mixing config (spin channels, gamma-only and angle flags)
 * @param omega cell volume
 * @param tpiba 2*pi/lattice constant
 * @return pooled Hartree inner product
 */
double inner_product_recip_hartree(const std::complex<double>* rhog1,
                                   const std::complex<double>* rhog2,
                                   const ModulePW::PW_Basis& rhopw,
                                   const MixingConfig& cfg,
                                   const double omega,
                                   const double tpiba);

} // namespace module_charge

#endif // CHG_DRHO_H

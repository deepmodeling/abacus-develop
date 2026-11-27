#ifndef CHG_TAU_H
#define CHG_TAU_H

// Mixing of the kinetic energy density (tau) in reciprocal space.
// The implementation lives in chg_tau.cpp; this header exposes the
// free function so that Charge_Mixing::mix_rho_recip can call it
// without pulling in the broader chg_rho_detail.h helpers.

#include "charge.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/module_mixing/plain_mixing.h"

namespace module_charge
{
namespace detail
{

/**
 * @brief Mix kinetic energy density in reciprocal space.
 *        Handles the double-grid split/merge for the smooth and
 *        high-frequency parts, DIIS mixing of the smooth part, and
 *        plain mixing of the high-frequency part.
 * @param chr pointer to Charge object (must have kin_r/kin_r_save)
 * @param nspin number of spins
 * @param double_grid whether double grid is used
 * @param rhopw smooth grid
 * @param rhodpw dense grid (same as rhopw when double_grid is off)
 * @param mixing DIIS mixing object
 * @param tau_mdata mixing data for tau
 * @param mixing_highf plain mixing for high-frequency part (may be null when double_grid is off)
 */
void mix_tau_recip(Charge* chr,
                   const int nspin,
                   const bool double_grid,
                   ModulePW::PW_Basis* rhopw,
                   ModulePW::PW_Basis* rhodpw,
                   Base_Mixing::Mixing* mixing,
                   Base_Mixing::Mixing_Data& tau_mdata,
                   Base_Mixing::Plain_Mixing* mixing_highf);

} // namespace detail
} // namespace module_charge

#endif // CHG_TAU_H

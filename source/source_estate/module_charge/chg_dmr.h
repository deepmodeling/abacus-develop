#ifndef CHG_DMR_H
#define CHG_DMR_H

// Stateless real-space density-matrix (DMR) mixing kernels extracted from
// Charge_Mixing. The mixing history (Mixing_Data) and the Mixing strategy
// object remain owned by Charge_Mixing and are passed in explicitly; the
// functions do not read Charge_Mixing members or PARAM/GlobalV. The kernels
// work on raw contiguous double buffers (one per spin channel) and do not
// depend on DensityMatrix or HContainer; the caller extracts those buffers
// from its LCAO containers.

#include <vector>

#include "chg_mix_cfg.h"

namespace Base_Mixing
{
class Mixing;
class Mixing_Data;
} // namespace Base_Mixing

namespace module_charge
{

/**
 * @brief Allocate the mixing buffer for the real-space density matrix and
 *        clear its history.
 *
 * The buffer cannot be allocated in Charge_Mixing::set_mixing(): its length
 * nnr (number of non-zero R-matrix elements) is only known after
 * DensityMatrix::init_dmr(), which runs later in beforescf().
 *
 * @param mixing    mixing strategy object, non-null
 * @param mdata     mixing history buffer for DMR, resized and reset in place
 * @param nnr       number of real-space density-matrix elements per spin, > 0
 * @param cfg       mixing config (nspin and scf_thr_type select the path)
 */
void init_mixing_dmr(Base_Mixing::Mixing* mixing,
                     Base_Mixing::Mixing_Data& mdata,
                     const int nnr,
                     const MixingConfig& cfg);

/**
 * @brief Mix the real-space density matrix (LCAO calculations only).
 *
 * For nspin == 1/4 the single spin channel is mixed directly; for nspin == 2
 * the up/down channels are transformed into charge/magnetization channels,
 * mixed with independent betas, and transformed back.
 *
 * @param dmr_out  writable DMR buffers, one per spin channel, each of length
 *                 nnr; mixed results are written back through these pointers
 * @param dmr_in   DMR buffers saved at the previous mixing step, one per spin
 *                 channel, each of length nnr (read-only)
 * @param nnr      number of DMR elements per spin channel, > 0
 * @param mixing   mixing strategy object, non-null
 * @param mdata    DMR mixing history buffer
 * @param cfg      mixing config (nspin and the two mixing betas)
 */
void mix_dmr(const std::vector<double*>& dmr_out,
             const std::vector<const double*>& dmr_in,
             const int nnr,
             Base_Mixing::Mixing* mixing,
             Base_Mixing::Mixing_Data& mdata,
             const MixingConfig& cfg);

} // namespace module_charge

#endif // CHG_DMR_H

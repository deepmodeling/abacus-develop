#ifndef USPP_DENSITY_H_
#define USPP_DENSITY_H_

#include "source_base/module_device/device.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/vnl_pw.h"

#include <complex>
#include <memory>
#include <vector>

namespace elecstate
{
/** @brief Size of the canonical [spin][atom][packed projector pair] storage. */
int uspp_becsum_size(const UnitCell& ucell, const pseudopot_cell_vnl& ppcell, const int nspin);

/** @brief Accumulate explicitly weighted collinear USPP states with reusable device buffers. */
template <typename T, typename Device>
class UsppProjector
{
  public:
    /** @brief Bind cell/projectors and allocate storage for a block of nbands states. */
    UsppProjector(const UnitCell& ucell, const pseudopot_cell_vnl& ppcell, const int nbands);
    ~UsppProjector();

    /**
     * @brief Add weighted projector products to the caller's packed coefficients.
     * @param coefficients Device-resident PW coefficients, accessed with the supplied stride.
     * @param stride Storage leading dimension between successive states.
     * @param npw Number of active plane waves on this pool rank.
     * @param spin Scalar/collinear density channel receiving the contribution.
     * @param weights One explicit weight for each state in the block.
     * @param becsum Host output; existing contributions are retained.
     * All ranks in the owning POOL_WORLD must participate in the overlap reduction.
     * The current contract has no spinor component layout or cross-spin products.
     * Supporting nspin=4 requires extending this interface and constructing charge/magnetization projector weights.
     */
    void accumulate(const int ik,
                    const T* coefficients,
                    const int stride,
                    const int npw,
                    const int spin,
                    const std::vector<double>& weights,
                    std::vector<double>* becsum);

  private:
    const UnitCell& ucell_;
    const pseudopot_cell_vnl& ppcell_;
    const int nbands_;
    class Workspace;
    std::unique_ptr<Workspace> work_;
};

/**
 * @brief Add USPP augmentation to caller-owned dense-grid reciprocal densities.
 * becsum uses i <= j pairs, with the off-diagonal multiplicity already included.
 * The Q(G) tables supply volume normalization; state weights are dimensionless.
 * The current scalar/collinear contraction does not acquire spinor support by setting nspin=4.
 * Spinor support requires charge/magnetization coefficients and a compatible Q/projector representation,
 * including review of the packed-pair convention for SOC.
 */
void add_uspp_density(const UnitCell& ucell,
                      const pseudopot_cell_vnl& ppcell,
                      const ModulePW::PW_Basis& basis,
                      const int nspin,
                      const std::vector<double>& becsum,
                      std::complex<double>** rhog);
} // namespace elecstate

#endif

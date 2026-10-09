#ifndef TD_HISTORY_PW_H
#define TD_HISTORY_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_psi/psi.h"

#include <memory>

namespace pw
{

/** @brief Non-owning device pointers to the smooth local potentials. */
template <typename Real>
struct TDPotentialView
{
    const Real* veff;
    const Real* vofk;
};

/** @brief Converged time-step history and predictor-corrector potential storage. */
template <typename T, typename Device>
class TDHistoryPW
{
  public:
    /** @brief Prepare storage after the potential grids have been updated. */
    void prepare(const elecstate::Potential& potential, const bool needs_ked);

    /** @brief Use the current potential for the predictor, then average with the saved potential. */
    TDPotentialView<typename GetTypeReal<T>::type> prop_pot(elecstate::Potential& potential,
                                                          const bool needs_ked,
                                                          const bool corrector);

    /** @brief Save the converged state; the caller owns the convergence decision. */
    void save(const psi::Psi<T, Device>& current, elecstate::Potential& potential, const bool needs_ked);

    /** @brief Access the fixed origin of every corrector solve in the next step. */
    const psi::Psi<T, Device>& previous() const;

  private:
    using Real = typename GetTypeReal<T>::type;
    std::unique_ptr<psi::Psi<T, Device>> previous_;
    ct::Tensor veff_prev_;
    ct::Tensor veff_mid_;
    ct::Tensor vofk_prev_;
    ct::Tensor vofk_mid_;
};

} // namespace pw
#endif

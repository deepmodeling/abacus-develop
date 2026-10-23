#ifndef TD_HISTORY_PW_H
#define TD_HISTORY_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_psi/psi.h"

#include <memory>

namespace pw
{

/** @brief Accepted endpoint wavefunctions and potential snapshots. */
template <typename T, typename Device>
class TDHistoryPW
{
  public:
    /** @brief Prepare storage after the potential grids have been updated. */
    void prepare(const elecstate::Potential& potential);

    /** @brief Borrow the immutable accepted left-endpoint potential. */
    const typename GetTypeReal<T>::type* left_potential() const;
    unsigned long long generation() const
    {
        return generation_;
    }

    /** @brief Save the accepted endpoint; the caller owns the convergence decision. */
    void save(const psi::Psi<T, Device>& current, elecstate::Potential& potential);

    /** @brief Access the fixed origin of every corrector solve in the next step. */
    const psi::Psi<T, Device>& previous() const;

  private:
    using Real = typename GetTypeReal<T>::type;
    std::unique_ptr<psi::Psi<T, Device>> previous_;
    ct::Tensor veff_prev_;
    unsigned long long generation_ = 0;
};

} // namespace pw
#endif

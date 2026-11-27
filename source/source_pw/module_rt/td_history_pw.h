#ifndef TD_HISTORY_PW_H
#define TD_HISTORY_PW_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_estate/fp_energy.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_psi/psi.h"

#include <memory>
#include <vector>

namespace pw
{

/** @brief Electronic data at an accepted endpoint; orbitals and potential stay in TDHistoryPW. */
struct TDElecSnap
{
    elecstate::fenergy energy;
    ModuleBase::matrix eigenvalues;
    std::vector<double> density;
};

/** @brief Copy energy, band expectations and density without advancing orbital history. */
TDElecSnap capture_td(const elecstate::fenergy& energy, const ModuleBase::matrix& eigenvalues, const Charge& charge);

/** @brief Restore electronic data on the unchanged density grid and refresh its mixing backup. */
void restore_td(const TDElecSnap& state, elecstate::fenergy& energy, ModuleBase::matrix& eigenvalues, Charge& charge);

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

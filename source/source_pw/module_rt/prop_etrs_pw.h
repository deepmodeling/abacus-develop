#ifndef PW_PROP_ETRS_PW_H
#define PW_PROP_ETRS_PW_H
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/module_linalg/exp_action.h"
#include "source_pw/module_rt/td_prop_pw.h"

namespace pw
{
template <typename T, typename Device>
class ETRSPropPW final : public PWProp<T, Device>
{
  private:
    const hsolver::diag_comm_info pool_comm_;
    std::unique_ptr<hsolver::ExpAction<T, Device>> exponential_;
    ct::Tensor half_psi_;
    const int npol_;
    const std::vector<int> global_k_;
    std::ostream& log_;
    bool half_ready_ = false;
    unsigned long long left_gen_ = 0;
    void report(const hsolver::ExpResult& result, int ik, int elec_step, int scf_iter, bool log_detail) const;

  public:
    ETRSPropPW(const hsolver::ExpControls& exp_ctrl,
               int npol,
               const std::vector<int>& global_k,
               const hsolver::diag_comm_info& comm,
               std::ostream& log);
    void begin_trial() override;
    void discard_trial() override;
    void invalidate_basis() override;
    void evolve(TDHamiltPW<T, Device>& hamilt,
                const TDHistoryPW<T, Device>& history,
                elecstate::Potential& potential,
                const elecstate::TDFieldManager& field,
                psi::Psi<T, Device>* current,
                int elec_step,
                int scf_iter,
                bool log_detail) override;
};
} // namespace pw
#endif

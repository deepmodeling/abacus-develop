#ifndef PW_PROP_ETRS_PW_H
#define PW_PROP_ETRS_PW_H
#include "source_hsolver/diag_comm_info.h"
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
    hsolver::TDOrbOpts orb_opts_;
    std::ostream& log_;
    bool half_ready_ = false;
    unsigned long long left_gen_ = 0;
    void report(const hsolver::ExpResult& result, const PropArgs<T, Device>& args, int ik) const;

  public:
    ETRSPropPW(const TDPropOpts& prop_opts, const hsolver::diag_comm_info& comm, std::ostream& log);
    void begin_trial() override;
    void discard_trial() override;
    void invalidate_basis() override;
    void evolve(PropArgs<T, Device>& args, psi::Psi<T, Device>* current) override;
};
} // namespace pw
#endif

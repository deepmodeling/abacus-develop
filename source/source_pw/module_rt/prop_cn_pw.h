#ifndef PW_PROP_CN_PW_H
#define PW_PROP_CN_PW_H
#include "source_pw/module_rt/td_prop_pw.h"

namespace pw
{
template <typename T, typename Device>
class CNPropPW final : public PWProp<T, Device>
{
  private:
    hsolver::CNSolverPW<T, Device> solver_;
    ct::Tensor veff_mid_;

  public:
    CNPropPW(const ModulePW::PW_Basis_K& basis,
             const hsolver::PWLinearOptions& cn_opts,
             const hsolver::diag_comm_info& comm,
             std::ostream& log);
    void begin_trial() override
    {
    }
    void discard_trial() override
    {
    }
    void invalidate_basis() override;
    void evolve(PropArgs<T, Device>& args, psi::Psi<T, Device>* current) override;
};
} // namespace pw
#endif

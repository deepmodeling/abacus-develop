#ifndef HSOLVER_EXP_TAYLOR_H
#define HSOLVER_EXP_TAYLOR_H
#include "source_hsolver/module_linalg/exp_action.h"

namespace hsolver
{
template <typename T, typename Device>
class ExpTaylor final : public ExpAction<T, Device>
{
  private:
    ct::Tensor term_;
    ct::Tensor image_;
    ExpResult evaluate(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* half, T* full);

  public:
    ExpTaylor(const ExpControls& exp_ctrl, const diag_comm_info& comm);
    ExpResult apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output) override;
    ExpResult apply_pair(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* half, T* full) override;
};
} // namespace hsolver
#endif

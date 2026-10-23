#ifndef HSOLVER_EXP_LANCZOS_H
#define HSOLVER_EXP_LANCZOS_H
#include "source_hsolver/module_rt_pw/exp_action.h"

namespace hsolver
{
template <typename T, typename Device>
class ExpLanczos final : public ExpAction<T, Device>
{
  private:
    ct::Tensor basis_;
    ct::Tensor image_;
    ct::Tensor packed_;
    ct::Tensor hpacked_;
    ExpResult attempt(const HSOperator<T, Device>& op,
                 const OrbLayout& layout,
                 double alpha,
                 double rel_tol,
                 const T* input,
                 T* output);
    ExpResult interval(const HSOperator<T, Device>& op,
                  const OrbLayout& layout,
                  double alpha,
                  double rel_tol,
                  const T* input,
                  T* output,
                  int depth);

  public:
    ExpLanczos(const ExpOpts& exp_opts, const diag_comm_info& comm);
    ExpResult apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output) override;
};
} // namespace hsolver
#endif

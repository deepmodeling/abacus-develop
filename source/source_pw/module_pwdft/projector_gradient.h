#ifndef PW_PROJECTOR_GRADIENT_H
#define PW_PROJECTOR_GRADIENT_H
#include "source_pw/module_pwdft/vnl_pw.h"
#include "source_base/module_container/ATen/core/tensor.h"

namespace hamilt
{
/** @brief One-k-point gradient workspace with versioned radial tables. */
template <typename Real, typename Device>
class ProjectorGradient
{
  private:
    using Complex = std::complex<Real>;
    ct::Tensor radial_;
    ct::Tensor derivative_;
    ct::Tensor metadata_;
    ct::Tensor momentum_;
    ct::Tensor structure_;
    size_t version_ = 0;
    std::vector<Real> host_q_;
    std::vector<int> host_metadata_;

    template <typename Value>
    void reserve(ct::Tensor* tensor, const int64_t count);

  public:
    /** @brief Refresh the shifted gradient in native precision and device storage. */
    void calculate(pseudopot_cell_vnl* pp, const UnitCell& cell, const ModulePW::PW_Basis_K& basis,
                   const int ik, const ModuleBase::Vector3<double>& A, Complex* output);
};
} // namespace hamilt
#endif

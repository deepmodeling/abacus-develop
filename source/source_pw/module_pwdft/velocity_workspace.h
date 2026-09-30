#ifndef HAMILT_VELOCITY_WORKSPACE_H
#define HAMILT_VELOCITY_WORKSPACE_H

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace hamilt
{
/** @brief Reusable projection storage and device contraction for a velocity operator. */
template <typename Real, typename Device>
class VelocityWorkspace
{
  private:
    using Complex = std::complex<Real>;
    ct::Tensor coefficients_;
    std::vector<Complex> host_;
    int capacity_ = 0;

  public:
    /** @brief Reserve four input and four output projector blocks. */
    Complex* prepare(const int count);

    /** @brief Reduce projections and contract D on their native device. */
    void contract(const UnitCell& cell, const pseudopot_cell_vnl& pp, const int spin,
                  const int bands, const Real scale, const ModulePW::PW_Basis_K& basis, Complex* buffer);
};
} // namespace hamilt
#endif

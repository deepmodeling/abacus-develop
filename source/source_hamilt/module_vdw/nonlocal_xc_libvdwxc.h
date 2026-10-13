#ifndef ABACUS_NONLOCAL_XC_LIBVDWXC_H
#define ABACUS_NONLOCAL_XC_LIBVDWXC_H

#include "nonlocal_xc.h"

#include <array>

#ifdef __LIBVDWXC
// libvdwxc's public C header includes <complex.h> inside an extern "C"
// block. Pre-include the standard header here so libc++ does not parse its
// C++ templates with C linkage.
#include <complex.h>
#include <vdwxc.h>
#endif

namespace vdw
{

#ifdef __LIBVDWXC
/**
 * @brief Thin serial/PW adapter for libvdwxc's rVV10 evaluator.
 *
 * libvdwxc owns the nonlocal kernel, interpolation and convolution. ABACUS
 * owns density construction, FFT gradients/divergence, units and XC
 * bookkeeping. The first adapter intentionally rejects distributed layouts;
 * the libvdwxc MPI grid/communicator contract must be validated separately.
 */
class LibvdwxcNonlocalXC final : public NonlocalXC
{
  public:
    LibvdwxcNonlocalXC(const ModulePW::PW_Basis& pw, double b, double c);
    ~LibvdwxcNonlocalXC() override;

    NonlocalXCResult evaluate(const ModulePW::PW_Basis& pw,
                              const std::vector<double>& total_density,
                              const std::vector<double>& valence_density) const override;

  private:
    vdwxc_data data_ = nullptr;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    double omega_ = 0.0;
    std::array<double, 9> cell_{};
};
#endif

} // namespace vdw

#endif // ABACUS_NONLOCAL_XC_LIBVDWXC_H

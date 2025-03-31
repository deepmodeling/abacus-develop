#include "pw_basis.h"
#include "module_base/timer.h"
#include "module_basis/module_pw/kernels/pw_op.h"
namespace ModulePW
{
#if (defined(__CUDA) || defined(__ROCM))
template <typename FPTYPE>
void PW_Basis::real2recip_gpu(const FPTYPE* in,
                             std::complex<FPTYPE>* out,
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    base_device::DEVICE_GPU* ctx;
    // base_device::memory::synchronize_memory_op<std::complex<FPTYPE>,
    //                                            base_device::DEVICE_GPU,
    //                                            base_device::DEVICE_GPU>()(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
    //                                                                       in,
    //                                                                       this->nrxx);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                                  this->nxyz,
                                                                  add,
                                                                  factor,
                                                                  this->ig2isz,
                                                                  this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                  out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}
template <typename FPTYPE>
void PW_Basis::real2recip_gpu(const std::complex<FPTYPE>* in,
                             std::complex<FPTYPE>* out,
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    base_device::DEVICE_GPU* ctx;
    base_device::memory::synchronize_memory_op<std::complex<FPTYPE>,
                                               base_device::DEVICE_GPU,
                                               base_device::DEVICE_GPU>()(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                          in,
                                                                          this->nrxx);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2isz,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);
    ModuleBase::timer::tick(this->classname, "real_to_recip gpu");
}

template <typename FPTYPE>
void PW_Basis::recip2real_gpu(const std::complex<FPTYPE>* in,
                             FPTYPE* out,
                             const bool add,
                             const FPTYPE factor) const
{
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    base_device::DEVICE_GPU* ctx;
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<FPTYPE>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);

    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                        this->ig2isz,
                                                        in,
                                                        this->fft_bundle.get_auxr_3d_data<FPTYPE>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                    this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                  add,
                                                                  factor,
                                                                  this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                  out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}
template <typename FPTYPE>
    void PW_Basis::recip2real_gpu(const std::complex<FPTYPE> *in,
                                 std::complex<FPTYPE> *out,
                                 const bool add,
                                 const FPTYPE factor) const
{
    base_device::DEVICE_GPU* ctx;
    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
    assert(this->gamma_only == false);
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<double>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);

    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                         this->ig2isz,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<FPTYPE>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                    this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);

    ModuleBase::timer::tick(this->classname, "recip_to_real gpu");
}
template void PW_Basis::real2recip_gpu<double>(const double* in,
                                              std::complex<double>* out,
                                              const bool add,
                                              const double factor) const;
template void PW_Basis::real2recip_gpu<float>(const float* in,
                                              std::complex<float>* out,
                                              const bool add,
                                              const float factor) const;

template void PW_Basis::real2recip_gpu<double>(const std::complex<double>* in,
                                              std::complex<double>* out,
                                              const bool add,
                                              const double factor) const;
template void PW_Basis::real2recip_gpu<float>(const std::complex<float>* in,
                                              std::complex<float>* out,
                                              const bool add,
                                              const float factor) const;

template void PW_Basis::recip2real_gpu<double>(const std::complex<double>* in,
                                              double* out,
                                              const bool add,
                                              const double factor) const;
template void PW_Basis::recip2real_gpu<float>(const std::complex<float>* in,
                                              float* out,
                                              const bool add,
                                              const float factor) const;

template void PW_Basis::recip2real_gpu<double>(const std::complex<double>* in,
                                              std::complex<double>* out,
                                              const bool add,
                                              const double factor) const;
template void PW_Basis::recip2real_gpu<float>(const std::complex<float>* in,
                                              std::complex<float>* out,
                                              const bool add,
                                              const float factor) const;

#endif
} // namespace ModulePW
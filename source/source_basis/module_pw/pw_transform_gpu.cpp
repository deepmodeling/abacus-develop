#include "source_base/timer.h"
#include "source_basis/module_pw/kernels/pw_op.h"
#include "pw_basis.h"

#if defined(__ROCM) && !defined(__CUDA)
#error "PW GPU batch transforms are not implemented for ROCm in this merge."
#endif

namespace ModulePW
{
#if (defined(__CUDA) || defined(__ROCM))
template <typename FPTYPE>
void PW_Basis::real2recip_gpu(const FPTYPE* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip gpu");
    assert(this->poolnproc == 1);
    const size_t size = this->nrxx;
    base_device::memory::cast_memory_op<std::complex<FPTYPE>, FPTYPE,base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        in,
        size);

    this->fft_bundle.fft3D_forward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                                   this->nxyz,
                                                                   add,
                                                                   factor,
                                                                   this->ig2ixyz_gpu,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);
    ModuleBase::timer::end(this->classname, "real_to_recip gpu");
}
template <typename FPTYPE>
void PW_Basis::real2recip_gpu(const std::complex<FPTYPE>* in,
                              std::complex<FPTYPE>* out,
                              const bool add,
                              const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip gpu");
    assert(this->poolnproc == 1);
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
                                                                   this->ig2ixyz_gpu,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);
    ModuleBase::timer::end(this->classname, "real_to_recip gpu");
}

template <typename FPTYPE>
void PW_Basis::recip2real_gpu(const std::complex<FPTYPE>* in, FPTYPE* out, const bool add, const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real gpu");
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<FPTYPE>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                         this->ig2ixyz_gpu,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<FPTYPE>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                    this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);

    ModuleBase::timer::end(this->classname, "recip_to_real gpu");
}
template <typename FPTYPE>
void PW_Basis::recip2real_gpu(const std::complex<FPTYPE>* in,
                              std::complex<FPTYPE>* out,
                              const bool add,
                              const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real gpu");
    assert(this->poolnproc == 1);
    // ModuleBase::GlobalFunc::ZEROS(fft_bundle.get_auxr_3d_data<double>(), this->nxyz);
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
        0,
        this->nxyz);

    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>()(npw,
                                                         this->ig2ixyz_gpu,
                                                         in,
                                                         this->fft_bundle.get_auxr_3d_data<FPTYPE>());
    this->fft_bundle.fft3D_backward(this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                    this->fft_bundle.get_auxr_3d_data<FPTYPE>());

    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>()(this->nrxx,
                                                                   add,
                                                                   factor,
                                                                   this->fft_bundle.get_auxr_3d_data<FPTYPE>(),
                                                                   out);

    ModuleBase::timer::end(this->classname, "recip_to_real gpu");
}

// Batch FFT transforms (GPU only)
template <typename FPTYPE, typename Device,
          typename std::enable_if<std::is_same<Device, base_device::DEVICE_GPU>::value, int>::type>
void PW_Basis::real_to_recip_batch(const Device* ctx,
                                   const std::complex<FPTYPE>* in_batch,
                                   std::complex<FPTYPE>* out_batch,
                                   int batch_count,
                                   const bool add,
                                   const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "real_to_recip_batch gpu");

    // Check if batch FFT available
    if (!this->fft_bundle.is_batch_fft_available<FPTYPE>() ||
        batch_count > this->fft_bundle.get_batch_size<FPTYPE>())
    {
        // Fallback to sequential transforms
        for (int ib = 0; ib < batch_count; ++ib)
        {
            this->real2recip_gpu(
                in_batch + ib * this->nrxx,
                out_batch + ib * this->npw,
                add,
                factor);
        }
        ModuleBase::timer::end(this->classname, "real_to_recip_batch gpu");
        return;
    }

    // Batch FFT path
    assert(this->gamma_only == false);

    // Get batch output buffer from FFT bundle
    std::complex<FPTYPE>* batch_out = this->fft_bundle.get_batch_output_buffer<FPTYPE>();
    const int nxyz = this->nxyz;

    // Perform batch forward FFT directly on input (out-of-place FFT doesn't modify input)
    // const_cast is safe because cuFFT out-of-place C2C transform only reads from input
    this->fft_bundle.fft3D_forward_batch(const_cast<std::complex<FPTYPE>*>(in_batch), batch_out, batch_count);

#if defined(__CUDA)
    // Extract results for each transform in batch.
    set_real_to_recip_output_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        this->npw,
        nxyz,
        this->npw,      // out_stride = npw (same as count for PW_Basis)
        batch_count,
        add,
        factor,
        this->ig2ixyz_gpu,
        batch_out,
        out_batch);
#endif

    ModuleBase::timer::end(this->classname, "real_to_recip_batch gpu");
}

template <typename FPTYPE, typename Device,
          typename std::enable_if<std::is_same<Device, base_device::DEVICE_GPU>::value, int>::type>
void PW_Basis::recip_to_real_batch(const Device* ctx,
                                   const std::complex<FPTYPE>* in_batch,
                                   std::complex<FPTYPE>* out_batch,
                                   int batch_count,
                                   const bool add,
                                   const FPTYPE factor) const
{
    ModuleBase::timer::start(this->classname, "recip_to_real_batch gpu");

    // Check if batch FFT available
    if (!this->fft_bundle.is_batch_fft_available<FPTYPE>() ||
        batch_count > this->fft_bundle.get_batch_size<FPTYPE>())
    {
        // Fallback to sequential transforms
        for (int ib = 0; ib < batch_count; ++ib)
        {
            this->recip2real_gpu(
                in_batch + ib * this->npw,
                out_batch + ib * this->nrxx,
                add,
                factor);
        }
        ModuleBase::timer::end(this->classname, "recip_to_real_batch gpu");
        return;
    }

    // Batch FFT path
    assert(this->gamma_only == false);

    // Get batch buffers
    std::complex<FPTYPE>* batch_in = this->fft_bundle.get_batch_input_buffer<FPTYPE>();
    std::complex<FPTYPE>* batch_out = this->fft_bundle.get_batch_output_buffer<FPTYPE>();
    const int nxyz = this->nxyz;

    // Zero batch input buffer
    base_device::memory::set_memory_op<std::complex<FPTYPE>, base_device::DEVICE_GPU>()(
        batch_in,
        0,
        batch_count * nxyz);

    // Set up 3D FFT input grids for all batch elements using batched kernel
    set_3d_fft_box_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        this->npw,
        nxyz,
        this->npw,          // in_stride = npw (same as count for PW_Basis)
        this->ig2ixyz_gpu,
        batch_count,
        in_batch,
        batch_in);

    // Perform batch backward FFT
    this->fft_bundle.fft3D_backward_batch(batch_in, batch_out, batch_count);

#if defined(__CUDA)
    // Extract results with the scalar factor shared by all batch elements.
    set_recip_to_real_output_op<FPTYPE, base_device::DEVICE_GPU>().operator_batch(
        this->nrxx,     // Per-element size (not total!)
        nxyz,           // Input stride
        batch_count,
        add,
        factor,
        batch_out,
        out_batch);
#endif

    ModuleBase::timer::end(this->classname, "recip_to_real_batch gpu");
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

// Template instantiations for batch FFT methods
template void PW_Basis::real_to_recip_batch<float, base_device::DEVICE_GPU, 0>(
    const base_device::DEVICE_GPU*, const std::complex<float>*,
    std::complex<float>*, int, const bool, const float) const;

template void PW_Basis::real_to_recip_batch<double, base_device::DEVICE_GPU, 0>(
    const base_device::DEVICE_GPU*, const std::complex<double>*,
    std::complex<double>*, int, const bool, const double) const;

template void PW_Basis::recip_to_real_batch<float, base_device::DEVICE_GPU, 0>(
    const base_device::DEVICE_GPU*, const std::complex<float>*,
    std::complex<float>*, int, const bool, const float) const;

template void PW_Basis::recip_to_real_batch<double, base_device::DEVICE_GPU, 0>(
    const base_device::DEVICE_GPU*, const std::complex<double>*,
    std::complex<double>*, int, const bool, const double) const;

#endif
} // namespace ModulePW

#include "fft_bundle.h"

#include "source_base/module_device/device.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/tool_quit.h"

#include <cassert>
#if defined(__CUDA)
#include "fft_cuda.h"
#endif
#if defined(__ROCM)
#include "fft_rocm.h"
#endif
#if defined(__DSP)
#include "fft_dsp.h"
#endif
#include <stdexcept>
template <typename FFT_BASE, typename... Args>
std::unique_ptr<FFT_BASE> make_unique(Args&&... args)
{
    return std::unique_ptr<FFT_BASE>(new FFT_BASE(std::forward<Args>(args)...));
}
namespace ModuleBase
{
FFT_Bundle::~FFT_Bundle()
{
    this->clear();
}

void FFT_Bundle::setfft(std::string device_in, std::string precision_in)
{
    this->device = device_in;
    this->precision = precision_in;
}

void FFT_Bundle::initfft(int nx_in,
                         int ny_in,
                         int nz_in,
                         int lixy_in,
                         int rixy_in,
                         int ns_in,
                         int nplane_in,
                         int nproc_in,
                         bool gamma_only_in,
                         bool xprime_in,
                         bool mpifft_in)
{
    assert(this->device == "cpu" || this->device == "gpu" || this->device == "dsp");
    assert(this->precision == "single" || this->precision == "double" || this->precision == "mixing");

    if (this->precision == "single" || this->precision == "mixing")
    {
        float_flag = true;
        if (this->precision == "mixing")
        {
            double_flag = true;
        }
#if not defined(__ENABLE_FLOAT_FFTW)
        if (this->device == "cpu")
        {
            ModuleBase::WARNING_QUIT("FFT_Bundle", "Please enable float fftw in the cmake to use float fft");
        }
#endif
    }
    else if (this->precision == "double")
    {
        double_flag = true;
    }else{
        ModuleBase::WARNING_QUIT("FFT_Bundle", "Please set the precision to single or double or mixing");
    }
#if defined(__DSP)
    if (device == "dsp")
    {
        if (float_flag)
        {
            ModuleBase::WARNING_QUIT("device", "now dsp fft is not supported for the float type");
        }
        auto dsp_fft = make_unique<FFT_DSP<double>>();
        dsp_fft->cluster_id = this->dsp_cluster_id_;
        fft_double = std::move(dsp_fft);
        fft_double->initfft(nx_in, ny_in, nz_in);
    }else
#endif
    if (device == "cpu")
    {
        if (float_flag)
        {
            fft_float = make_unique<FFT_CPU<float>>(this->fft_mode);
            fft_float
                ->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in, nplane_in, nproc_in, gamma_only_in, xprime_in);
        }
        if (double_flag)
        {
            fft_double = make_unique<FFT_CPU<double>>(this->fft_mode);
            fft_double
                ->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in, nplane_in, nproc_in, gamma_only_in, xprime_in);
        }
    }else if (device == "gpu")
    {
#if defined(__ROCM)
        fft_float = make_unique<FFT_ROCM<float>>();
        fft_float->initfft(nx_in, ny_in, nz_in);
        fft_double = make_unique<FFT_ROCM<double>>();
        fft_double->initfft(nx_in, ny_in, nz_in);
#elif defined(__CUDA)
        fft_float = make_unique<FFT_CUDA<float>>();
        fft_float->initfft(nx_in, ny_in, nz_in);
        fft_double = make_unique<FFT_CUDA<double>>();
        fft_double->initfft(nx_in, ny_in, nz_in);
#endif

#if defined(__ENABLE_FLOAT_FFTW)
        if (float_flag)
        {
            fft_float_cpu = make_unique<FFT_CPU<float>>(this->fft_mode);
            fft_float_cpu
                ->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in, nplane_in, nproc_in, gamma_only_in, xprime_in);
        }
#endif
        if (double_flag)
        {
            fft_double_cpu = make_unique<FFT_CPU<double>>(this->fft_mode);
            fft_double_cpu
                ->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in, nplane_in, nproc_in, gamma_only_in, xprime_in);
        }
    }else{
        ModuleBase::WARNING_QUIT("FFT_Bundle", "Please set the device to cpu or gpu or dsp");
    }
}

void FFT_Bundle::setupFFT()
{
    if (double_flag)
    {
        fft_double->setupFFT();
        if (fft_double_cpu != nullptr)
        {
            fft_double_cpu->setupFFT();
        }
    }
    if (float_flag)
    {
        fft_float->setupFFT();
        if (fft_float_cpu != nullptr)
        {
            fft_float_cpu->setupFFT();
        }
    }
}

void FFT_Bundle::clearFFT()
{
    if (double_flag)
    {
        fft_double->cleanFFT();
        if (fft_double_cpu != nullptr)
        {
            fft_double_cpu->cleanFFT();
        }
    }
    if (float_flag)
    {
        fft_float->cleanFFT();
        if (fft_float_cpu != nullptr)
        {
            fft_float_cpu->cleanFFT();
        }
    }
}
void FFT_Bundle::clear()
{
    this->clearFFT();
    if (double_flag)
    {
        fft_double->clear();
        if (fft_double_cpu != nullptr)
        {
            fft_double_cpu->clear();
        }
    }
    if (float_flag)
    {
        fft_float->clear();
        if (fft_float_cpu != nullptr)
        {
            fft_float_cpu->clear();
        }
    }
}

void FFT_Bundle::resource_handler(const int flag) const
{
    if (this->device=="dsp")
    {
        if (double_flag)
        {
            fft_double->resource_handler(flag);
        }
        if (float_flag)
        {
            fft_float->resource_handler(flag);
        }
    }
}
template <>
void FFT_Bundle::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    // Use CPU FFT for xy operations (required for non-templated CPU code paths)
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftxyfor(in, out);
    }
    else
    {
        fft_float->fftxyfor(in, out);
    }
}
template <>
void FFT_Bundle::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftxyfor(in, out);
    }
    else
    {
        fft_double->fftxyfor(in, out);
    }
}

template <>
void FFT_Bundle::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftzfor(in, out);
    }
    else
    {
        fft_float->fftzfor(in, out);
    }
}
template <>
void FFT_Bundle::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftzfor(in, out);
    }
    else
    {
        fft_double->fftzfor(in, out);
    }
}

template <>
void FFT_Bundle::fftxybac(std::complex<float>* in, std::complex<float>* out) const
{
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftxybac(in, out);
    }
    else
    {
        fft_float->fftxybac(in, out);
    }
}
template <>
void FFT_Bundle::fftxybac(std::complex<double>* in, std::complex<double>* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftxybac(in, out);
    }
    else
    {
        fft_double->fftxybac(in, out);
    }
}

template <>
void FFT_Bundle::fftzbac(std::complex<float>* in, std::complex<float>* out) const
{
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftzbac(in, out);
    }
    else
    {
        fft_float->fftzbac(in, out);
    }
}
template <>
void FFT_Bundle::fftzbac(std::complex<double>* in, std::complex<double>* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftzbac(in, out);
    }
    else
    {
        fft_double->fftzbac(in, out);
    }
}

template <>
void FFT_Bundle::fftxyr2c(float* in, std::complex<float>* out) const
{
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftxyr2c(in, out);
    }
    else
    {
        fft_float->fftxyr2c(in, out);
    }
}
template <>
void FFT_Bundle::fftxyr2c(double* in, std::complex<double>* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftxyr2c(in, out);
    }
    else
    {
        fft_double->fftxyr2c(in, out);
    }
}

template <>
void FFT_Bundle::fftxyc2r(std::complex<float>* in, float* out) const
{
    if (fft_float_cpu != nullptr)
    {
        fft_float_cpu->fftxyc2r(in, out);
    }
    else
    {
        fft_float->fftxyc2r(in, out);
    }
}
template <>
void FFT_Bundle::fftxyc2r(std::complex<double>* in, double* out) const
{
    if (fft_double_cpu != nullptr)
    {
        fft_double_cpu->fftxyc2r(in, out);
    }
    else
    {
        fft_double->fftxyc2r(in, out);
    }
}

template <>
void FFT_Bundle::fft3D_forward(std::complex<float>* in,
                               std::complex<float>* out) const
{
    fft_float->fft3D_forward(in, out);
}
template <>
void FFT_Bundle::fft3D_forward(std::complex<double>* in,
                               std::complex<double>* out) const
{
    fft_double->fft3D_forward(in, out);
}

template <>
void FFT_Bundle::fft3D_backward(std::complex<float>* in,
                                std::complex<float>* out) const
{
    fft_float->fft3D_backward(in, out);
}
template <>
void FFT_Bundle::fft3D_backward(std::complex<double>* in,
                                std::complex<double>* out) const
{
    fft_double->fft3D_backward(in, out);
}

// access the real space data - always use CPU FFT for non-templated code paths
template <>
float* FFT_Bundle::get_rspace_data() const
{
    if (fft_float_cpu != nullptr)
    {
        return fft_float_cpu->get_rspace_data();
    }
    return fft_float->get_rspace_data();
}
template <>
double* FFT_Bundle::get_rspace_data() const
{
    if (fft_double_cpu != nullptr)
    {
        return fft_double_cpu->get_rspace_data();
    }
    return fft_double->get_rspace_data();
}

template <>
std::complex<float>* FFT_Bundle::get_auxr_data() const
{
    if (fft_float_cpu != nullptr)
    {
        return fft_float_cpu->get_auxr_data();
    }
    return fft_float->get_auxr_data();
}
template <>
std::complex<double>* FFT_Bundle::get_auxr_data() const
{
    if (fft_double_cpu != nullptr)
    {
        return fft_double_cpu->get_auxr_data();
    }
    return fft_double->get_auxr_data();
}

template <>
std::complex<float>* FFT_Bundle::get_auxg_data() const
{
    if (fft_float_cpu != nullptr)
    {
        return fft_float_cpu->get_auxg_data();
    }
    return fft_float->get_auxg_data();
}
template <>
std::complex<double>* FFT_Bundle::get_auxg_data() const
{
    if (fft_double_cpu != nullptr)
    {
        return fft_double_cpu->get_auxg_data();
    }
    return fft_double->get_auxg_data();
}

// get_auxr_3d_data - uses primary (GPU) FFT for GPU operations
template <>
std::complex<float>* FFT_Bundle::get_auxr_3d_data() const
{
    return fft_float->get_auxr_3d_data();
}
template <>
std::complex<double>* FFT_Bundle::get_auxr_3d_data() const
{
    return fft_double->get_auxr_3d_data();
}

// ============================================================================
// Batch FFT Implementation
// ============================================================================

void FFT_Bundle::setupBatchFFT()
{
    // Only setup batch FFT for GPU devices with a true batch size.
    if (this->device != "gpu" || this->batch_size <= 1)
    {
        return;
    }

#if defined(__CUDA)
    if (double_flag && fft_double != nullptr)
    {
        // Cast to FFT_CUDA to access batch FFT methods
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr)
        {
            fft_cuda_double->setupBatchFFT(this->batch_size);
        }
    }
    if (float_flag && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr)
        {
            fft_cuda_float->setupBatchFFT(this->batch_size);
        }
    }
#endif
    // ROCm batch FFT is not supported yet.
}

template <>
void FFT_Bundle::fft3D_forward_batch(std::complex<float>* in_batch,
                                     std::complex<float>* out_batch,
                                     int batch_count) const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr && fft_cuda_float->is_batch_fft_ready())
        {
            fft_cuda_float->fft3D_forward_batch(in_batch, out_batch, batch_count);
            return;
        }
    }
#endif
    throw std::runtime_error("FFT_Bundle::fft3D_forward_batch: batch FFT is unavailable for float precision");
}

template <>
void FFT_Bundle::fft3D_forward_batch(std::complex<double>* in_batch,
                                     std::complex<double>* out_batch,
                                     int batch_count) const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr && fft_cuda_double->is_batch_fft_ready())
        {
            fft_cuda_double->fft3D_forward_batch(in_batch, out_batch, batch_count);
            return;
        }
    }
#endif
    throw std::runtime_error("FFT_Bundle::fft3D_forward_batch: batch FFT is unavailable for double precision");
}

template <>
void FFT_Bundle::fft3D_backward_batch(std::complex<float>* in_batch,
                                      std::complex<float>* out_batch,
                                      int batch_count) const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr && fft_cuda_float->is_batch_fft_ready())
        {
            fft_cuda_float->fft3D_backward_batch(in_batch, out_batch, batch_count);
            return;
        }
    }
#endif
    throw std::runtime_error("FFT_Bundle::fft3D_backward_batch: batch FFT is unavailable for float precision");
}

template <>
void FFT_Bundle::fft3D_backward_batch(std::complex<double>* in_batch,
                                      std::complex<double>* out_batch,
                                      int batch_count) const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr && fft_cuda_double->is_batch_fft_ready())
        {
            fft_cuda_double->fft3D_backward_batch(in_batch, out_batch, batch_count);
            return;
        }
    }
#endif
    throw std::runtime_error("FFT_Bundle::fft3D_backward_batch: batch FFT is unavailable for double precision");
}

template <>
bool FFT_Bundle::is_batch_fft_available<float>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr)
        {
            return fft_cuda_float->is_batch_fft_ready();
        }
    }
#endif
    return false;
}

template <>
bool FFT_Bundle::is_batch_fft_available<double>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr)
        {
            return fft_cuda_double->is_batch_fft_ready();
        }
    }
#endif
    return false;
}

template <>
int FFT_Bundle::get_batch_size<float>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        // Only return batch_size if batch FFT is actually initialized
        if (fft_cuda_float != nullptr && fft_cuda_float->is_batch_fft_ready())
        {
            return fft_cuda_float->get_batch_size();
        }
    }
#endif
    return 0;
}

template <>
int FFT_Bundle::get_batch_size<double>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        // Only return batch_size if batch FFT is actually initialized
        if (fft_cuda_double != nullptr && fft_cuda_double->is_batch_fft_ready())
        {
            return fft_cuda_double->get_batch_size();
        }
    }
#endif
    return 0;
}

template <>
std::complex<float>* FFT_Bundle::get_batch_input_buffer<float>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr)
        {
            return fft_cuda_float->get_batch_input_buffer();
        }
    }
#endif
    return nullptr;
}

template <>
std::complex<double>* FFT_Bundle::get_batch_input_buffer<double>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr)
        {
            return fft_cuda_double->get_batch_input_buffer();
        }
    }
#endif
    return nullptr;
}

template <>
std::complex<float>* FFT_Bundle::get_batch_output_buffer<float>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_float != nullptr)
    {
        auto fft_cuda_float = dynamic_cast<FFT_CUDA<float>*>(fft_float.get());
        if (fft_cuda_float != nullptr)
        {
            return fft_cuda_float->get_batch_output_buffer();
        }
    }
#endif
    return nullptr;
}

template <>
std::complex<double>* FFT_Bundle::get_batch_output_buffer<double>() const
{
#if defined(__CUDA)
    if (this->device == "gpu" && fft_double != nullptr)
    {
        auto fft_cuda_double = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
        if (fft_cuda_double != nullptr)
        {
            return fft_cuda_double->get_batch_output_buffer();
        }
    }
#endif
    return nullptr;
}

} // namespace ModuleBase

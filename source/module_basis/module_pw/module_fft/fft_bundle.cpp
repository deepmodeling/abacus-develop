#include <cassert>
#include "fft_bundle.h"

#include "module_base/module_device/device.h"
#include "module_base/module_device/memory_op.h"
#if defined(__CUDA)
#include "fft_cuda.h"
#endif
#if defined(__ROCM)
#include "fft_rocm.h"
#endif

template<typename FFT_BASE, typename... Args>
std::unique_ptr<FFT_BASE> make_unique(Args &&... args)
{
    return std::unique_ptr<FFT_BASE>(new FFT_BASE(std::forward<Args>(args)...));
}
namespace ModulePW
{
FFT_Bundle::~FFT_Bundle()
{
    this->clear();
}

void FFT_Bundle::setfft(std::string device_in,std::string precision_in)
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
                         bool xprime_in , 
                         bool mpifft_in)
{
    assert(this->device=="cpu" || this->device=="gpu");
    assert(this->precision=="single" || this->precision=="double" || this->precision=="mixing");

    if (this->precision=="single")
    {
        #ifndef __ENABLE_FLOAT_FFTW
        float_define = false;
        #endif
        float_flag = float_define;
        double_flag = true;
    }
    if (this->precision=="double")
    {
        double_flag = true;
    }

    if (device=="cpu")
    {
        fft_float = make_unique<FFT_CPU<float>>(this->fft_mode);
        fft_double = make_unique<FFT_CPU<double>>(this->fft_mode);
        if (float_flag)
        {
            fft_float->initfft(nx_in,
                               ny_in,
                               nz_in,
                               lixy_in,
                               rixy_in,
                               ns_in,
                               nplane_in,
                               nproc_in,
                               gamma_only_in,
                               xprime_in);
        }
        if (double_flag)
        {
            fft_double->initfft(nx_in,
                                ny_in,
                                nz_in,
                                lixy_in,
                                rixy_in,
                                ns_in,
                                nplane_in,
                                nproc_in,
                                gamma_only_in,
                                xprime_in);
        }
    }
    if (device=="gpu")
    {
        #if defined(__ROCM)
            fft_float = make_unique<FFT_ROCM<float>>();
            fft_float->initfft(nx_in,ny_in,nz_in);
            fft_double = make_unique<FFT_ROCM<double>>();
            fft_double->initfft(nx_in,ny_in,nz_in);
        #elif defined(__CUDA)
            fft_float = make_unique<FFT_CUDA<float>>();
            fft_float->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in,
                               nplane_in, nproc_in, gamma_only_in, xprime_in);
            fft_double = make_unique<FFT_CUDA<double>>();
            fft_double->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in,
                                nplane_in, nproc_in, gamma_only_in, xprime_in);
        #endif
        // Host-memory interface: always provide a double CPU backend, and a
        // float one when float FFTs are enabled. Single precision still needs
        // a double host backend for host-only post-processing.
        if (float_flag)
        {
            host_fft_float = make_unique<FFT_CPU<float>>(this->fft_mode);
            host_fft_float->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in,
                                    nplane_in, nproc_in, gamma_only_in, xprime_in);
        }
        host_fft_double = make_unique<FFT_CPU<double>>(this->fft_mode);
        host_fft_double->initfft(nx_in, ny_in, nz_in, lixy_in, rixy_in, ns_in,
                                 nplane_in, nproc_in, gamma_only_in, xprime_in);
    }

}

void FFT_Bundle::setupFFT()
{
    if (double_flag){fft_double->setupFFT();}
    if (float_flag) {fft_float->setupFFT();}
    if (host_fft_double != nullptr) { host_fft_double->setupFFT(); }
    if (host_fft_float != nullptr) { host_fft_float->setupFFT(); }
}

void FFT_Bundle::clearFFT()
{
    if (double_flag){fft_double->cleanFFT();}
    if (float_flag) {fft_float->cleanFFT();}
    if (host_fft_double != nullptr) { host_fft_double->cleanFFT(); }
    if (host_fft_float != nullptr) { host_fft_float->cleanFFT(); }
}
void FFT_Bundle::clear()
{
    this->clearFFT();
    if (double_flag){fft_double->clear();}
    if (float_flag) {fft_float->clear();}
    if (host_fft_double != nullptr) { host_fft_double->clear(); }
    if (host_fft_float != nullptr) { host_fft_float->clear(); }
}

template <> void 
FFT_Bundle::fftxyfor(std::complex<float>* in, 
                     std::complex<float>* out) 
const {get_host_fft_float()->fftxyfor(in,out);}
template <> void 
FFT_Bundle::fftxyfor(std::complex<double>* in,
                     std::complex<double>* out) 
const {get_host_fft_double()->fftxyfor(in,out);}


template <> void 
FFT_Bundle::fftzfor(std::complex<float>* in, 
                    std::complex<float>* out) 
const {get_host_fft_float()->fftzfor(in,out);}
template <> void 
FFT_Bundle::fftzfor(std::complex<double>* in, 
                    std::complex<double>* out) 
const {get_host_fft_double()->fftzfor(in,out);}

template <> void 
FFT_Bundle::fftxybac(std::complex<float>* in, 
                     std::complex<float>* out) 
const {get_host_fft_float()->fftxybac(in,out);}
template <> void 
FFT_Bundle::fftxybac(std::complex<double>* in, 
                     std::complex<double>* out) 
const {get_host_fft_double()->fftxybac(in,out);}

template <> void 
FFT_Bundle::fftzbac(std::complex<float>* in, 
                    std::complex<float>* out) 
const {get_host_fft_float()->fftzbac(in,out);}
template <> void 
FFT_Bundle::fftzbac(std::complex<double>* in, 
                    std::complex<double>* out) 
const {get_host_fft_double()->fftzbac(in,out);}

template <> void 
FFT_Bundle::fftxyr2c(float* in, 
                     std::complex<float>* out) 
const {get_host_fft_float()->fftxyr2c(in,out);}
template <> void 
FFT_Bundle::fftxyr2c(double* in, 
                    std::complex<double>* out) 
const {get_host_fft_double()->fftxyr2c(in,out);}

template <> void 
FFT_Bundle::fftxyc2r(std::complex<float>* in, 
                    float* out) 
const {get_host_fft_float()->fftxyc2r(in,out);}
template <> void 
FFT_Bundle::fftxyc2r(std::complex<double>* in, 
                    double* out) 
const {get_host_fft_double()->fftxyc2r(in,out);}

template <> void  
FFT_Bundle::fft3D_forward(const base_device::DEVICE_GPU* ctx, 
                                std::complex<float>* in, 
                                std::complex<float>* out)
const {fft_float->fft3D_forward(in, out);}
template <> void  
FFT_Bundle::fft3D_forward(const base_device::DEVICE_GPU* ctx, 
                            std::complex<double>* in, 
                            std::complex<double>* out) 
const {fft_double->fft3D_forward(in, out);}

template <> void  
FFT_Bundle::fft3D_backward(const base_device::DEVICE_GPU* ctx, 
                                std::complex<float>* in, 
                                std::complex<float>* out) 
const {fft_float->fft3D_backward(in, out);}
template <> void  
FFT_Bundle::fft3D_backward(const base_device::DEVICE_GPU* ctx, 
                                 std::complex<double>* in, 
                                 std::complex<double>* out) 
const {fft_double->fft3D_backward(in, out);}

#if defined(__CUFFTDX)
int FFT_Bundle::fused_veff_grid_size() const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->fused_veff_grid_size();
    }
#endif
    return 0;
}

int FFT_Bundle::fused_veff_stick_size() const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->fused_veff_stick_size();
    }
#endif
    return 0;
}

int FFT_Bundle::spinor_batch_size() const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->spinor_batch_size();
    }
#endif
    return 1;
}

std::complex<double>* FFT_Bundle::get_spinor_stick_data(int bands) const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->get_spinor_stick_data(bands);
    }
#endif
    return nullptr;
}

void FFT_Bundle::transform_spinor_sticks(std::complex<double>* up, std::complex<double>* down,
                                         const int* reverse, bool inverse, int bands) const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        cuda->transform_spinor_sticks(up, down, reverse, inverse, bands);
    }
#endif
}

void FFT_Bundle::apply_spinor_potential(const double* potential, std::complex<double>* up,
                                        std::complex<double>* down, const int* reverse, int bands) const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        cuda->apply_spinor_potential(potential, up, down, reverse, bands);
    }
#endif
}

int FFT_Bundle::lixy() const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->lixy();
    }
#endif
    return 0;
}

int FFT_Bundle::rixy() const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        return cuda->rixy();
    }
#endif
    return 0;
}
#endif

void FFT_Bundle::fft3D_forward_pw(std::complex<double>* in, std::complex<double>* out) const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        cuda->fft3D_forward_pw(in, out);
        return;
    }
#endif
    fft_double->fft3D_forward(in, out);
}

void FFT_Bundle::fft3D_backward_pw(std::complex<double>* in, std::complex<double>* out) const
{
#if defined(__CUDA)
    auto* cuda = dynamic_cast<FFT_CUDA<double>*>(fft_double.get());
    if (cuda != nullptr)
    {
        cuda->fft3D_backward_pw(in, out);
        return;
    }
#endif
    fft_double->fft3D_backward(in, out);
}

// access the real space data
template <> float* 
FFT_Bundle::get_rspace_data()  const {return get_host_fft_float()->get_rspace_data();}
template <> double* 
FFT_Bundle::get_rspace_data()  const {return get_host_fft_double()->get_rspace_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxr_data()    const {return get_host_fft_float()->get_auxr_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxr_data()    const {return get_host_fft_double()->get_auxr_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxg_data()    const {return get_host_fft_float()->get_auxg_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxg_data()    const {return get_host_fft_double()->get_auxg_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxr_3d_data() const {return fft_float->get_auxr_3d_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxr_3d_data() const {return fft_double->get_auxr_3d_data();}
}
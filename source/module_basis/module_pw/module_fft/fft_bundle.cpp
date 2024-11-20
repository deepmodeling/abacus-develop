#include <cassert>
#include "fft_bundle.h"
#include "fft_cpu.h"
#include "module_base/module_device/device.h"
#if defined(__CUDA)
#include "fft_cuda.h"
#endif
#if defined(__ROCM)
#include "fft_rcom.h"
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
            fft_float = new FFT_RCOM<float>();
            fft_float->initfft(nx_in,ny_in,nz_in);
            fft_double = new FFT_RCOM<double>();
            fft_double->initfft(nx_in,ny_in,nz_in);
        #elif defined(__CUDA)
            fft_float = make_unique<FFT_CUDA<float>>();
            fft_float->initfft(nx_in,ny_in,nz_in);
            fft_double = make_unique<FFT_CUDA<double>>();
            fft_double->initfft(nx_in,ny_in,nz_in);
        #endif
    }

}

void FFT_Bundle::setupFFT()
{
    if (double_flag){fft_double->setupFFT();}
    if (float_flag) {fft_float->setupFFT();}
}

void FFT_Bundle::clearFFT()
{
    if (double_flag){fft_double->cleanFFT();}
    if (float_flag) {fft_float->cleanFFT();}
}
void FFT_Bundle::clear()
{
    this->clearFFT();
    if (double_flag){fft_double->clear();}
    if (float_flag) {fft_float->clear();}
}

template <> void 
FFT_Bundle::fftxyfor(std::complex<float>* in, 
                     std::complex<float>* out) 
const {fft_float->fftxyfor(in,out);}
template <> void 
FFT_Bundle::fftxyfor(std::complex<double>* in,
                     std::complex<double>* out) 
const {fft_double->fftxyfor(in,out);}


template <> void 
FFT_Bundle::fftzfor(std::complex<float>* in, 
                    std::complex<float>* out) 
const {fft_float->fftzfor(in,out);}
template <> void 
FFT_Bundle::fftzfor(std::complex<double>* in, 
                    std::complex<double>* out) 
const {fft_double->fftzfor(in,out);}

template <> void 
FFT_Bundle::fftxybac(std::complex<float>* in, 
                     std::complex<float>* out) 
const {fft_float->fftxybac(in,out);}
template <> void 
FFT_Bundle::fftxybac(std::complex<double>* in, 
                     std::complex<double>* out) 
const {fft_double->fftxybac(in,out);}

template <> void 
FFT_Bundle::fftzbac(std::complex<float>* in, 
                    std::complex<float>* out) 
const {fft_float->fftzbac(in,out);}
template <> void 
FFT_Bundle::fftzbac(std::complex<double>* in, 
                    std::complex<double>* out) 
const {fft_double->fftzbac(in,out);}

template <> void 
FFT_Bundle::fftxyr2c(float* in, 
                     std::complex<float>* out) 
const {fft_float->fftxyr2c(in,out);}
template <> void 
FFT_Bundle::fftxyr2c(double* in, 
                    std::complex<double>* out) 
const {fft_double->fftxyr2c(in,out);}

template <> void 
FFT_Bundle::fftxyc2r(std::complex<float>* in, 
                    float* out) 
const {fft_float->fftxyc2r(in,out);}
template <> void 
FFT_Bundle::fftxyc2r(std::complex<double>* in, 
                    double* out) 
const {fft_double->fftxyc2r(in,out);}

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

// access the real space data
template <> float* 
FFT_Bundle::get_rspace_data()  const {return fft_float->get_rspace_data();}
template <> double* 
FFT_Bundle::get_rspace_data()  const {return fft_double->get_rspace_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxr_data()    const {return fft_float->get_auxr_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxr_data()    const {return fft_double->get_auxr_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxg_data()    const {return fft_float->get_auxg_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxg_data()    const {return fft_double->get_auxg_data();}

template <> std::complex<float>* 
FFT_Bundle::get_auxr_3d_data() const {return fft_float->get_auxr_3d_data();}
template <> std::complex<double>* 
FFT_Bundle::get_auxr_3d_data() const {return fft_double->get_auxr_3d_data();}
}
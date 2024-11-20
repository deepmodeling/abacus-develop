#include "fft_cuda.h"
#include "module_base/module_device/memory_op.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
namespace ModulePW
{
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in, 
                               int ny_in, 
                               int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
}
template <>
void FFT_CUDA<float>::setupFFT()
{
    cufftPlan3d(&c_handle, this->nx, this->ny, this->nz, CUFFT_C2C);
    resmem_cd_op()(gpu_ctx, this->c_auxr_3d, this->nx * this->ny * this->nz);
        
}
template <>  
void FFT_CUDA<double>::setupFFT()
{
    cufftPlan3d(&z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z);
    resmem_zd_op()(gpu_ctx, this->z_auxr_3d, this->nx * this->ny * this->nz);
}
template <>
void FFT_CUDA<float>::cleanFFT()
{
    if (c_handle)
    {
        cufftDestroy(c_handle);
        c_handle = {};
    }
}
template <>
void FFT_CUDA<double>::cleanFFT()
{
    if (z_handle)
    {
        cufftDestroy(z_handle);
        z_handle = {};
    }
}
template <>
void FFT_CUDA<float>::clear()
{
    this->cleanFFT();
    if (c_auxr_3d != nullptr)
    {
        delmem_cd_op()(gpu_ctx, c_auxr_3d);
        c_auxr_3d = nullptr;
    }
}
template <>
void FFT_CUDA<double>::clear()
{
    this->cleanFFT();
    if (z_auxr_3d != nullptr)
    {
        delmem_zd_op()(gpu_ctx, z_auxr_3d);
        z_auxr_3d = nullptr;
    }
}

template <>
void FFT_CUDA<float>::fft3D_forward(std::complex<float>* in, 
                                    std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle, 
                             reinterpret_cast<cufftComplex*>(in), 
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<double>::fft3D_forward(std::complex<double>* in, 
                                     std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle, 
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out), 
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<float>::fft3D_backward(std::complex<float>* in, 
                                     std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle, 
                             reinterpret_cast<cufftComplex*>(in), 
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_INVERSE));
}

template <>
void FFT_CUDA<double>::fft3D_backward(std::complex<double>* in, 
                                      std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle, 
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out), 
                             CUFFT_INVERSE));
}
template <> std::complex<float>* 
FFT_CUDA<float>::get_auxr_3d_data()  const {return this->c_auxr_3d;}
template <> std::complex<double>* 
FFT_CUDA<double>::get_auxr_3d_data() const {return this->z_auxr_3d;}
}// namespace ModulePW
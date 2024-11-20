#include "fft_rcom.h"
#include "module_base/module_device/memory_op.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
namespace ModulePW
{
template <typename FPTYPE>
void FFT_RCOM<FPTYPE>::initfft(int nx_in, 
                               int ny_in, 
                               int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
}
template <>
void FFT_RCOM<float>::setupFFT()
{
    hipfftPlan3d(&c_handle, this->nx, this->ny, this->nz, HIPFFT_C2C);
    resmem_cd_op()(gpu_ctx, this->c_auxr_3d, this->nx * this->ny * this->nz);
        
}
template <>  
void FFT_RCOM<double>::setupFFT()
{
    hipfftPlan3d(&z_handle, this->nx, this->ny, this->nz, HIPFFT_Z2Z);
    resmem_zd_op()(gpu_ctx, this->z_auxr_3d, this->nx * this->ny * this->nz);
}
template <>
void FFT_RCOM<float>::cleanFFT()
{
    if (c_handle)
    {
        hipfftDestroy(c_handle);
        c_handle = {};
    }
}
template <>
void FFT_RCOM<double>::cleanFFT()
{
    if (z_handle)
    {
        hipfftDestroy(z_handle);
        z_handle = {};
    }
}
template <>
void FFT_RCOM<float>::clear()
{
    this->cleanFFT();
    if (c_auxr_3d != nullptr)
    {
        delmem_cd_op()(gpu_ctx, c_auxr_3d);
        c_auxr_3d = nullptr;
    }
}
template <>
void FFT_RCOM<double>::clear()
{
    this->cleanFFT();
    if (z_auxr_3d != nullptr)
    {
        delmem_zd_op()(gpu_ctx, z_auxr_3d);
        z_auxr_3d = nullptr;
    }
}
template <>
void FFT_RCOM<float>::fft3D_forward(std::complex<float>* in, 
                                    std::complex<float>* out) const
{
    CHECK_CUFFT(hipfftExecC2C(this->c_handle, 
                              reinterpret_cast<hipfftComplex*>(in),
                              reinterpret_cast<hipfftComplex*>(out), 
                              HIPFFT_FORWARD));
}
template <>
void FFT_RCOM<double>::fft3D_forward(std::complex<double>* in, 
                                     std::complex<double>* out) const
{
    CHECK_CUFFT(hipfftExecZ2Z(this->z_handle, 
                              reinterpret_cast<hipfftDoubleComplex*>(in),
                              reinterpret_cast<hipfftDoubleComplex*>(out), 
                              HIPFFT_FORWARD));
}
template <>
void FFT_RCOM<float>::fft3D_backward(std::complex<float>* in, 
                                     std::complex<float>* out) const
{
    CHECK_CUFFT(hipfftExecC2C(this->c_handle, 
                              reinterpret_cast<hipfftComplex*>(in),
                              reinterpret_cast<hipfftComplex*>(out), 
                              HIPFFT_BACKWARD));
}
template <>
void FFT_RCOM<double>::fft3D_backward(std::complex<double>* in, 
                                      std::complex<double>* out) const
{
    CHECK_CUFFT(hipfftExecZ2Z(this->z_handle, 
                              reinterpret_cast<hipfftDoubleComplex*>(in),
                              reinterpret_cast<hipfftDoubleComplex*>(out), 
                              HIPFFT_BACKWARD));
}
template <> std::complex<float>* 
FFT_RCOM<float>::get_auxr_3d_data()  const {return this->c_auxr_3d;}
template <> std::complex<double>* 
FFT_RCOM<double>::get_auxr_3d_data() const {return this->z_auxr_3d;}
}// namespace ModulePW
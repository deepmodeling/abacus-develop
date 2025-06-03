#include "fft_cuda_batch.h"

#include "module_base/module_device/memory_op.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModulePW
{
template <typename FPTYPE>
void FFT_CUDA_BATCH<FPTYPE>::initfft(int nx_in, int ny_in, int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->batch = 10;
}
template <>
void FFT_CUDA_BATCH<float>::setupFFT()
{
    int rank = 3;                  
    int n[3] = {this->nx, this->ny, this->nz};     
    const int size = this->nx* this->ny *this->nz; 
    cufftPlanMany(&c_handle, rank, n,
                  n, 1, size, 
                  n, 1, size,
                  CUFFT_Z2Z, this->batch);      
    resmem_cd_op()(this->c_auxr_3d, this->nx * this->ny * this->nz * this->batch);
}
template <>
void FFT_CUDA_BATCH<double>::setupFFT()
{
    std::cout<<"the nx ,ny,nz,batch is: "
             <<this->nx<<" "<<this->ny<<" "<<this->nz<<" "<<this->batch<<std::endl;
    int rank = 3;                
    int n[3] = {this->nx, this->ny, this->nz};       
    const int size = this->nx* this->ny *this->nz; 
    cufftPlanMany(&z_handle, rank, n,
                  n, 1, size, 
                  n, 1, size, 
                  CUFFT_Z2Z, this->batch);       
    cudaMalloc(&z_auxr_3d, this->nx * this->ny * this->nz * this->batch * sizeof(std::complex<double>));
    // resmem_zd_op()(this->z_auxr_3d, this->nx * this->ny * this->nz * this->batch);
}
template <>
void FFT_CUDA_BATCH<float>::cleanFFT()
{
    if (c_handle)
    {
        cufftDestroy(c_handle);
        c_handle = {};
    }
}
template <>
void FFT_CUDA_BATCH<double>::cleanFFT()
{
    if (z_handle)
    {
        cufftDestroy(z_handle);
        z_handle = {};
    }
}
template <>
void FFT_CUDA_BATCH<float>::clear()
{
    this->cleanFFT();
    if (c_auxr_3d != nullptr)
    {
        delmem_cd_op()(c_auxr_3d);
        c_auxr_3d = nullptr;
    }
}
template <>
void FFT_CUDA_BATCH<double>::clear()
{
    this->cleanFFT();
    if (z_auxr_3d != nullptr)
    {
        delmem_zd_op()(z_auxr_3d);
        z_auxr_3d = nullptr;
    }
}

template <>
void FFT_CUDA_BATCH<float>::fft3D_forward(std::complex<float>* in, std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA_BATCH<double>::fft3D_forward(std::complex<double>* in, std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA_BATCH<float>::fft3D_backward(std::complex<float>* in, std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_INVERSE));
}

template <>
void FFT_CUDA_BATCH<double>::fft3D_backward(std::complex<double>* in, std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_INVERSE));
}
template <>
std::complex<float>* FFT_CUDA_BATCH<float>::get_auxr_3d_data() const
{
    return this->c_auxr_3d;
}
template <>
std::complex<double>* FFT_CUDA_BATCH<double>::get_auxr_3d_data() const
{
    return this->z_auxr_3d;
}

template FFT_CUDA_BATCH<float>::FFT_CUDA_BATCH();
template FFT_CUDA_BATCH<float>::~FFT_CUDA_BATCH();
template FFT_CUDA_BATCH<double>::FFT_CUDA_BATCH();
template FFT_CUDA_BATCH<double>::~FFT_CUDA_BATCH();
} // namespace ModulePW
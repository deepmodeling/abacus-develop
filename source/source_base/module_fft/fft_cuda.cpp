#include "fft_cuda.h"

#include "source_base/module_device/memory_op.h"
#include "source_pw/module_pwdft/global.h"

namespace ModuleBase
{
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in, int ny_in, int nz_in, int batch_size)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->batch_size = batch_size;
}
template <>
void FFT_CUDA<float>::setupFFT()
{
    if (this->batch_size){
        int rank = 3; // this means the dimension is 3
        int n[3] = {this->nx, this->ny, this->nz};
        int inembed[3] = {this->nx, this->ny, this->nz};
        int onembed[3] = {this->nx, this->ny, this->nz};
        int istride = 1, ostride = 1;
        size_t N = static_cast<size_t>(this->nx) * this->ny * this->nz;
        int idist = N;
        int odist = N;
        cufftPlanMany(&c_handle, rank, n,
                      inembed, istride, idist,
                      onembed, ostride, odist,
                      CUFFT_C2C, this->batch_size)
    }
    else{
        cufftPlan3d(&c_handle, this->nx, this->ny, this->nz, CUFFT_C2C);
        resmem_cd_op()(this->c_auxr_3d, this->nx * this->ny * this->nz);
    }
    
}
template <>
void FFT_CUDA<double>::setupFFT()
{
    if (this->batch_size){
        int rank = 3; // this means the dimension is 3
        int n[3] = {this->nx, this->ny, this->nz};
        int inembed[3] = {this->nx, this->ny, this->nz};
        int onembed[3] = {this->nx, this->ny, this->nz};
        int istride = 1, ostride = 1;
        size_t N = static_cast<size_t>(this->nx) * this->ny * this->nz;
        int idist = N;
        int odist = N;
        cufftPlanMany(&z_handle, rank, n,
                      inembed, istride, idist,
                      onembed, ostride, odist,
                      CUFFT_Z2Z, this->batch_size)
    }
    else{
        cufftPlan3d(&z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z);
        resmem_zd_op()(this->z_auxr_3d, this->nx * this->ny * this->nz);
    }
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
        delmem_cd_op()(c_auxr_3d);
        c_auxr_3d = nullptr;
    }
}
template <>
void FFT_CUDA<double>::clear()
{
    this->cleanFFT();
    if (z_auxr_3d != nullptr)
    {
        delmem_zd_op()(z_auxr_3d);
        z_auxr_3d = nullptr;
    }
}

template <>
void FFT_CUDA<float>::fft3D_forward(std::complex<float>* in, std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<double>::fft3D_forward(std::complex<double>* in, std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<float>::fft3D_backward(std::complex<float>* in, std::complex<float>* out) const
{
    CHECK_CUFFT(cufftExecC2C(this->c_handle,
                             reinterpret_cast<cufftComplex*>(in),
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_INVERSE));
}

template <>
void FFT_CUDA<double>::fft3D_backward(std::complex<double>* in, std::complex<double>* out) const
{
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle,
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out),
                             CUFFT_INVERSE));
}
template <>
std::complex<float>* FFT_CUDA<float>::get_auxr_3d_data() const
{
    return this->c_auxr_3d;
}
template <>
std::complex<double>* FFT_CUDA<double>::get_auxr_3d_data() const
{
    return this->z_auxr_3d;
}

template FFT_CUDA<float>::FFT_CUDA();
template FFT_CUDA<float>::~FFT_CUDA();
template FFT_CUDA<double>::FFT_CUDA();
template FFT_CUDA<double>::~FFT_CUDA();
} // namespace ModuleBase
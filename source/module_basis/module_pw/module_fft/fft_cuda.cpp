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
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in,
                               int ny_in,
                               int nz_in,
                               int lixy_in,
                               int rixy_in,
                               int /*ns_in*/,
                               int /*nplane_in*/,
                               int nproc_in,
                               bool gamma_only_in,
                               bool xprime_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->lixy = lixy_in;
    this->rixy = rixy_in;
    this->nproc = nproc_in;
    this->gamma_only = gamma_only_in;
    this->xprime = xprime_in;
}
template <>
void FFT_CUDA<float>::setupFFT()
{
}
template <>  
void FFT_CUDA<double>::setupFFT()
{
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
    for (auto* handle : {&z_x_handle, &z_yz_left_handle, &z_yz_right_handle})
    {
        if (*handle)
        {
            cufftDestroy(*handle);
            *handle = {};
        }
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
    if (!this->c_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->c_handle, this->nx, this->ny, this->nz, CUFFT_C2C));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_handle, 
                             reinterpret_cast<cufftComplex*>(in), 
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<double>::fft3D_forward(std::complex<double>* in, 
                                     std::complex<double>* out) const
{
    if (!this->z_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle, 
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out), 
                             CUFFT_FORWARD));
}
template <>
void FFT_CUDA<float>::fft3D_backward(std::complex<float>* in, 
                                     std::complex<float>* out) const
{
    if (!this->c_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->c_handle, this->nx, this->ny, this->nz, CUFFT_C2C));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_handle, 
                             reinterpret_cast<cufftComplex*>(in), 
                             reinterpret_cast<cufftComplex*>(out),
                             CUFFT_INVERSE));
}

template <>
void FFT_CUDA<double>::fft3D_backward(std::complex<double>* in, 
                                      std::complex<double>* out) const
{
    if (!this->z_handle)
    {
        CHECK_CUFFT(cufftPlan3d(&this->z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_handle, 
                             reinterpret_cast<cufftDoubleComplex*>(in),
                             reinterpret_cast<cufftDoubleComplex*>(out), 
                             CUFFT_INVERSE));
}
template <typename FPTYPE>
bool FFT_CUDA<FPTYPE>::use_pruned_gpu_fft() const
{
    return false;
}
template <>
bool FFT_CUDA<double>::use_pruned_gpu_fft() const
{
    // Small grids are dominated by the extra kernel launches. Sparse large
    // grids avoid the y-z transforms in the unoccupied middle x slabs.
    return this->xprime && !this->gamma_only && this->nproc == 1
           && this->nx * this->ny * this->nz >= 1048576
           && 4 * (this->lixy + 1 + this->nx - this->rixy) <= 3 * this->nx;
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::init_gpu_pruned_plan() const
{
}
template <>
void FFT_CUDA<double>::init_gpu_pruned_plan() const
{
    if (this->z_x_handle)
    {
        return;
    }
    int x_size = this->nx;
    int yz_size[2] = {this->ny, this->nz};
    const int plane = this->ny * this->nz;
    CHECK_CUFFT(cufftPlanMany(&this->z_x_handle, 1, &x_size, &x_size, plane, 1,
                             &x_size, plane, 1, CUFFT_Z2Z, plane));
    CHECK_CUFFT(cufftPlanMany(&this->z_yz_left_handle, 2, yz_size, yz_size, 1, plane,
                             yz_size, 1, plane, CUFFT_Z2Z, this->lixy + 1));
    if (this->rixy < this->nx)
    {
        CHECK_CUFFT(cufftPlanMany(&this->z_yz_right_handle, 2, yz_size, yz_size, 1, plane,
                                 yz_size, 1, plane, CUFFT_Z2Z, this->nx - this->rixy));
    }
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::fft3D_forward_pw(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const
{
    this->fft3D_forward(in, out);
}
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::fft3D_backward_pw(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const
{
    this->fft3D_backward(in, out);
}
template <>
void FFT_CUDA<double>::fft3D_forward_pw(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->use_pruned_gpu_fft())
    {
        this->fft3D_forward(in, out);
        return;
    }
    this->init_gpu_pruned_plan();
    auto* input = reinterpret_cast<cufftDoubleComplex*>(in);
    auto* output = reinterpret_cast<cufftDoubleComplex*>(out);
    CHECK_CUFFT(cufftExecZ2Z(this->z_x_handle, input, output, CUFFT_FORWARD));
    CHECK_CUFFT(cufftExecZ2Z(this->z_yz_left_handle, output, output, CUFFT_FORWARD));
    if (this->z_yz_right_handle)
    {
        auto* right = output + this->rixy * this->ny * this->nz;
        CHECK_CUFFT(cufftExecZ2Z(this->z_yz_right_handle, right, right, CUFFT_FORWARD));
    }
}
template <>
void FFT_CUDA<double>::fft3D_backward_pw(std::complex<double>* in, std::complex<double>* out) const
{
    if (!this->use_pruned_gpu_fft())
    {
        this->fft3D_backward(in, out);
        return;
    }
    this->init_gpu_pruned_plan();
    auto* input = reinterpret_cast<cufftDoubleComplex*>(in);
    auto* output = reinterpret_cast<cufftDoubleComplex*>(out);
    CHECK_CUFFT(cufftExecZ2Z(this->z_yz_left_handle, input, input, CUFFT_INVERSE));
    if (this->z_yz_right_handle)
    {
        auto* right = input + this->rixy * this->ny * this->nz;
        CHECK_CUFFT(cufftExecZ2Z(this->z_yz_right_handle, right, right, CUFFT_INVERSE));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_x_handle, input, output, CUFFT_INVERSE));
}

template <> std::complex<float>* 
FFT_CUDA<float>::get_auxr_3d_data()  const
{
    if (this->c_auxr_3d == nullptr)
    {
        resmem_cd_op()(gpu_ctx, this->c_auxr_3d, this->nx * this->ny * this->nz);
    }
    return this->c_auxr_3d;
}
template <> std::complex<double>* 
FFT_CUDA<double>::get_auxr_3d_data() const
{
    if (this->z_auxr_3d == nullptr)
    {
        resmem_zd_op()(gpu_ctx, this->z_auxr_3d, this->nx * this->ny * this->nz);
    }
    return this->z_auxr_3d;
}

template void FFT_CUDA<float>::fft3D_forward_pw(std::complex<float>*, std::complex<float>*) const;
template void FFT_CUDA<float>::fft3D_backward_pw(std::complex<float>*, std::complex<float>*) const;
template FFT_CUDA<float>::FFT_CUDA();
template FFT_CUDA<float>::~FFT_CUDA();
template FFT_CUDA<double>::FFT_CUDA();
template FFT_CUDA<double>::~FFT_CUDA();
}// namespace ModulePW
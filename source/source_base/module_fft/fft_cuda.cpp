#include "fft_cuda.h"

#include <cassert>
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/device_check.h"

namespace ModuleBase
{
template <typename FPTYPE>
void FFT_CUDA<FPTYPE>::initfft(int nx_in, int ny_in, int nz_in)
{
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
}
template <>
void FFT_CUDA<float>::setupFFT()
{
    cufftPlan3d(&c_handle, this->nx, this->ny, this->nz, CUFFT_C2C);
    resmem_cd_op()(this->c_auxr_3d, this->nx * this->ny * this->nz);
}
template <>
void FFT_CUDA<double>::setupFFT()
{
    cufftPlan3d(&z_handle, this->nx, this->ny, this->nz, CUFFT_Z2Z);
    resmem_zd_op()(this->z_auxr_3d, this->nx * this->ny * this->nz);
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

// ============================================================================
// Batch FFT cleanup (must be defined before clear() which calls it)
// ============================================================================
template <>
void FFT_CUDA<float>::cleanBatchFFT()
{
    if (c_batch_handle != 0)
    {
        cufftDestroy(c_batch_handle);
        c_batch_handle = 0;
    }
    if (c_auxr_batch_in != nullptr)
    {
        delmem_cd_op()(c_auxr_batch_in);
        c_auxr_batch_in = nullptr;
    }
    if (c_auxr_batch_out != nullptr)
    {
        delmem_cd_op()(c_auxr_batch_out);
        c_auxr_batch_out = nullptr;
    }
}

template <>
void FFT_CUDA<double>::cleanBatchFFT()
{
    if (z_batch_handle != 0)
    {
        cufftDestroy(z_batch_handle);
        z_batch_handle = 0;
    }
    if (z_auxr_batch_in != nullptr)
    {
        delmem_zd_op()(z_auxr_batch_in);
        z_auxr_batch_in = nullptr;
    }
    if (z_auxr_batch_out != nullptr)
    {
        delmem_zd_op()(z_auxr_batch_out);
        z_auxr_batch_out = nullptr;
    }
}

template <>
void FFT_CUDA<float>::clear()
{
    this->cleanFFT();
    this->cleanBatchFFT();
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
    this->cleanBatchFFT();
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

// ============================================================================
// Batch FFT Implementation
// ============================================================================

template <>
void FFT_CUDA<float>::setupBatchFFT(int batch_size_in)
{
    if (c_batch_handle != 0 || batch_size_in <= 1)
    {
        // Already initialized
        return;
    }
    if (this->c_handle == 0)
    {
        this->setupFFT();
    }

    this->batch_size = batch_size_in;

    const int rank = 3;
    int n[3] = {this->nx, this->ny, this->nz};
    const int idist = this->nx * this->ny * this->nz;
    const int odist = this->nx * this->ny * this->nz;
    const int istride = 1;
    const int ostride = 1;

    // Create batch FFT plan using cufftPlanMany
    CHECK_CUFFT(cufftPlanMany(&c_batch_handle,
                              rank, n,
                              nullptr, istride, idist,  // input parameters
                              nullptr, ostride, odist,  // output parameters
                              CUFFT_C2C,
                              this->batch_size));

    // Allocate batch buffers on device
    const size_t batch_buffer_size = this->batch_size * this->nx * this->ny * this->nz;
    resmem_cd_op()(this->c_auxr_batch_in, batch_buffer_size, "FFT_CUDA::c_batch_in");
    resmem_cd_op()(this->c_auxr_batch_out, batch_buffer_size, "FFT_CUDA::c_batch_out");
}

template <>
void FFT_CUDA<double>::setupBatchFFT(int batch_size_in)
{
    if (z_batch_handle != 0 || batch_size_in <= 1)
    {
        // Already initialized
        return;
    }
    if (this->z_handle == 0)
    {
        this->setupFFT();
    }

    this->batch_size = batch_size_in;

    const int rank = 3;
    int n[3] = {this->nx, this->ny, this->nz};
    const int idist = this->nx * this->ny * this->nz;
    const int odist = this->nx * this->ny * this->nz;
    const int istride = 1;
    const int ostride = 1;

    // Create batch FFT plan using cufftPlanMany
    CHECK_CUFFT(cufftPlanMany(&z_batch_handle,
                              rank, n,
                              nullptr, istride, idist,  // input parameters
                              nullptr, ostride, odist,  // output parameters
                              CUFFT_Z2Z,
                              this->batch_size));

    // Allocate batch buffers on device
    const size_t batch_buffer_size = this->batch_size * this->nx * this->ny * this->nz;
    resmem_zd_op()(this->z_auxr_batch_in, batch_buffer_size, "FFT_CUDA::z_batch_in");
    resmem_zd_op()(this->z_auxr_batch_out, batch_buffer_size, "FFT_CUDA::z_batch_out");
}

template <>
void FFT_CUDA<float>::fft3D_forward_batch(std::complex<float>* in_batch,
                                          std::complex<float>* out_batch,
                                          int batch_count) const
{
    // Validate batch_count - programming error if out of range
    assert(batch_count > 0 && "batch_count must be positive");
    assert(batch_count <= this->batch_size && "batch_count exceeds allocated batch_size");
    const std::size_t nxyz = static_cast<std::size_t>(this->nx) * this->ny * this->nz;

    if (batch_count == this->batch_size)
    {
        CHECK_CUFFT(cufftExecC2C(this->c_batch_handle,
                                 reinterpret_cast<cufftComplex*>(in_batch),
                                 reinterpret_cast<cufftComplex*>(out_batch),
                                 CUFFT_FORWARD));
        return;
    }

    const std::size_t valid_size = static_cast<std::size_t>(batch_count) * nxyz;
    if (in_batch != this->c_auxr_batch_in)
    {
        CHECK_CUDA(cudaMemcpy(this->c_auxr_batch_in,
                              in_batch,
                              valid_size * sizeof(std::complex<float>),
                              cudaMemcpyDeviceToDevice));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_batch_handle,
                             reinterpret_cast<cufftComplex*>(this->c_auxr_batch_in),
                             reinterpret_cast<cufftComplex*>(this->c_auxr_batch_out),
                             CUFFT_FORWARD));
    if (out_batch != this->c_auxr_batch_out)
    {
        CHECK_CUDA(cudaMemcpy(out_batch,
                              this->c_auxr_batch_out,
                              valid_size * sizeof(std::complex<float>),
                              cudaMemcpyDeviceToDevice));
    }
}

template <>
void FFT_CUDA<double>::fft3D_forward_batch(std::complex<double>* in_batch,
                                           std::complex<double>* out_batch,
                                           int batch_count) const
{
    // Validate batch_count - programming error if out of range
    assert(batch_count > 0 && "batch_count must be positive");
    assert(batch_count <= this->batch_size && "batch_count exceeds allocated batch_size");
    const std::size_t nxyz = static_cast<std::size_t>(this->nx) * this->ny * this->nz;

    if (batch_count == this->batch_size)
    {
        CHECK_CUFFT(cufftExecZ2Z(this->z_batch_handle,
                                 reinterpret_cast<cufftDoubleComplex*>(in_batch),
                                 reinterpret_cast<cufftDoubleComplex*>(out_batch),
                                 CUFFT_FORWARD));
        return;
    }

    const std::size_t valid_size = static_cast<std::size_t>(batch_count) * nxyz;
    if (in_batch != this->z_auxr_batch_in)
    {
        CHECK_CUDA(cudaMemcpy(this->z_auxr_batch_in,
                              in_batch,
                              valid_size * sizeof(std::complex<double>),
                              cudaMemcpyDeviceToDevice));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_batch_handle,
                             reinterpret_cast<cufftDoubleComplex*>(this->z_auxr_batch_in),
                             reinterpret_cast<cufftDoubleComplex*>(this->z_auxr_batch_out),
                             CUFFT_FORWARD));
    if (out_batch != this->z_auxr_batch_out)
    {
        CHECK_CUDA(cudaMemcpy(out_batch,
                              this->z_auxr_batch_out,
                              valid_size * sizeof(std::complex<double>),
                              cudaMemcpyDeviceToDevice));
    }
}

template <>
void FFT_CUDA<float>::fft3D_backward_batch(std::complex<float>* in_batch,
                                           std::complex<float>* out_batch,
                                           int batch_count) const
{
    // Validate batch_count - programming error if out of range
    assert(batch_count > 0 && "batch_count must be positive");
    assert(batch_count <= this->batch_size && "batch_count exceeds allocated batch_size");
    const std::size_t nxyz = static_cast<std::size_t>(this->nx) * this->ny * this->nz;

    if (batch_count == this->batch_size)
    {
        CHECK_CUFFT(cufftExecC2C(this->c_batch_handle,
                                 reinterpret_cast<cufftComplex*>(in_batch),
                                 reinterpret_cast<cufftComplex*>(out_batch),
                                 CUFFT_INVERSE));
        return;
    }

    const std::size_t valid_size = static_cast<std::size_t>(batch_count) * nxyz;
    if (in_batch != this->c_auxr_batch_in)
    {
        CHECK_CUDA(cudaMemcpy(this->c_auxr_batch_in,
                              in_batch,
                              valid_size * sizeof(std::complex<float>),
                              cudaMemcpyDeviceToDevice));
    }
    CHECK_CUFFT(cufftExecC2C(this->c_batch_handle,
                             reinterpret_cast<cufftComplex*>(this->c_auxr_batch_in),
                             reinterpret_cast<cufftComplex*>(this->c_auxr_batch_out),
                             CUFFT_INVERSE));
    if (out_batch != this->c_auxr_batch_out)
    {
        CHECK_CUDA(cudaMemcpy(out_batch,
                              this->c_auxr_batch_out,
                              valid_size * sizeof(std::complex<float>),
                              cudaMemcpyDeviceToDevice));
    }
}

template <>
void FFT_CUDA<double>::fft3D_backward_batch(std::complex<double>* in_batch,
                                            std::complex<double>* out_batch,
                                            int batch_count) const
{
    // Validate batch_count - programming error if out of range
    assert(batch_count > 0 && "batch_count must be positive");
    assert(batch_count <= this->batch_size && "batch_count exceeds allocated batch_size");
    const std::size_t nxyz = static_cast<std::size_t>(this->nx) * this->ny * this->nz;

    if (batch_count == this->batch_size)
    {
        CHECK_CUFFT(cufftExecZ2Z(this->z_batch_handle,
                                 reinterpret_cast<cufftDoubleComplex*>(in_batch),
                                 reinterpret_cast<cufftDoubleComplex*>(out_batch),
                                 CUFFT_INVERSE));
        return;
    }

    const std::size_t valid_size = static_cast<std::size_t>(batch_count) * nxyz;
    if (in_batch != this->z_auxr_batch_in)
    {
        CHECK_CUDA(cudaMemcpy(this->z_auxr_batch_in,
                              in_batch,
                              valid_size * sizeof(std::complex<double>),
                              cudaMemcpyDeviceToDevice));
    }
    CHECK_CUFFT(cufftExecZ2Z(this->z_batch_handle,
                             reinterpret_cast<cufftDoubleComplex*>(this->z_auxr_batch_in),
                             reinterpret_cast<cufftDoubleComplex*>(this->z_auxr_batch_out),
                             CUFFT_INVERSE));
    if (out_batch != this->z_auxr_batch_out)
    {
        CHECK_CUDA(cudaMemcpy(out_batch,
                              this->z_auxr_batch_out,
                              valid_size * sizeof(std::complex<double>),
                              cudaMemcpyDeviceToDevice));
    }
}

template <>
bool FFT_CUDA<float>::is_batch_fft_ready() const
{
    return (c_batch_handle != 0 && c_auxr_batch_in != nullptr && c_auxr_batch_out != nullptr);
}

template <>
bool FFT_CUDA<double>::is_batch_fft_ready() const
{
    return (z_batch_handle != 0 && z_auxr_batch_in != nullptr && z_auxr_batch_out != nullptr);
}

template <>
int FFT_CUDA<float>::get_batch_size() const
{
    return this->batch_size;
}

template <>
int FFT_CUDA<double>::get_batch_size() const
{
    return this->batch_size;
}

template <>
std::complex<float>* FFT_CUDA<float>::get_batch_input_buffer() const
{
    return this->c_auxr_batch_in;
}

template <>
std::complex<double>* FFT_CUDA<double>::get_batch_input_buffer() const
{
    return this->z_auxr_batch_in;
}

template <>
std::complex<float>* FFT_CUDA<float>::get_batch_output_buffer() const
{
    return this->c_auxr_batch_out;
}

template <>
std::complex<double>* FFT_CUDA<double>::get_batch_output_buffer() const
{
    return this->z_auxr_batch_out;
}

// Template instantiations (must be at the end after all specializations)
template FFT_CUDA<float>::FFT_CUDA();
template FFT_CUDA<float>::~FFT_CUDA();
template FFT_CUDA<double>::FFT_CUDA();
template FFT_CUDA<double>::~FFT_CUDA();
} // namespace ModuleBase

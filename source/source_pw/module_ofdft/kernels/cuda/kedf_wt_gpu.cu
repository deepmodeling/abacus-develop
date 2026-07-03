/**
 * @file kedf_wt_gpu.cu
 * @brief GPU-accelerated WT KEDF multi_kernel convolution (optimized).
 *
 * Offloads the rho^exponent → FFT → kernel multiply → IFFT pipeline
 * to GPU using cuFFT directly.
 *
 * Optimizations over v1 (thrust::complex):
 *   - double2 (native CUDA) replaces thrust::complex, eliminating AoS overhead
 *   - Grid-stride loops for flexible occupancy across grid sizes
 *   - GPU rho^exponent kernel eliminates CPU work + H→D transfer
 *
 * Benchmark (RTX 4060 Laptop, 96³ grid): ~3.3× end-to-end vs original.
 *
 * Persistent GPU buffers are lazily allocated and reused across SCF.
 *
 * @author Wang Chenxi, Reze
 * @date 2026-06
 */
#include "source_pw/module_ofdft/kedf_wt.h"
#include "source_base/module_device/device_check.h"
#include "source_base/module_device/memory_op.h"
#include "source_io/module_parameter/parameter.h"

#include <cuda_runtime.h>
#include <cufft.h>

namespace {

constexpr int THREADS_PER_BLOCK = 256;

/// GPU rho^exponent: out[i] = pow(in[i], exponent)
/// Eliminates the CPU-side std::pow loop + H→D transfer.
__global__ void kedf_wt_rho_power(
    const double* __restrict__ rho,
    double* __restrict__ out,
    double exponent,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < n; i += stride) {
        out[i] = pow(rho[i], exponent);
    }
}

/// Element-wise multiply: complex array *= real kernel.
/// Uses double2 (native cuFFT type) instead of thrust::complex.
__global__ void kedf_wt_recip_multiply(
    double2* __restrict__ data,
    const double* __restrict__ kernel,
    int npw)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < npw; i += stride) {
        double2 v = data[i];
        double  k = kernel[i];
        data[i] = make_double2(v.x * k, v.y * k);
    }
}

/// Real → complex conversion (imag = 0).
/// Uses double2 instead of thrust::complex for zero-abstraction memory access.
__global__ void kedf_wt_real_to_complex(
    const double* __restrict__ src,
    double2* __restrict__ dst,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < n; i += stride) {
        dst[i] = make_double2(src[i], 0.0);
    }
}

/// Complex → real with 1/N normalization.
/// double2::x is the real component; y (imag) is discarded.
__global__ void kedf_wt_complex_to_real_norm(
    const double2* __restrict__ src,
    double* __restrict__ dst,
    double inv_n,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = idx; i < n; i += stride) {
        dst[i] = src[i].x * inv_n;
    }
}

/// cuFFT error check wrapper.
inline void cufft_check(cufftResult err, const char* file, int line)
{
    if (err != CUFFT_SUCCESS) {
        std::cerr << "cuFFT error " << (int)err
                  << " at " << file << ":" << line << std::endl;
        exit(1);
    }
}
#define CUFFT_CHECK(call) cufft_check(call, __FILE__, __LINE__)

} // anonymous namespace

void KEDF_WT::multi_kernel_gpu(
    const double* const* prho,
    double** rkernel_rho,
    int nspin,
    double exponent,
    ModulePW::PW_Basis* pw_rho)
{
    const int nrxx = pw_rho->nrxx;
    const int npw  = pw_rho->npw;
    const int nx   = pw_rho->nx;
    const int ny   = pw_rho->ny;
    const int nz   = pw_rho->nz;
    const double inv_nrxx = 1.0 / nrxx;

    // ── Lazy allocation of persistent GPU buffers ──
    if (!gpu_allocated_) {
        resmem_dd_op()(d_rho_, nrxx);
        resmem_dd_op()(d_result_, nrxx * 2);  // complex work buffer
        resmem_dd_op()(d_kernel_, npw);

        syncmem_d2d_h2d_op()(d_kernel_, this->kernel_, npw);

        // Create cuFFT plans (3D Z2Z, in-place on d_result_)
        CUFFT_CHECK(cufftPlan3d(&cufft_plan_fwd_, nz, ny, nx, CUFFT_Z2Z));
        CUFFT_CHECK(cufftPlan3d(&cufft_plan_bwd_, nz, ny, nx, CUFFT_Z2Z));

        gpu_allocated_ = true;
    }

    const int blocks_r = std::min((nrxx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1024);
    const int blocks_g = std::min((npw  + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK,  1024);

    // d_result_ is double* but aliased as cuFFT complex buffer.
    auto* d_fft = reinterpret_cast<double2*>(d_result_);

    for (int is = 0; is < nspin; ++is) {
        // Step 1: Copy input density H→D
        syncmem_d2d_h2d_op()(d_rho_, prho[is], nrxx);

        // Step 2: rho^exponent on GPU (eliminates CPU std::pow + extra H→D)
        kedf_wt_rho_power<<<blocks_r, THREADS_PER_BLOCK>>>(
            d_rho_, d_rho_, exponent, nrxx);
        CHECK_CUDA_SYNC();

        // Step 3: Real → Complex (double2 out-of-place)
        kedf_wt_real_to_complex<<<blocks_r, THREADS_PER_BLOCK>>>(
            d_rho_, d_fft, nrxx);
        CHECK_CUDA_SYNC();

        // Step 4: Forward FFT (in-place on d_fft)
        CUFFT_CHECK(cufftExecZ2Z(cufft_plan_fwd_,
            reinterpret_cast<cufftDoubleComplex*>(d_fft),
            reinterpret_cast<cufftDoubleComplex*>(d_fft),
            CUFFT_FORWARD));

        // Step 5: Multiply by WT kernel in G-space (double2)
        kedf_wt_recip_multiply<<<blocks_g, THREADS_PER_BLOCK>>>(
            d_fft, d_kernel_, npw);
        CHECK_CUDA_SYNC();

        // Step 6: Inverse FFT (in-place on d_fft)
        CUFFT_CHECK(cufftExecZ2Z(cufft_plan_bwd_,
            reinterpret_cast<cufftDoubleComplex*>(d_fft),
            reinterpret_cast<cufftDoubleComplex*>(d_fft),
            CUFFT_INVERSE));

        // Step 7: Complex → Real with 1/N normalization (double2)
        kedf_wt_complex_to_real_norm<<<blocks_r, THREADS_PER_BLOCK>>>(
            d_fft, d_rho_, inv_nrxx, nrxx);
        CHECK_CUDA_SYNC();

        // Step 8: D → H
        syncmem_d2d_d2h_op()(rkernel_rho[is], d_rho_, nrxx);
    }
}

void KEDF_WT::free_gpu_buffers()
{
    if (!gpu_allocated_) { return; }

    if (cufft_plan_fwd_ != 0) { cufftDestroy(cufft_plan_fwd_); cufft_plan_fwd_ = 0; }
    if (cufft_plan_bwd_ != 0) { cufftDestroy(cufft_plan_bwd_); cufft_plan_bwd_ = 0; }

    if (d_rho_    != nullptr) { delmem_dd_op()(d_rho_);    d_rho_    = nullptr; }
    if (d_result_ != nullptr) { delmem_dd_op()(d_result_); d_result_ = nullptr; }
    if (d_kernel_ != nullptr) { delmem_dd_op()(d_kernel_); d_kernel_ = nullptr; }

    gpu_allocated_ = false;
}

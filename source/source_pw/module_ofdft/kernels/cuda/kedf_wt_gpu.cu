/**
 * @file kedf_wt_gpu.cu
 * @brief GPU-accelerated WT KEDF multi_kernel convolution.
 *
 * Offloads the rho^exponent → FFT → kernel multiply → IFFT pipeline
 * to GPU using cuFFT directly (bypassing PW_Basis template API to
 * avoid cross-target template instantiation issues).
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
#include <thrust/complex.h>

namespace {

constexpr int THREADS_PER_BLOCK = 256;

/// Element-wise multiply: complex array *= real kernel.
__global__ void kedf_wt_recip_multiply(
    thrust::complex<double>* __restrict__ data,
    const double* __restrict__ kernel,
    int npw)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= npw) { return; }
    data[idx] *= kernel[idx];
}

/// Real → complex conversion (imag = 0).
__global__ void kedf_wt_real_to_complex(
    const double* __restrict__ src,
    thrust::complex<double>* __restrict__ dst,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) { return; }
    dst[idx] = thrust::complex<double>(src[idx], 0.0);
}

/// Complex → real with 1/N normalization (cuFFT inverse doesn't normalize).
__global__ void kedf_wt_complex_to_real_norm(
    const thrust::complex<double>* __restrict__ src,
    double* __restrict__ dst,
    double inv_n,
    int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) { return; }
    dst[idx] = src[idx].real() * inv_n;
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
        resmem_dd_op()(d_result_, nrxx * 2);  // reused as complex work buffer
        resmem_dd_op()(d_kernel_, npw);

        syncmem_d2d_h2d_op()(d_kernel_, this->kernel_, npw);

        // Create cuFFT plans (3D Z2Z, in-place)
        CUFFT_CHECK(cufftPlan3d(&cufft_plan_fwd_, nz, ny, nx, CUFFT_Z2Z));
        CUFFT_CHECK(cufftPlan3d(&cufft_plan_bwd_, nz, ny, nx, CUFFT_Z2Z));

        gpu_allocated_ = true;
    }

    const int blocks = (npw + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    // d_result_ is double* but reused as cuFFT complex buffer (nrxx*2 = nrxx complex doubles).
    auto* d_fft = reinterpret_cast<cufftDoubleComplex*>(d_result_);

    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        // Step 1: rho^exponent on CPU
        for (int ir = 0; ir < nrxx; ++ir) {
            rkernel_rho[is][ir] = std::pow(prho[is][ir], exponent);
        }

        // Step 2: H → D
        syncmem_d2d_h2d_op()(d_rho_, rkernel_rho[is], nrxx);

        // Step 3: Real → Complex on GPU
        kedf_wt_real_to_complex<<<blocks, THREADS_PER_BLOCK>>>(
            d_rho_,
            reinterpret_cast<thrust::complex<double>*>(d_fft),
            nrxx);
        CHECK_CUDA_SYNC();

        // Step 4: Forward FFT (in-place on d_fft)
        CUFFT_CHECK(cufftExecZ2Z(cufft_plan_fwd_, d_fft, d_fft, CUFFT_FORWARD));

        // Step 5: Multiply by WT kernel in G-space
        kedf_wt_recip_multiply<<<blocks, THREADS_PER_BLOCK>>>(
            reinterpret_cast<thrust::complex<double>*>(d_fft),
            d_kernel_,
            npw);
        CHECK_CUDA_SYNC();

        // Step 6: Inverse FFT (in-place on d_fft)
        CUFFT_CHECK(cufftExecZ2Z(cufft_plan_bwd_, d_fft, d_fft, CUFFT_INVERSE));

        // Step 7: Complex → Real with 1/N normalization
        kedf_wt_complex_to_real_norm<<<blocks, THREADS_PER_BLOCK>>>(
            reinterpret_cast<thrust::complex<double>*>(d_fft),
            d_result_,  // back to double* usage
            inv_nrxx,
            nrxx);
        CHECK_CUDA_SYNC();

        // Step 8: D → H
        syncmem_d2d_d2h_op()(rkernel_rho[is], d_result_, nrxx);
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

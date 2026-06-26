/**
 * @file diago_kernels.cu
 * @brief CUDA kernel implementations for eigenvalue solver GPU acceleration.
 *
 * Implements highly optimized CUDA kernels for:
 *   1. Batched dot products with warp-level reductions
 *   2. Coalesced vector-prececonditioner division
 *   3. Fused CG gradient computation
 *   4. Schmidt orthogonalization
 *   5. Wavefunction subspace update (via cuBLAS)
 *   6. Band energy computation
 *   7. Preconditioner application
 *
 * Design principles:
 *   - Warp-level parallel reduction for dot products
 *   - Coalesced global memory access (threads in a warp access consecutive
 *     bands at the same basis index for batched operations)
 *   - Kernel fusion to minimize global memory round-trips
 *   - CUDA stream support for overlapping computation and data transfer
 */

#include "diago_kernels.cuh"

#include <cstdio>
#include <cassert>

namespace hsolver {
namespace cuda {

// ============================================================================
// Device constants
// ============================================================================
static constexpr int WARP_SIZE = 32;
static constexpr int MAX_THREADS_PER_BLOCK = 256;
static constexpr int FULL_MASK = 0xffffffff;

// ============================================================================
// Global resource handles (one per process)
// ============================================================================
static cublasHandle_t g_cublas_handle = nullptr;
static cudaStream_t  g_compute_stream = nullptr;
static cudaStream_t  g_transfer_stream = nullptr;
static bool          g_resources_initialized = false;

void init_diago_cuda_resources()
{
    if (g_resources_initialized) return;

    CHECK_CUBLAS(cublasCreate(&g_cublas_handle));
    CHECK_CUDA(cudaStreamCreate(&g_compute_stream));
    CHECK_CUDA(cudaStreamCreate(&g_transfer_stream));

    g_resources_initialized = true;
}

void destroy_diago_cuda_resources()
{
    if (!g_resources_initialized) return;

    if (g_cublas_handle)   { cublasDestroy(g_cublas_handle); g_cublas_handle = nullptr; }
    if (g_compute_stream)  { cudaStreamDestroy(g_compute_stream); g_compute_stream = nullptr; }
    if (g_transfer_stream) { cudaStreamDestroy(g_transfer_stream); g_transfer_stream = nullptr; }

    g_resources_initialized = false;
}

// ============================================================================
// Helper: warp-level reduction (sum)
// ============================================================================
template <typename Real>
__device__ __forceinline__
Real warp_reduce_sum(Real val)
{
    #pragma unroll
    for (int offset = WARP_SIZE / 2; offset > 0; offset >>= 1)
    {
        val += __shfl_down_sync(FULL_MASK, val, offset);
    }
    return val;
}

// ============================================================================
// Helper: block-level reduction (sum)
// ============================================================================
template <typename Real, int BLOCK_SIZE>
__device__ __forceinline__
Real block_reduce_sum(Real val, Real* shared, int tid)
{
    int lane = tid % WARP_SIZE;
    int warp_id = tid / WARP_SIZE;

    // Step 1: warp-level reduction
    val = warp_reduce_sum<Real>(val);

    // Step 2: write warp result to shared memory
    if (lane == 0) {
        shared[warp_id] = val;
    }
    __syncthreads();

    // Step 3: first warp reduces warp results
    constexpr int num_warps = BLOCK_SIZE / WARP_SIZE;
    if (warp_id == 0)
    {
        val = (lane < num_warps) ? shared[lane] : static_cast<Real>(0);
        val = warp_reduce_sum<Real>(val);
    }
    __syncthreads();

    return val;
}

// ============================================================================
// Kernel 1: Batched Dot Product
// ============================================================================

template <typename Real>
__global__ void batched_dot_real_kernel(
    Real* __restrict__ result,
    const thrust::complex<Real>* __restrict__ psi_L,
    const thrust::complex<Real>* __restrict__ psi_R,
    int n_basis,
    int ld_psi,
    int n_band)
{
    // Each block handles one band
    int band_idx = blockIdx.x;
    if (band_idx >= n_band) return;

    const int tid = threadIdx.x;
    constexpr int BLOCK_SIZE = MAX_THREADS_PER_BLOCK;

    __shared__ Real shared[BLOCK_SIZE / WARP_SIZE];

    Real sum = static_cast<Real>(0);

    // Grid-stride loop for basis elements (coalesced since each thread
    // processes different bands, same basis index — but here we do
    // one block per band, stride over basis)
    const thrust::complex<Real>* L = psi_L + band_idx * ld_psi;
    const thrust::complex<Real>* R = psi_R + band_idx * ld_psi;

    for (int i = tid; i < n_basis; i += BLOCK_SIZE)
    {
        // dot = conj(L[i]) * R[i]
        sum += L[i].real() * R[i].real() + L[i].imag() * R[i].imag();
    }

    sum = block_reduce_sum<Real, BLOCK_SIZE>(sum, shared, tid);

    if (tid == 0) {
        result[band_idx] = sum;
    }
}

template <typename Real>
void batched_dot_real_op(
    Real* result,
    const thrust::complex<Real>* psi_L,
    const thrust::complex<Real>* psi_R,
    int n_basis,
    int ld_psi,
    int n_band,
    cudaStream_t stream)
{
    int blocks = n_band;
    int threads = MAX_THREADS_PER_BLOCK;
    batched_dot_real_kernel<Real><<<blocks, threads, 0, stream>>>(
        result, psi_L, psi_R, n_basis, ld_psi, n_band);
}

// ============================================================================
// Kernel 2: Batched Vector-Prececonditioner Division (coalesced access)
// ============================================================================

template <typename Real>
__global__ void batched_div_preconditioner_kernel(
    thrust::complex<Real>* __restrict__ out,
    const thrust::complex<Real>* __restrict__ in_vec,
    const Real* __restrict__ prec,
    int n_basis,
    int ld,
    int n_band)
{
    // Coalesced access: threads in a warp access consecutive bands at
    // the same basis index.
    int basis_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (basis_idx >= n_basis) return;

    Real inv_prec = static_cast<Real>(1.0) / prec[basis_idx];

    for (int b = 0; b < n_band; ++b)
    {
        int idx = b * ld + basis_idx;
        out[idx].real(in_vec[idx].real() * inv_prec);
        out[idx].imag(in_vec[idx].imag() * inv_prec);
    }
}

template <typename Real>
void batched_div_preconditioner_op(
    thrust::complex<Real>* out,
    const thrust::complex<Real>* in_vec,
    const Real* prec,
    int n_basis,
    int ld,
    int n_band,
    cudaStream_t stream)
{
    int threads = MAX_THREADS_PER_BLOCK;
    int blocks = (n_basis + threads - 1) / threads;
    batched_div_preconditioner_kernel<Real><<<blocks, threads, 0, stream>>>(
        out, in_vec, prec, n_basis, ld, n_band);
}

// ============================================================================
// Kernel 3: Combined CG Gradient Computation (kernel fusion)
// ============================================================================

template <typename Real>
__global__ void calc_grad_cg_kernel(
    thrust::complex<Real>* __restrict__ grad,
    thrust::complex<Real>* __restrict__ pphi,
    const thrust::complex<Real>* __restrict__ hphi,
    const thrust::complex<Real>* __restrict__ sphi,
    const Real* __restrict__ prec,
    int n_basis,
    Real* __restrict__ global_dot_result) // [2]: eh, es
{
    const int tid = threadIdx.x;
    constexpr int BLOCK_SIZE = MAX_THREADS_PER_BLOCK;

    __shared__ Real shared_eh[BLOCK_SIZE / WARP_SIZE];
    __shared__ Real shared_es[BLOCK_SIZE / WARP_SIZE];

    Real eh_local = static_cast<Real>(0);
    Real es_local = static_cast<Real>(0);

    // Step 1: divide by preconditioner and accumulate dot products
    for (int i = tid; i < n_basis; i += BLOCK_SIZE)
    {
        Real inv_p = static_cast<Real>(1.0) / prec[i];

        // grad[i] = hphi[i] / prec[i]
        Real g_re = hphi[i].real() * inv_p;
        Real g_im = hphi[i].imag() * inv_p;
        grad[i].real(g_re);
        grad[i].imag(g_im);

        // pphi[i] = sphi[i] / prec[i]
        Real p_re = sphi[i].real() * inv_p;
        Real p_im = sphi[i].imag() * inv_p;
        pphi[i].real(p_re);
        pphi[i].imag(p_im);

        // eh += conj(sphi[i]) * grad[i] = sphi_re * g_re + sphi_im * g_im
        eh_local += sphi[i].real() * g_re + sphi[i].imag() * g_im;
        // es += conj(sphi[i]) * pphi[i] = sphi_re * p_re + sphi_im * p_im
        es_local += sphi[i].real() * p_re + sphi[i].imag() * p_im;
    }

    // Block reduction for eh and es
    int lane = tid % WARP_SIZE;
    int warp_id = tid / WARP_SIZE;

    eh_local = warp_reduce_sum<Real>(eh_local);
    es_local = warp_reduce_sum<Real>(es_local);

    if (lane == 0) {
        shared_eh[warp_id] = eh_local;
        shared_es[warp_id] = es_local;
    }
    __syncthreads();

    constexpr int num_warps = BLOCK_SIZE / WARP_SIZE;
    if (warp_id == 0)
    {
        eh_local = (lane < num_warps) ? shared_eh[lane] : static_cast<Real>(0);
        es_local = (lane < num_warps) ? shared_es[lane] : static_cast<Real>(0);
        eh_local = warp_reduce_sum<Real>(eh_local);
        es_local = warp_reduce_sum<Real>(es_local);

        if (lane == 0) {
            global_dot_result[0] = eh_local;
            global_dot_result[1] = es_local;
        }
    }
    __syncthreads();

    Real lambda = global_dot_result[0] / global_dot_result[1];

    // Step 2: grad[i] = grad[i] - lambda * pphi[i]
    for (int i = tid; i < n_basis; i += BLOCK_SIZE)
    {
        grad[i].real(grad[i].real() - lambda * pphi[i].real());
        grad[i].imag(grad[i].imag() - lambda * pphi[i].imag());
    }
}

template <typename Real>
void calc_grad_cg_op(
    thrust::complex<Real>* grad,
    thrust::complex<Real>* pphi,
    const thrust::complex<Real>* hphi,
    const thrust::complex<Real>* sphi,
    const Real* prec,
    int n_basis,
    cudaStream_t stream)
{
    // Allocate 2-element buffer for dot results (eh, es)
    Real* d_dot_buffer = nullptr;
    CHECK_CUDA(cudaMallocAsync(&d_dot_buffer, 2 * sizeof(Real), stream));

    int threads = MAX_THREADS_PER_BLOCK;
    calc_grad_cg_kernel<Real><<<1, threads, 0, stream>>>(
        grad, pphi, hphi, sphi, prec, n_basis, d_dot_buffer);

    CHECK_CUDA(cudaFreeAsync(d_dot_buffer, stream));
}

// ============================================================================
// Kernel 4: Schmidt Orthogonalization for CG
// ============================================================================

template <typename Real>
__global__ void schmidt_orth_cg_kernel_step1(
    thrust::complex<Real>* __restrict__ lagrange,
    const thrust::complex<Real>* __restrict__ psi,
    const thrust::complex<Real>* __restrict__ scg,
    int n_basis,
    int ld_psi,
    int m)
{
    // Each block handles one band j for computing lagrange[j]
    int j = blockIdx.x;
    if (j >= m) return;

    const int tid = threadIdx.x;
    constexpr int BLOCK_SIZE = MAX_THREADS_PER_BLOCK;

    __shared__ Real shared[BLOCK_SIZE / WARP_SIZE];

    Real sum_re = static_cast<Real>(0);
    Real sum_im = static_cast<Real>(0);

    const thrust::complex<Real>* psi_j = psi + j * ld_psi;

    for (int i = tid; i < n_basis; i += BLOCK_SIZE)
    {
        // lagrange[j] = sum conj(psi_j[i]) * scg[i]
        sum_re += psi_j[i].real() * scg[i].real() + psi_j[i].imag() * scg[i].imag();
        sum_im += psi_j[i].real() * scg[i].imag() - psi_j[i].imag() * scg[i].real();
    }

    int lane = tid % WARP_SIZE;
    int warp_id = tid / WARP_SIZE;

    sum_re = warp_reduce_sum<Real>(sum_re);
    sum_im = warp_reduce_sum<Real>(sum_im);

    if (lane == 0) {
        shared[warp_id] = sum_re;
    }
    __syncthreads();

    constexpr int num_warps = BLOCK_SIZE / WARP_SIZE;
    if (warp_id == 0)
    {
        sum_re = (lane < num_warps) ? shared[lane] : static_cast<Real>(0);
        sum_re = warp_reduce_sum<Real>(sum_re);
        if (lane == 0) lagrange[j].real(sum_re);
    }
    // For imaginary part (need separate shared memory buffer)
    // Simplified: use another shared buffer slot
    __syncthreads();
    if (lane == 0) shared[warp_id] = sum_im;
    __syncthreads();
    if (warp_id == 0)
    {
        sum_im = (lane < num_warps) ? shared[lane] : static_cast<Real>(0);
        sum_im = warp_reduce_sum<Real>(sum_im);
        if (lane == 0) lagrange[j].imag(sum_im);
    }
}

template <typename Real>
__global__ void schmidt_orth_cg_kernel_step2(
    thrust::complex<Real>* __restrict__ grad,
    thrust::complex<Real>* __restrict__ scg,
    const thrust::complex<Real>* __restrict__ psi,
    const thrust::complex<Real>* __restrict__ spsi,
    const thrust::complex<Real>* __restrict__ lagrange,
    int n_basis,
    int ld_psi,
    int m)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_basis) return;

    thrust::complex<Real> g_correction(static_cast<Real>(0), static_cast<Real>(0));
    thrust::complex<Real> s_correction(static_cast<Real>(0), static_cast<Real>(0));

    for (int j = 0; j < m; ++j)
    {
        thrust::complex<Real> lag = lagrange[j];
        const thrust::complex<Real>* psi_j = psi + j * ld_psi;
        const thrust::complex<Real>* spsi_j = spsi + j * ld_psi;

        // g_correction += lagrange[j] * psi[j][i]
        g_correction.real(g_correction.real() + lag.real() * psi_j[i].real()
                          - lag.imag() * psi_j[i].imag());
        g_correction.imag(g_correction.imag() + lag.real() * psi_j[i].imag()
                          + lag.imag() * psi_j[i].real());

        // s_correction += lagrange[j] * spsi[j][i]
        s_correction.real(s_correction.real() + lag.real() * spsi_j[i].real()
                          - lag.imag() * spsi_j[i].imag());
        s_correction.imag(s_correction.imag() + lag.real() * spsi_j[i].imag()
                          + lag.imag() * spsi_j[i].real());
    }

    grad[i].real(grad[i].real() - g_correction.real());
    grad[i].imag(grad[i].imag() - g_correction.imag());
    scg[i].real(scg[i].real() - s_correction.real());
    scg[i].imag(scg[i].imag() - s_correction.imag());
}

template <typename Real>
void schmidt_orth_cg_op(
    thrust::complex<Real>* grad,
    thrust::complex<Real>* scg,
    const thrust::complex<Real>* psi,
    const thrust::complex<Real>* spsi,
    thrust::complex<Real>* lagrange,
    int n_basis,
    int ld_psi,
    int m,
    cudaStream_t stream)
{
    // Step 1: compute lagrange multipliers (one block per band)
    int threads = MAX_THREADS_PER_BLOCK;
    schmidt_orth_cg_kernel_step1<Real><<<m, threads, 0, stream>>>(
        lagrange, psi, scg, n_basis, ld_psi, m);

    // Step 2: apply orthogonalization
    int blocks = (n_basis + threads - 1) / threads;
    schmidt_orth_cg_kernel_step2<Real><<<blocks, threads, 0, stream>>>(
        grad, scg, psi, spsi, lagrange, n_basis, ld_psi, m);
}

// ============================================================================
// Kernel 5: Wavefunction Subspace Update (cuBLAS GEMM)
// ============================================================================

template <typename T>
void subspace_update_op(
    T* psi_out,
    const T* psi,
    const T* vcc,
    int dim,
    int nbase,
    int nband,
    int ld_psi,
    int ld_vcc,
    cudaStream_t stream)
{
    cublasSetStream(g_cublas_handle, stream);

    if (std::is_same<T, thrust::complex<double>>::value)
    {
        const cuDoubleComplex alpha = {1.0, 0.0};
        const cuDoubleComplex beta  = {0.0, 0.0};
        CHECK_CUBLAS(cublasZgemm(g_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
            dim, nband, nbase, &alpha,
            reinterpret_cast<const cuDoubleComplex*>(psi), ld_psi,
            reinterpret_cast<const cuDoubleComplex*>(vcc), ld_vcc,
            &beta,
            reinterpret_cast<cuDoubleComplex*>(psi_out), ld_psi));
    }
    else if (std::is_same<T, thrust::complex<float>>::value)
    {
        const cuComplex alpha = {1.0f, 0.0f};
        const cuComplex beta  = {0.0f, 0.0f};
        CHECK_CUBLAS(cublasCgemm(g_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
            dim, nband, nbase, &alpha,
            reinterpret_cast<const cuComplex*>(psi), ld_psi,
            reinterpret_cast<const cuComplex*>(vcc), ld_vcc,
            &beta,
            reinterpret_cast<cuComplex*>(psi_out), ld_psi));
    }
    else
    {
        const double alpha = 1.0;
        const double beta  = 0.0;
        CHECK_CUBLAS(cublasDgemm(g_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
            dim, nband, nbase, &alpha,
            reinterpret_cast<const double*>(psi), ld_psi,
            reinterpret_cast<const double*>(vcc), ld_vcc,
            &beta,
            reinterpret_cast<double*>(psi_out), ld_psi));
    }
}

// ============================================================================
// Kernel 6: Band Energy Computation (Batched Rayleigh Quotients)
// ============================================================================

template <typename Real>
__global__ void compute_band_energies_kernel(
    Real* __restrict__ eigen,
    const thrust::complex<Real>* __restrict__ psi,
    const thrust::complex<Real>* __restrict__ hpsi,
    int n_basis,
    int ld_psi,
    int n_band)
{
    int band_idx = blockIdx.x;
    if (band_idx >= n_band) return;

    const int tid = threadIdx.x;
    constexpr int BLOCK_SIZE = MAX_THREADS_PER_BLOCK;

    __shared__ Real shared[BLOCK_SIZE / WARP_SIZE];

    const thrust::complex<Real>* p = psi + band_idx * ld_psi;
    const thrust::complex<Real>* hp = hpsi + band_idx * ld_psi;

    Real sum = static_cast<Real>(0);
    for (int i = tid; i < n_basis; i += BLOCK_SIZE)
    {
        sum += p[i].real() * hp[i].real() + p[i].imag() * hp[i].imag();
    }

    sum = block_reduce_sum<Real, BLOCK_SIZE>(sum, shared, tid);

    if (tid == 0) {
        eigen[band_idx] = sum;
    }
}

template <typename Real>
void compute_band_energies_op(
    Real* eigen,
    const thrust::complex<Real>* psi,
    const thrust::complex<Real>* hpsi,
    int n_basis,
    int ld_psi,
    int n_band,
    cudaStream_t stream)
{
    int blocks = n_band;
    int threads = MAX_THREADS_PER_BLOCK;
    compute_band_energies_kernel<Real><<<blocks, threads, 0, stream>>>(
        eigen, psi, hpsi, n_basis, ld_psi, n_band);
}

// ============================================================================
// Kernel 7: Preconditioner Application (Residual + Preconditioner)
// ============================================================================

template <typename Real>
__global__ void apply_preconditioner_kernel(
    thrust::complex<Real>* __restrict__ grad,
    const thrust::complex<Real>* __restrict__ hpsi,
    const thrust::complex<Real>* __restrict__ spsi,
    const Real* __restrict__ prec,
    Real eigen,
    int n_basis)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_basis) return;

    Real inv_p = static_cast<Real>(1.0) / prec[i];

    // grad[i] = (hpsi[i] - eigen * spsi[i]) / prec[i]
    grad[i].real((hpsi[i].real() - eigen * spsi[i].real()) * inv_p);
    grad[i].imag((hpsi[i].imag() - eigen * spsi[i].imag()) * inv_p);
}

template <typename Real>
void apply_preconditioner_op(
    thrust::complex<Real>* grad,
    const thrust::complex<Real>* hpsi,
    const thrust::complex<Real>* spsi,
    const Real* prec,
    Real eigen,
    int n_basis,
    cudaStream_t stream)
{
    int threads = MAX_THREADS_PER_BLOCK;
    int blocks = (n_basis + threads - 1) / threads;
    apply_preconditioner_kernel<Real><<<blocks, threads, 0, stream>>>(
        grad, hpsi, spsi, prec, eigen, n_basis);
}

// ============================================================================
// Explicit template instantiations
// ============================================================================

// Batched dot product
template void batched_dot_real_op<float>(
    float*, const thrust::complex<float>*, const thrust::complex<float>*,
    int, int, int, cudaStream_t);
template void batched_dot_real_op<double>(
    double*, const thrust::complex<double>*, const thrust::complex<double>*,
    int, int, int, cudaStream_t);

// Batched div preconditioner
template void batched_div_preconditioner_op<float>(
    thrust::complex<float>*, const thrust::complex<float>*,
    const float*, int, int, int, cudaStream_t);
template void batched_div_preconditioner_op<double>(
    thrust::complex<double>*, const thrust::complex<double>*,
    const double*, int, int, int, cudaStream_t);

// Calc grad CG
template void calc_grad_cg_op<float>(
    thrust::complex<float>*, thrust::complex<float>*,
    const thrust::complex<float>*, const thrust::complex<float>*,
    const float*, int, cudaStream_t);
template void calc_grad_cg_op<double>(
    thrust::complex<double>*, thrust::complex<double>*,
    const thrust::complex<double>*, const thrust::complex<double>*,
    const double*, int, cudaStream_t);

// Schmidt orth CG
template void schmidt_orth_cg_op<float>(
    thrust::complex<float>*, thrust::complex<float>*,
    const thrust::complex<float>*, const thrust::complex<float>*,
    thrust::complex<float>*, int, int, int, cudaStream_t);
template void schmidt_orth_cg_op<double>(
    thrust::complex<double>*, thrust::complex<double>*,
    const thrust::complex<double>*, const thrust::complex<double>*,
    thrust::complex<double>*, int, int, int, cudaStream_t);

// Subspace update
template void subspace_update_op<thrust::complex<float>>(
    thrust::complex<float>*, const thrust::complex<float>*,
    const thrust::complex<float>*, int, int, int, int, int, cudaStream_t);
template void subspace_update_op<thrust::complex<double>>(
    thrust::complex<double>*, const thrust::complex<double>*,
    const thrust::complex<double>*, int, int, int, int, int, cudaStream_t);
template void subspace_update_op<double>(
    double*, const double*, const double*, int, int, int, int, int, cudaStream_t);

// Band energies
template void compute_band_energies_op<float>(
    float*, const thrust::complex<float>*, const thrust::complex<float>*,
    int, int, int, cudaStream_t);
template void compute_band_energies_op<double>(
    double*, const thrust::complex<double>*, const thrust::complex<double>*,
    int, int, int, cudaStream_t);

// Apply preconditioner
template void apply_preconditioner_op<float>(
    thrust::complex<float>*, const thrust::complex<float>*,
    const thrust::complex<float>*, const float*, float, int, cudaStream_t);
template void apply_preconditioner_op<double>(
    thrust::complex<double>*, const thrust::complex<double>*,
    const thrust::complex<double>*, const double*, double, int, cudaStream_t);

} // namespace cuda
} // namespace hsolver

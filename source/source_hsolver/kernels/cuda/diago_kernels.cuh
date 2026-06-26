/**
 * @file diago_kernels.cuh
 * @brief CUDA kernel declarations for eigenvalue solver GPU acceleration.
 *
 * This module provides GPU-accelerated kernels for the most compute-intensive
 * operations in the eigenvalue solvers (CG, Davidson, BPCG):
 *   - Batched dot products (band-parallel reductions)
 *   - Vector-prececonditioner division
 *   - Combined gradient computation
 *   - Schmidt orthogonalization
 *   - Wavefunction subspace update (Davidson)
 *
 * All kernels use warp-level reductions and coalesced memory access patterns.
 */

#ifndef MODULE_HSOLVER_DIAGO_KERNELS_CUH
#define MODULE_HSOLVER_DIAGO_KERNELS_CUH

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <thrust/complex.h>

#include "source_base/module_device/device_check.h"

namespace hsolver {
namespace cuda {

// ============================================================================
// Handle management for CUDA streams and cuBLAS
// ============================================================================

/**
 * @brief Initialize CUDA solver resources (streams, cuBLAS handle).
 * Thread-safe. Only creates resources on the first call.
 */
void init_diago_cuda_resources();

/**
 * @brief Destroy CUDA solver resources.
 * Thread-safe.
 */
void destroy_diago_cuda_resources();

// ============================================================================
// Kernel 1: Batched Dot Product
// ============================================================================

/**
 * @brief Compute dot products for multiple bands in a single kernel launch.
 *
 * For each band m (0..n_band-1), computes:
 *   result[m] = sum_{i=0}^{n_basis-1} conj(psi_L[m * ld_psi + i]) * psi_R[m * ld_psi + i]
 *
 * @tparam Real  Real floating point type (float or double).
 * @param result    Output array of size n_band (device pointer).
 * @param psi_L     Left operand blockvector (device pointer, column-major layout).
 * @param psi_R     Right operand blockvector (device pointer, column-major layout).
 * @param n_basis   Number of basis elements per band.
 * @param ld_psi    Leading dimension (stride between consecutive bands).
 * @param n_band    Number of bands to process.
 * @param stream    CUDA stream for asynchronous execution.
 */
template <typename Real>
void batched_dot_real_op(
    Real* result,
    const thrust::complex<Real>* psi_L,
    const thrust::complex<Real>* psi_R,
    int n_basis,
    int ld_psi,
    int n_band,
    cudaStream_t stream = 0);

// ============================================================================
// Kernel 2: Vector-Prececonditioner Division
// ============================================================================

/**
 * @brief Divide multiple vectors by preconditioner elements.
 *
 * For each band m and basis i:
 *   out[m * ld + i] = in_vec[m * ld + i] / prec[i]
 *
 * This kernel uses coalesced reads: threads within a warp access consecutive
 * bands at the same basis index, achieving coalesced global memory access.
 *
 * @tparam T       Complex type (thrust::complex<float> or thrust::complex<double>).
 * @param out      Output array (device pointer).
 * @param in_vec   Input array (device pointer).
 * @param prec     Preconditioner array of size n_basis (device pointer, real type).
 * @param n_basis  Number of basis elements.
 * @param ld       Leading dimension (stride between bands).
 * @param n_band   Number of bands.
 * @param stream   CUDA stream.
 */
template <typename Real>
void batched_div_preconditioner_op(
    thrust::complex<Real>* out,
    const thrust::complex<Real>* in_vec,
    const Real* prec,
    int n_basis,
    int ld,
    int n_band,
    cudaStream_t stream = 0);

// ============================================================================
// Kernel 3: Combined CG Gradient Computation
// ============================================================================

/**
 * @brief Combined gradient computation for the CG eigenvalue solver.
 *
 * Performs in a single kernel launch:
 *   1. grad[i] = hphi[i] / prec[i]
 *   2. pphi[i] = sphi[i] / prec[i]
 *   3. eh = sum conj(sphi[i]) * grad[i]
 *   4. es = sum conj(sphi[i]) * pphi[i]
 *   5. lambda = eh / es
 *   6. grad[i] = grad[i] - lambda * pphi[i]
 *
 * This fusion eliminates multiple kernel launches and reduces global memory traffic.
 *
 * @tparam Real  Real floating point type.
 * @param grad      In/out gradient vector (device pointer, n_basis elements).
 * @param pphi      Output PS|psi> vector (device pointer, n_basis elements).
 * @param hphi      Input H|psi> vector (device pointer).
 * @param sphi      Input S|psi> vector (device pointer).
 * @param prec      Preconditioner array (device pointer).
 * @param n_basis   Number of basis elements.
 * @param stream    CUDA stream.
 */
template <typename Real>
void calc_grad_cg_op(
    thrust::complex<Real>* grad,
    thrust::complex<Real>* pphi,
    const thrust::complex<Real>* hphi,
    const thrust::complex<Real>* sphi,
    const Real* prec,
    int n_basis,
    cudaStream_t stream = 0);

// ============================================================================
// Kernel 4: Schmidt Orthogonalization for CG
// ============================================================================

/**
 * @brief Schmidt orthogonalization for the CG solver on GPU.
 *
 * Given the gradient vector g and the existing psi vectors (0..m-1),
 * performs:
 *   1. scg = S * g  (handled externally via spsi_func)
 *   2. lagrange[j] = sum_{i} conj(psi[j * ld + i]) * scg[i]  for j=0..m-1
 *   3. g[i] = g[i] - sum_{j} lagrange[j] * psi[j * ld + i]
 *   4. scg[i] = scg[i] - sum_{j} lagrange[j] * spsi[j * ld + i]
 *
 * @tparam Real  Real floating point type.
 * @param grad      In/out gradient vector (device pointer).
 * @param scg       In/out S|grad> vector (device pointer).
 * @param psi       Existing psi blockvector (device pointer).
 * @param spsi      Existing S|psi> blockvector (device pointer).
 * @param lagrange  Output Lagrange multipliers (device pointer, size m).
 * @param n_basis   Number of basis elements.
 * @param ld_psi    Leading dimension of psi and spsi.
 * @param m         Number of existing bands to orthogonalize against.
 * @param stream    CUDA stream.
 */
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
    cudaStream_t stream = 0);

// ============================================================================
// Kernel 5: Wavefunction Subspace Update (Davidson)
// ============================================================================

/**
 * @brief Update the wavefunction subspace for the Davidson solver.
 *
 * Computes: psi_out = psi * vcc  (matrix multiplication)
 *   psi_out(dim, nband) = psi(dim, nbase) * vcc(nbase, nband)
 *
 * Uses cuBLAS gemm for performance.
 *
 * @tparam T       Complex type.
 * @param psi_out  Output wavefunction (device pointer).
 * @param psi      Input basis vectors (device pointer).
 * @param vcc      Eigenvectors of reduced Hamiltonian (device pointer).
 * @param dim      Dimension of each basis vector.
 * @param nbase    Current dimension of reduced basis.
 * @param nband    Number of bands.
 * @param ld_psi   Leading dimension of psi.
 * @param ld_vcc   Leading dimension of vcc.
 * @param stream   CUDA stream.
 */
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
    cudaStream_t stream);

// ============================================================================
// Kernel 6: Band Energy Computation (Batched Eigenvalue Update)
// ============================================================================

/**
 * @brief Compute band energies (Rayleigh quotients) for all bands.
 *
 * For each band m: eigen[m] = sum conj(psi[m * ld + i]) * hpsi[m * ld + i]
 *
 * Uses batched dot product reduction.
 *
 * @tparam Real  Real floating point type.
 * @param eigen     Output eigenvalues (device pointer, size n_band).
 * @param psi       Wavefunction blockvector (device pointer).
 * @param hpsi      H|psi> blockvector (device pointer).
 * @param n_basis   Number of basis elements.
 * @param ld_psi    Leading dimension.
 * @param n_band    Number of bands.
 * @param stream    CUDA stream.
 */
template <typename Real>
void compute_band_energies_op(
    Real* eigen,
    const thrust::complex<Real>* psi,
    const thrust::complex<Real>* hpsi,
    int n_basis,
    int ld_psi,
    int n_band,
    cudaStream_t stream = 0);

// ============================================================================
// Kernel 7: Preconditioner Application (Batched AXPY variant)
// ============================================================================

/**
 * @brief Apply preconditioner and compute gradient: g = K^{-1} * (H - e*S) * psi
 *
 * This combines the residual computation and preconditioner application:
 *   g[i] = (hpsi[i] - e * spsi[i]) / prec[i]
 *
 * @tparam Real  Real floating point type.
 * @param grad      Output gradient (device pointer).
 * @param hpsi      H|psi> (device pointer).
 * @param spsi      S|psi> (device pointer).
 * @param prec      Preconditioner (device pointer).
 * @param eigen     Eigenvalue e.
 * @param n_basis   Number of basis elements.
 * @param stream    CUDA stream.
 */
template <typename Real>
void apply_preconditioner_op(
    thrust::complex<Real>* grad,
    const thrust::complex<Real>* hpsi,
    const thrust::complex<Real>* spsi,
    const Real* prec,
    Real eigen,
    int n_basis,
    cudaStream_t stream = 0);

} // namespace cuda
} // namespace hsolver

#endif // MODULE_HSOLVER_DIAGO_KERNELS_CUH

/**
 * @file diago_cg_gpu.h
 * @brief GPU-accelerated CG eigenvalue solver helper.
 *
 * This module provides GPU-accelerated implementations of the key
 * operations in the CG eigenvalue solver:
 *   - calc_grad_gpu: GPU version of gradient computation
 *   - orth_grad_gpu: GPU version of Schmidt orthogonalization
 *   - update_psi_gpu: GPU version of wavefunction update
 *   - band_energies_gpu: GPU batched band energy computation
 *
 * The GPU path is automatically enabled when __CUDA is defined and
 * the solver is instantiated with a GPU device type.
 *
 * @note All GPU operations use CUDA streams for asynchronous execution
 *       and overlap of computation with data transfer.
 */

#ifndef MODULE_HSOLVER_DIAGO_CG_GPU_H
#define MODULE_HSOLVER_DIAGO_CG_GPU_H

#include "source_hsolver/kernels/cuda/diago_kernels.cuh"

#include <complex>
#include <cuda_runtime.h>

namespace hsolver {

/**
 * @class DiagoCG_GPUHelper
 * @brief GPU-accelerated operations for the CG eigenvalue solver.
 *
 * This class encapsulates all GPU-accelerated operations used by the
 * CG solver. It manages GPU memory allocation/deallocation for temporary
 * buffers and provides optimized fused kernels.
 *
 * @tparam T  Complex data type (std::complex<float> or std::complex<double>).
 */
template <typename T>
class DiagoCG_GPUHelper
{
public:
    using Real = typename GetTypeReal<T>::type;

    DiagoCG_GPUHelper()
    {
        init_diago_cuda_resources();
    }

    ~DiagoCG_GPUHelper()
    {
        free_buffers();
    }

    // ---- GPU memory management ----

    /**
     * @brief Ensure temporary GPU buffers are allocated for the given problem size.
     * @param n_basis  Number of basis elements per band.
     * @param n_band   Number of bands.
     */
    void ensure_buffers(int n_basis, int n_band)
    {
        if (n_basis_ != n_basis || n_band_ != n_band)
        {
            free_buffers();
            n_basis_ = n_basis;
            n_band_ = n_band;

            // Allocate workspace for CG operations
            size_t complex_bytes = sizeof(T) * n_basis;
            size_t real_bytes = sizeof(Real) * n_basis;
            size_t lagrange_bytes = sizeof(T) * n_band;

            CHECK_CUDA(cudaMalloc(&d_grad_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_pphi_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_hphi_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_sphi_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_scg_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_cg_, complex_bytes));
            CHECK_CUDA(cudaMalloc(&d_prec_, real_bytes));
            CHECK_CUDA(cudaMalloc(&d_lagrange_, lagrange_bytes));
            CHECK_CUDA(cudaMalloc(&d_dot_workspace_, 2 * sizeof(Real)));
        }
    }

    /**
     * @brief Copy preconditioner data to GPU.
     * @param prec  Host preconditioner array of size n_basis.
     */
    void upload_preconditioner(const Real* prec)
    {
        CHECK_CUDA(cudaMemcpyAsync(d_prec_, prec, sizeof(Real) * n_basis_,
                                   cudaMemcpyHostToDevice, gpu_stream()));
    }

    // ---- GPU-accelerated CG operations ----

    /**
     * @brief GPU-accelerated gradient computation for CG solver.
     *
     * Performs:
     *   1. grad = hphi / prec
     *   2. pphi = sphi / prec
     *   3. eh = <sphi | grad>   (Rayleigh quotient numerator)
     *   4. es = <sphi | pphi>   (Rayleigh quotient denominator)
     *   5. lambda = eh / es
     *   6. grad = grad - lambda * pphi
     *
     * @param grad   Output gradient vector (device pointer, size n_basis).
     * @param pphi   Output preconditioned S|psi> vector (device pointer).
     * @param hphi   Input H|psi> vector (device pointer).
     * @param sphi   Input S|psi> vector (device pointer).
     * @param prec   Preconditioner array (device pointer, size n_basis).
     * @param n_basis Number of basis elements.
     */
    void calc_grad(T* grad, T* pphi,
                   const T* hphi, const T* sphi,
                   const Real* prec, int n_basis)
    {
        auto* cuda_stream = gpu_stream();

        calc_grad_cg_op<Real>(
            reinterpret_cast<thrust::complex<Real>*>(grad),
            reinterpret_cast<thrust::complex<Real>*>(pphi),
            reinterpret_cast<const thrust::complex<Real>*>(hphi),
            reinterpret_cast<const thrust::complex<Real>*>(sphi),
            prec, n_basis, cuda_stream);
    }

    /**
     * @brief GPU-accelerated Schmidt orthogonalization.
     *
     * Orthogonalizes the gradient vector against existing psi bands (0..m-1).
     *
     * @param grad      In/out gradient vector (device pointer).
     * @param scg       In/out S|grad> vector (device pointer).
     * @param psi       Existing psi blockvector (device pointer, size ld_psi * m).
     * @param spsi      Existing S|psi> blockvector (device pointer).
     * @param lagrange  Output Lagrange multipliers (device pointer, size m).
     * @param n_basis   Number of basis elements.
     * @param ld_psi    Leading dimension of psi/spsi.
     * @param m         Number of existing bands.
     */
    void schmidt_orth(T* grad, T* scg,
                      const T* psi, const T* spsi,
                      T* lagrange,
                      int n_basis, int ld_psi, int m)
    {
        schmidt_orth_cg_op<Real>(
            reinterpret_cast<thrust::complex<Real>*>(grad),
            reinterpret_cast<thrust::complex<Real>*>(scg),
            reinterpret_cast<const thrust::complex<Real>*>(psi),
            reinterpret_cast<const thrust::complex<Real>*>(spsi),
            reinterpret_cast<thrust::complex<Real>*>(lagrange),
            n_basis, ld_psi, m, gpu_stream());
    }

    /**
     * @brief GPU batched computation of band energies (Rayleigh quotients).
     *
     * For each band m: eigen[m] = <psi[m] | hpsi[m]>
     *
     * @param eigen     Output eigenvalues (device pointer, size n_band).
     * @param psi       Wavefunction blockvector (device pointer).
     * @param hpsi      H|psi> blockvector (device pointer).
     * @param n_basis   Number of basis elements.
     * @param ld_psi    Leading dimension.
     * @param n_band    Number of bands.
     */
    void compute_band_energies(Real* eigen,
                               const T* psi, const T* hpsi,
                               int n_basis, int ld_psi, int n_band)
    {
        compute_band_energies_op<Real>(
            eigen,
            reinterpret_cast<const thrust::complex<Real>*>(psi),
            reinterpret_cast<const thrust::complex<Real>*>(hpsi),
            n_basis, ld_psi, n_band, gpu_stream());
    }

    /**
     * @brief GPU application of preconditioner to compute gradient.
     *
     * Computes: grad[i] = (hpsi[i] - eigen * spsi[i]) / prec[i]
     *
     * @param grad   Output gradient (device pointer).
     * @param hpsi   H|psi> (device pointer).
     * @param spsi   S|psi> (device pointer).
     * @param prec   Preconditioner (device pointer).
     * @param eigen  Eigenvalue.
     * @param n_basis Number of basis elements.
     */
    void apply_preconditioner(T* grad,
                              const T* hpsi, const T* spsi,
                              const Real* prec, Real eigen,
                              int n_basis)
    {
        apply_preconditioner_op<Real>(
            reinterpret_cast<thrust::complex<Real>*>(grad),
            reinterpret_cast<const thrust::complex<Real>*>(hpsi),
            reinterpret_cast<const thrust::complex<Real>*>(spsi),
            prec, eigen, n_basis, gpu_stream());
    }

    /**
     * @brief Synchronize the GPU stream (wait for all GPU work to complete).
     */
    void synchronize()
    {
        CHECK_CUDA(cudaStreamSynchronize(gpu_stream()));
    }

    // ---- Accessors for device buffers ----

    T*      d_grad()     const { return d_grad_; }
    T*      d_pphi()     const { return d_pphi_; }
    T*      d_hphi()     const { return d_hphi_; }
    T*      d_sphi()     const { return d_sphi_; }
    T*      d_scg()      const { return d_scg_; }
    T*      d_cg()       const { return d_cg_; }
    Real*   d_prec()     const { return d_prec_; }
    T*      d_lagrange() const { return d_lagrange_; }

private:
    static cudaStream_t gpu_stream()
    {
        // Use the default stream for simplicity; could be specialized
        return 0;
    }

    void free_buffers()
    {
        auto free_if = [](auto*& ptr) {
            if (ptr) { cudaFree(ptr); ptr = nullptr; }
        };
        free_if(d_grad_);
        free_if(d_pphi_);
        free_if(d_hphi_);
        free_if(d_sphi_);
        free_if(d_scg_);
        free_if(d_cg_);
        free_if(d_prec_);
        free_if(d_lagrange_);
        free_if(d_dot_workspace_);
    }

    int n_basis_ = 0;
    int n_band_  = 0;

    T*     d_grad_     = nullptr;
    T*     d_pphi_     = nullptr;
    T*     d_hphi_     = nullptr;
    T*     d_sphi_     = nullptr;
    T*     d_scg_      = nullptr;
    T*     d_cg_       = nullptr;
    Real*  d_prec_     = nullptr;
    T*     d_lagrange_ = nullptr;
    Real*  d_dot_workspace_ = nullptr;
};

} // namespace hsolver

#endif // MODULE_HSOLVER_DIAGO_CG_GPU_H

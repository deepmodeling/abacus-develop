#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/macros.h"
#include "source_base/module_device/device.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"

#include <ATen/core/tensor_types.h>

#include <complex>
#include <functional>
#include <vector>

namespace hsolver
{

/**
 * @class DiagoPPCG
 * @brief A class for diagonalization using the Projected Preconditioned Conjugate Gradient (PPCG) method.
 *
 * PPCG extends the standard band-by-band CG by constructing a small subspace (2D or 3D) per band
 * from the current eigenvector, the preconditioned residual, and the previous conjugate direction.
 * A small generalized eigenvalue problem is solved in this subspace to update each band.
 * Optionally supports a blocked variant where bands are grouped and a single larger subspace
 * eigenvalue problem is solved per block.
 *
 * @tparam T The floating-point type used for calculations (default: std::complex<double>).
 * @tparam Device The device used for calculations (e.g., DEVICE_CPU or DEVICE_GPU).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoPPCG
{
  private:
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;

  public:
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /**
     * @brief Constructor for DiagoPPCG class.
     *
     * @param precondition_in Pointer to the preconditioner array with [dim: n_basis].
     */
    explicit DiagoPPCG(const Real* precondition_in);

    /**
     * @brief Destructor — frees all device and host allocations.
     */
    ~DiagoPPCG();

    /**
     * @brief Initialize the class before diagonalization.
     *
     * This function allocates all the related variables, such as hpsi, w, p, etc.,
     * before the diag call.
     *
     * @param nband The number of bands of all processes.
     * @param nband_l The number of bands of current process.
     * @param nbasis The number of basis functions. Leading dimension of psi.
     * @param ndim The number of valid dimension of psi.
     */
    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    /**
     * @brief Diagonalize the Hamiltonian using the PPCG method.
     *
     * @param hpsi_func A function computing the product of the Hamiltonian matrix H
     * and a wavefunction blockvector X.
     * @param psi_in Pointer to input wavefunction psi matrix with [dim: n_basis x n_band, column major].
     * @param eigenvalue_in Pointer to the eigen array with [dim: n_band].
     * @param ethr_band Convergence threshold for each band.
     * @return The number of iterations taken.
     */
    int diag(const HPsiFunc& hpsi_func,
             T* psi_in,
             Real* eigenvalue_in,
             const std::vector<double>& ethr_band);

  private:
    /// the number of bands of all processes
    int n_band = 0;
    /// the number of bands of current process
    int n_band_l = 0;
    /// the number of cols of the input psi
    int n_basis = 0;
    /// valid dimension of psi
    int n_dim = 0;
    /// number of extra bands for convergence acceleration (n_work = n_band_l + n_extra)
    int n_extra = 0;
    /// total working bands: n_band_l + n_extra
    int n_work = 0;

    /// Pointer to the preconditioner array (does not own memory).
    const Real* precondition = nullptr;
    /// Device-side copy of the preconditioner (GPU only).
    Real* d_precondition = nullptr;

    /// Device context
    Device* ctx = {};
    base_device::AbacusDevice_t device = {};

    // ---- device-side working arrays (n_work × n_basis) ----
    T* hpsi = nullptr;      ///< H|psi>
    T* w = nullptr;         ///< preconditioned residual W = -K^{-1} R
    T* hw = nullptr;        ///< H|w>
    T* p = nullptr;         ///< conjugate directions
    T* hp = nullptr;        ///< H|p>
    T* p_new = nullptr;     ///< updated p for next iteration
    T* hp_new = nullptr;    ///< H|p_new>
    T* hpsi_new = nullptr;  ///< updated H|psi>
    T* work = nullptr;      ///< workspace for rotations / intermediates

    /// device-side eigenvalues / errors [dim: n_work]
    Real* d_eigen = nullptr;
    Real* d_err = nullptr;

    /// host-side mirrors (for MPI reduce, convergence check, output)
    Real* h_eigen = nullptr;
    Real* h_err = nullptr;

    // ---- control state (host only, small) ----
    std::vector<char> is_locked;       ///< convergence lock flags
    std::vector<int> converge_count;   ///< consecutive convergence counters
    std::vector<int> block_sizes;      ///< block sizes for blocked variant

  public:
    /**
     * @brief Set the block sizes for the blocked PPCG variant.
     *
     * When set, update_vectors_from_ppcg_subspace switches from per-band (2D/3D)
     * subspace diagonalization to a blocked approach where bands within each block
     * are solved jointly in a 3k_i-dimensional subspace.
     *
     * @param sizes Vector of block sizes. An empty vector disables the blocked variant.
     */
    void set_block_sizes(const std::vector<int>& sizes)
    {
        this->block_sizes = sizes;
    }
    /**
     * @brief Set the number of extra bands used for convergence acceleration.
     *
     * Extra bands (n_extra) are added to the working set beyond n_band_l.
     * They participate in orthonormalization and subspace construction,
     * helping to accelerate convergence of the physical bands.
     *
     * @param n Number of extra bands.
     */
    void set_n_extra(const int n)
    {
        this->n_extra = n;
    }

  private:
    // ------------------------------------------------------------------
    //  memory-operation aliases
    // ------------------------------------------------------------------
    using resmem_op   = base_device::memory::resize_memory_op<T, Device>;
    using delmem_op   = base_device::memory::delete_memory_op<T, Device>;
    using setmem_op   = base_device::memory::set_memory_op<T, Device>;
    using syncmem_op  = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_d2h = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_h2d = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;

    using resmem_real_op    = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_real_op    = base_device::memory::delete_memory_op<Real, Device>;
    using setmem_real_op    = base_device::memory::set_memory_op<Real, Device>;
    using syncmem_real_h2d  = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_real_d2h  = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;

    using resmem_real_h = base_device::memory::resize_memory_op<Real, base_device::DEVICE_CPU>;
    using delmem_real_h = base_device::memory::delete_memory_op<Real, base_device::DEVICE_CPU>;

    // ------------------------------------------------------------------
    //  basic vector operations (operate on n_dim elements)
    // ------------------------------------------------------------------

    /// lhs^H * rhs  with MPI reduction.
    T inner_product(const T* lhs, const T* rhs) const;
    /// L2 norm.
    Real vector_norm(const T* vec) const;
    /// vec *= alpha, pad region zeroed.
    void scale_vector(T* vec, const Real alpha) const;
    /// y += alpha * x.
    void axpy_vector(T* y, const T* x, const T alpha) const;
    /// Copy n_basis elements.
    void copy_vector(T* dst, const T* src) const;
    /// Zero-fill n_basis elements.
    void zero_vector(T* vec) const;

    // ------------------------------------------------------------------
    //  higher-level operations
    // ------------------------------------------------------------------

    /// MPI-parallel convergence check.
    bool test_error(const std::vector<double>& ethr_band) const;
    /// hpsi_out = H |psi_in>
    void calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, T* hpsi_out) const;
    /// Modified Gram-Schmidt orthonormalization.
    void modified_gram_schmidt(T* psi_in, T* hpsi_in) const;
    /// Cholesky-based orthonormalization (more robust).
    void orth_cholesky(T* psi_in, T* hpsi_in);
    /// Check || <psi|psi> - I ||_F < 1e-1.
    bool check_orthonormality(T* psi_in) const;
    /// block_out = block * coeff  (gemm).
    void rotate_block(T* block, const T* coeff, T* workspace) const;
    /// Rayleigh-Ritz: Hsub = psi^H hpsi, diagonalize, rotate.
    void rayleigh_ritz(T* psi_in, T* hpsi_in);
    /// Compute preconditioned residuals and Rayleigh quotients.
    void calc_preconditioned_residual(T* psi_in);
    /// v_i -= sum_j <x_j|v_i> x_j  for each v in block.
    void project_to_orthogonal_complement(T* psi_in, T* block) const;
    /// Solve 2×2 / 3×3 generalized eigenproblem.
    bool solve_small_problem(const int active_dim, T* hsmall, T* ssmall, T* coeff, Real* eval) const;
    /// Per-band PPCG subspace update.
    void update_vectors_from_ppcg_subspace(T* psi_in);
    /// Block-diagonal PPCG subspace update.
    void update_vectors_blocked(T* psi_in);
};

} // namespace hsolver

#endif

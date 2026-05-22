#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/macros.h"
#include "source_base/module_device/types.h"

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
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, std::complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /**
     * @brief Constructor for DiagoPPCG class.
     *
     * @param precondition_in Pointer to the preconditioner array with [dim: n_basis].
     */
    explicit DiagoPPCG(const Real* precondition_in);

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
     * On GPU devices, falls back to DiagoBPCG. On CPU, runs the PPCG iteration:
     * each step computes the preconditioned residual, updates band locking,
     * constructs a per-band (or per-block) subspace, solves a small generalized
     * eigenvalue problem, and periodically re-orthonormalizes via Cholesky.
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
    /// @note prec[dim: n_basis]
    const Real* precondition = nullptr;

    /// H|psi> matrix [dim: n_basis x n_work, column major]
    std::vector<T> hpsi;
    /// Preconditioned residual vectors W = -K * R [dim: n_basis x n_work, column major]
    std::vector<T> w;
    /// H|w> matrix [dim: n_basis x n_work, column major]
    std::vector<T> hw;
    /// Conjugate direction vectors P [dim: n_basis x n_work, column major]
    std::vector<T> p;
    /// H|p> matrix [dim: n_basis x n_work, column major]
    std::vector<T> hp;
    /// Updated conjugate direction vectors for next iteration
    std::vector<T> p_new;
    /// H|p_new> matrix for next iteration
    std::vector<T> hp_new;
    /// Updated H|psi> matrix for next iteration
    std::vector<T> hpsi_new;
    /// Workspace buffer for vector rotations and intermediate results
    std::vector<T> work;
    /// Computed eigenvalues [dim: n_work]
    std::vector<Real> eigen;
    /// Residual norm for each band [dim: n_work]
    std::vector<Real> err;

    /// Convergence lock flag for each band [dim: n_work]
    std::vector<bool> is_locked;
    /// Consecutive convergence counter for each band [dim: n_work]
    std::vector<int> converge_count;

    /// Block sizes for the blocked PPCG variant; empty means per-band mode
    std::vector<int> block_sizes;

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
    /// @name Basic vector operations (operate on n_dim elements)
    /// @{

    /**
     * @brief Compute the inner product of two vectors: sum conj(lhs[i]) * rhs[i].
     * @note Includes MPI reduction across pool processes.
     */
    T inner_product(const T* lhs, const T* rhs) const;
    /// Compute the L2 norm of a vector.
    Real vector_norm(const T* vec) const;
    /// In-place scale a vector by a real scalar: vec *= alpha.
    void scale_vector(T* vec, const Real alpha) const;
    /// Compute y += alpha * x.
    void axpy_vector(T* y, const T* x, const T alpha) const;
    /// Copy n_basis elements from src to dst.
    void copy_vector(T* dst, const T* src) const;
    /// Zero-fill n_basis elements of vec.
    void zero_vector(T* vec) const;

    /// @}

    /**
     * @brief Check whether all bands satisfy the convergence threshold.
     *
     * @param ethr_band Convergence threshold for each band [dim: n_band].
     * @return true if any band (across all MPI ranks) is not converged, false if all converged.
     */
    bool test_error(const std::vector<double>& ethr_band) const;

    /**
     * @brief Apply the H operator to psi and obtain the hpsi matrix.
     *
     * @note hpsi_out = H|psi_in>
     *
     * @param hpsi_func A function computing the product of the Hamiltonian matrix H
     * and a wavefunction blockvector X.
     * @param psi_in Input wavefunction [dim: n_basis x n_work, column major].
     * @param hpsi_out Output H|psi> matrix [dim: n_basis x n_work, column major].
     */
    void calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, std::vector<T>& hpsi_out) const;

    /**
     * @brief Orthonormalize psi and hpsi using Modified Gram-Schmidt.
     *
     * @note psi_in and hpsi_in are modified in-place, column by column.
     * Aborts if linear dependence is detected (norm <= 1e-14).
     */
    void modified_gram_schmidt(T* psi_in, std::vector<T>& hpsi_in) const;

    /**
     * @brief Orthonormalize psi and hpsi using Cholesky decomposition of the overlap matrix.
     *
     * Computes S = <psi|psi>, factorizes S = L * L^H, then rotates vectors by L^{-1}.
     * More numerically robust than Gram-Schmidt for large block sizes or near-linear-dependence.
     */
    void orth_cholesky(T* psi_in, std::vector<T>& hpsi_in);

    /**
     * @brief Verify orthonormality of the working vectors.
     *
     * @return true if the Frobenius norm of (S - I) < 1e-6, false otherwise.
     */
    bool check_orthonormality(T* psi_in) const;

    /**
     * @brief Rotate a block of vectors by a coefficient matrix: block_out = block * coeff.
     *
     * @param block Input/output block of vectors [dim: n_basis x n_work, column major].
     * @param coeff Rotation coefficient matrix [dim: n_work x n_work, column major].
     * @param workspace Workspace buffer [dim: n_basis x n_work, column major].
     */
    void rotate_block(T* block, const std::vector<T>& coeff, std::vector<T>& workspace) const;

    /**
     * @brief Perform the Rayleigh-Ritz procedure.
     *
     * Builds the subspace Hamiltonian Hsub = <psi|H|psi>, diagonalizes it
     * via LAPACK zheevd, and rotates psi and hpsi by the eigenvectors.
     * On exit, eigenvalues are sorted ascending.
     */
    void rayleigh_ritz(T* psi_in, std::vector<T>& hpsi_in);

    /**
     * @brief Compute the preconditioned residual and eigenvalue for each band.
     *
     * For each non-locked band, computes:
     *   1. lambda_i = <x_i | H | x_i> (Rayleigh quotient as eigenvalue estimate)
     *   2. R_i = H x_i - lambda_i x_i (residual)
     *   3. w_i = -K^{-1} R_i (preconditioned residual)
     *
     * The residual norm is stored in err[ib] and reduced across MPI processes.
     * Locked bands have their w vector zeroed.
     */
    void calc_preconditioned_residual(T* psi_in);

    /**
     * @brief Project block vectors onto the orthogonal complement of the current subspace.
     *
     * For each vector v in block, subtracts its projection onto all current psi vectors:
     * v_i = v_i - sum_j <x_j | v_i> * x_j
     */
    void project_to_orthogonal_complement(T* psi_in, std::vector<T>& block) const;

    /**
     * @brief Solve a small generalized eigenvalue problem H * C = lambda * S * C.
     *
     * Uses LAPACK zhegvd. Falls back to the first basis vector on failure.
     *
     * @param active_dim Dimension of the small problem (2 or 3).
     * @param hsmall Subspace Hamiltonian matrix [dim: active_dim x active_dim, column major].
     * @param ssmall Subspace overlap matrix [dim: active_dim x active_dim, column major].
     * @param coeff Output eigenvector coefficients [dim: active_dim x active_dim, column major].
     * @param eval Output eigenvalues [dim: active_dim].
     * @return true on success, false if the generalized eigenproblem failed.
     */
    bool solve_small_problem(const int active_dim, T* hsmall, T* ssmall, T* coeff, Real* eval) const;

    /**
     * @brief Update psi, hpsi, p, hp from the per-band PPCG subspace.
     *
     * For each non-locked band, constructs a 2D or 3D subspace from {x_i, w_i, p_i},
     * solves a small generalized eigenvalue problem, and updates the working vectors
     * using the lowest eigenvector's coefficients.
     *
     * If block_sizes is set, delegates to update_vectors_blocked instead.
     */
    void update_vectors_from_ppcg_subspace(T* psi_in);

    /**
     * @brief Block-diagonal variant of the PPCG subspace update.
     *
     * Groups bands into blocks. For each block of size k_i, constructs a
     * 3k_i-dimensional subspace from {X_block, W_block, P_block}, solves
     * the generalized eigenvalue problem, and updates all bands in the block
     * simultaneously using the first k_i eigenvectors.
     */
    void update_vectors_blocked(T* psi_in);
};

} // namespace hsolver

#endif

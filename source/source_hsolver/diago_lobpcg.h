#ifndef DIAGO_LOBPCG_H_
#define DIAGO_LOBPCG_H_

#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h"
#include "source_hsolver/module_linalg/linear_algebra.h"

#include <complex>
#include <string>
#include <vector>

namespace hsolver
{

namespace dense
{
class ParallelDenseEigensolver;
class GpuDenseEigensolver;
struct SerialGvWorkspace;
} // namespace dense

enum class LobpcgLocking
{
    Hard,
    Soft
};

/** Valid rows per spinor component; the component stride is ndim / npol. */
struct LobpcgBasisLayout
{
    int npw;  ///< Local physical plane-wave rows in each spinor component.
    int npol; ///< Spinor components sharing the column storage stride ndim / npol.
};

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoLOBPCG
{
    friend class LobpcgProjectionTest;

  private:
    using Real = typename GetTypeReal<T>::type; ///< Real precision of the stored wavefunctions T.
    using Wide = std::complex<double>; ///< Accumulation and projected eigensolver precision.
    using CtDevice = typename ct::PsiToContainer<Device>::type;

  public:
    // Internal block: nmax >= nband > 0. Vectors and preconditioner have ndim entries.
    // Layout has npol > 0 dividing ndim and 0 <= npw <= ndim / npol.
    DiagoLOBPCG(const Real* precondition,
        const int nband,
        const int ndim,
        const int nmax,
        const LobpcgBasisLayout& layout,
        const diag_comm_info& comm,
        const int para_rr_method,
        const int para_rr_block_size,
        const bool use_overlap,
        const LobpcgLocking locking,
        dense::ParallelDenseEigensolver& rr_solver,
        dense::GpuDenseEigensolver& gpu_rr_solver);
    ~DiagoLOBPCG();

    // Start from a full-rank block of n_band_ vectors with ld_psi >= n_dim_.
    // max_iter and tolerance are positive; false indicates iteration-budget exhaustion.
    bool diag(const HSOperator<T, Device>& op,
        T* psi_in,
        const int ld_psi,
        Real* eigenvalue_in,
        const double tolerance,
        const int max_iter,
        const std::vector<double>& ethr_band);

    int iterations() const { return iterations_; }

  private:
    // Non-owning Device views; sx is null for the identity metric.
    struct VectorImages
    {
        T* x;  ///< Borrowed column block on Device.
        T* hx; ///< Matching H images; null when a transformation does not need them.
        T* sx; ///< Matching S images; null denotes S = I.
    };

    using SerialGvWorkspace = dense::SerialGvWorkspace;

    // One diag lifetime: root-owned K/G of retained [X,P] and reusable scratch.
    // All ranks need columns when assembling the new W panel.
    struct ProjectionCache
    {
        int columns = 0; ///< Width of the retained [X,P] block, also known on nonroot ranks.
        /// Root-owned projections of [X,P]; empty on other ranks.
        ct::Tensor k{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
        ct::Tensor g{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
        /// Host scratch for Q*C_P and the packed coefficient matrix [C_X,C_P].
        ct::Tensor h_product{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
        ct::Tensor structural_coefficients{ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {0}};
        /// Device scratch for the same retained-block products; host K/G remain authoritative.
        ct::Tensor device_matrix;
        ct::Tensor device_coefficients;
        ct::Tensor device_product;
        ct::Tensor device_result;
        ~ProjectionCache();
    };

    int n_band_ = 0; ///< Number of requested eigenpairs returned to the caller.
    int n_dim_ = 0; ///< Local column leading dimension, including spinor padding.
    int n_max_ = 0; ///< Internal block width: requested bands plus dimension-limited buffer bands.
    int iterations_ = 0; ///< Main-loop RR iterations in the last diag call; excludes initial RR.
    std::string error_message_; ///< Failure reason supplied by RR or recovery.

    const LobpcgBasisLayout layout_; ///< Physical rows and spinor count within each stored column.
    const diag_comm_info comm_; ///< Communicator/rank for distributed rows of one eigensolve.
    const int para_rr_method_ = 0; ///< RR backend: local (0) or ScaLAPACK (2).
    const int para_rr_bs_ = 0; ///< Distributed RR block size; zero selects automatically.
    const bool use_overlap_; ///< Store and apply explicit S images instead of S = I.
    const LobpcgLocking locking_; ///< Hard freezes the converged prefix; soft retains all X in RR.
    dense::ParallelDenseEigensolver& rr_solver_; ///< Borrowed, reusable distributed CPU RR resources.
    dense::GpuDenseEigensolver& gpu_rr_solver_; ///< Borrowed, reusable CUDA RR resources.

    ct::DataType r_type_ = ct::DataType::DT_INVALID; ///< Tensor element type corresponding to Real.
    ct::DataType t_type_ = ct::DataType::DT_INVALID; ///< Tensor element type corresponding to T.
    ct::DeviceType device_type_ = ct::DeviceType::UnKnown; ///< Wavefunction/scratch storage device.

    ct::Tensor prec_; ///< Device copy of the caller's diagonal preconditioner input, n_dim_ entries.
    // Large wavefunction blocks reside on Device; small spectra, RR matrices
    // and lock state reside on CPU on both backends.
    ct::Tensor eig_; ///< Host Ritz values for all n_max_ columns, including locked columns.
    ct::Tensor eig_device_; ///< Active Ritz values on Device for preconditioning.
    ct::Tensor evec_; ///< Device rotation scratch, n_dim_ by n_max_.
    ct::Tensor space_; ///< Device basis [X,P,W], capacity n_dim_ by 3*n_max_.
    ct::Tensor hspace_; ///< H images matching the columns of space_.
    ct::Tensor sspace_; ///< Explicit S images matching space_; unused for S = I.
    ct::Tensor ritz_coefficients_; ///< Host T coefficients C, capacity 3*n_max_ by n_max_.
    ct::Tensor e_red_; ///< Host eigenvalue output of the current reduced problem.
    ct::Tensor x_new_; ///< Device updated Ritz vectors; hard-locked prefix stays unchanged.
    ct::Tensor hx_new_; ///< Updated H images needed by residuals and the next trial basis.
    ct::Tensor sx_new_; ///< Updated S images; unused for S = I.
    ct::Tensor residual_; ///< Device residuals HX - SX*Lambda, active columns at their band offsets.
    ct::Tensor u_p_; ///< Host T history coefficients C_P, compacted/scaled together with P.

    const T zero_ = static_cast<T>(0.0);

    using syncmem_complex_2d_op = base_device::memory::synchronize_memory_2d_op<T, Device, Device>;
    using syncmem_var_h2d_op = ct::kernels::synchronize_memory<Real, CtDevice, ct::DEVICE_CPU>;
    using setmem_complex_op = ct::kernels::set_memory<T, CtDevice>;
    using syncmem_complex_h2d_op = ct::kernels::synchronize_memory<T, CtDevice, ct::DEVICE_CPU>;
    using copy_op = ct::kernels::blas_copy<T, CtDevice>;
    using axpy_op = ct::kernels::blas_axpy<T, CtDevice>;
    using nrm2_op = ct::kernels::blas_nrm2<T, CtDevice>;

    copy_op copy_op_;
    axpy_op axpy_op_;
    nrm2_op nrm2_op_;
    precondition_op<T, Device> precondition_op_;

    LinearAlgebra<T, Device> algebra_; ///< Distributed products and reusable conversion/rotation scratch.

    // Main algorithm operations.
    void initialize(T* psi_in, const int ld_psi);
    void form_residual(const VectorImages& vectors, const int count, const Real* eigenvalues, T* residual);
    int build_search(const HSOperator<T, Device>& op, const int locked, const int count, const int first_w);
    bool lock_converged(ct::Tensor& eig_prev, const ct::Tensor& tol_band, int& locked);
    // The caller supplies a nonempty active suffix; rank compaction follows rotation.
    int build_direction(const VectorImages& basis,
        const int dimension,
        const int q,
        const int active,
        const ct::Tensor& x_coefficients,
        const ct::Tensor& g);
    void write_results(T* psi_in, const int ld_psi, Real* eigenvalue_in);

    // Shared vector transformations and storage.
    void rotate_images(const VectorImages& source,
        const int source_columns,
        const int target_columns,
        const T* coefficients,
        const VectorImages& target,
        T* scratch);
    int normalize_columns(const VectorImages& block, const int columns, T* coefficients, const int coefficient_rows);
    const T* prepare_rotation_coefficients(const ct::Tensor& coefficients);
    VectorImages basis_images(const int first);
    void project_padding(T* x, const int columns) const;
    template <typename Value>
    void broadcast(Value* data, const int count) const;
    void allocate_workspace();

    // Projected matrices, generalized RR and retained directions.
    void form_projected_matrices(
        const VectorImages& basis, const int columns, const ProjectionCache* retained, ct::Tensor& k, ct::Tensor& g);
    bool solve_projected(const ct::Tensor& k,
        const ct::Tensor& g,
        const int dimension,
        const int wanted,
        SerialGvWorkspace& serial_workspace);
    void retain_pair(const ct::Tensor& k,
        const ct::Tensor& g,
        const int rows,
        const ct::Tensor& x_coefficients,
        const int first_x,
        const int x_columns,
        const int p_columns,
        const Real* ritz_values,
        const int zero_prefix,
        ProjectionCache& cache);

    // Physical-metric and rank recovery after a failed projected solve.
    bool recover_projected(const HSOperator<T, Device>& op,
        const VectorImages& basis,
        int& working_columns,
        const int wanted,
        ct::Tensor& projected,
        ct::Tensor& gram,
        SerialGvWorkspace& workspace);
    int repair_basis(const HSOperator<T, Device>& op, const int added);
    int rank_orthonormalize(const int m, T* x, T* sx);
};

} // namespace hsolver

#endif // DIAGO_LOBPCG_H_

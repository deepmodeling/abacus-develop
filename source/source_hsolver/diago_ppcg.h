#ifndef DIAGO_PPCG_H
#define DIAGO_PPCG_H

#include "ppcg_workspace.h"
#include "source_base/module_device/types.h"

#include <complex>
#include <iosfwd>
#include <type_traits>
#include <vector>

namespace hsolver
{

template <typename T, typename Device>
class HSOperator;

struct PPCGTraceContext
{
    int scf_iteration;
    int kpoint;
};

// -----------------------------------------------------------------------------
// DiagoPPCG: Projection Preconditioned solver
// -----------------------------------------------------------------------------
//
// Implements the block-subspace diagonalization strategy used by
// ks_solver=ppcg. The solver updates groups of active Ritz pairs in a bounded
// local subspace with configurable global Rayleigh-Ritz steps. Operator images
// are cached, with independent periodic and cancellation-triggered refreshes.
// -----------------------------------------------------------------------------

namespace base_device = ::base_device;

// Check workspace sizes before allocating. Dense kernels and pool
// reductions use int element counts.
bool ppcg_workspace_dimensions_valid(
    int ld_psi,
    int nband,
    int dim,
    int sbsize);

template <typename T, typename Device>
class DiagoPPCG
{
  public:
    // -------------------------------------------------------------------------
    // Type aliases
    // -------------------------------------------------------------------------
    using Real = typename std::conditional<std::is_same<T, std::complex<double>>::value, double, float>::type;

    using Workspace = typename PPCGWorkspaceType<T, Device>::type;

    // -------------------------------------------------------------------------
    // Constructor
    // rr_step is the explicit Ritz interval; PW INPUT defaults are supplied by
    // the caller. The internal H/S cache-refresh interval is independent of it.
    // -------------------------------------------------------------------------
    DiagoPPCG(
        const Real &diag_thr,
        const int &diag_iter_max,
        const int &sbsize,
        const int &rr_step);

    // -------------------------------------------------------------------------
    // Main entry point
    //
    // Returns average number of subspace iterations per band.
    // Optional ABACUS_PPCG_RESIDUAL_TRACE CSV: initial_rr and rayleigh_ritz
    // report global Ritz residuals. residual_refinement and final report the
    // residuals of the current eigenpairs.
    // block_sweep rows instead contain the current column Rayleigh quotients;
    // exclude those rows when plotting ordered eigenvalue errors.
    // -------------------------------------------------------------------------
    double diag(
        const HSOperator<T, Device> &op,
        int ld_psi,
        int nband,
        int dim,
        T *psi_in,
        Real *eigenvalue_in,
        const std::vector<double> &ethr_band,
        const Real *prec,
        std::ostream &log,
        const PPCGTraceContext &trace_context);

    bool converged() const { return converged_; }

    int active_band_count() const { return active_band_count_; }

    bool iteration_limit_reached() const { return iteration_limit_reached_; }

    Real max_residual() const { return max_residual_; }

    int first_active_band() const { return first_active_band_; }

    Real first_active_threshold() const { return first_active_threshold_; }

    int rayleigh_ritz_count() const { return rayleigh_ritz_count_; }

    int root_recovery_count() const { return root_recovery_count_; }

  private:
    friend class DiagoPPCGProjectionTest;
    friend class DiagoPPCGRefinementTest;
    // -------------------------------------------------------------------------
    // Data members
    // -------------------------------------------------------------------------
    int maxiter_;
    int sbsize_;
    int rr_step_;
    Real diag_thr_;

    // Problem dimensions (set in diag())
    int ld_psi_ = 0;
    int n_band_ = 0;
    int n_dim_ = 0;
    bool converged_ = false;
    int active_band_count_ = 0;
    bool iteration_limit_reached_ = false;
    Real max_residual_ = Real(0);
    int first_active_band_ = -1;
    Real first_active_threshold_ = Real(0);
    int rayleigh_ritz_count_ = 0;
    int root_recovery_count_ = 0;

    // Working storage (column-major: ld_psi_ rows, n_band_ columns).
    Workspace hpsi_;
    Workspace spsi_;
    Workspace w_;  // residual / preconditioned residual
    Workspace sw_; // S * w
    Workspace hw_; // H * w
    Workspace p_;  // previous search direction (LOBPCG)
    Workspace sp_; // S * p
    Workspace hp_; // H * p
    std::vector<T> rr_hsub_;
    std::vector<T> rr_ssub_;
    std::vector<Real> rr_eval_;
    typename PPCGWorkspaceType<Real, Device>::type device_prec_;

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------
    static int idx(int row, int col, int ld) { return row + col * ld; }

    void validate_input(
        const T *psi_in,
        const Real *eigenvalue_in,
        const std::vector<double> &ethr_band,
        const Real *prec) const;

    // Apply the shared H/S operator to a block.
    void apply_h(
        const HSOperator<T, Device> &op,
        T *psi_in,
        T *hpsi_out,
        int ncol) const;
    void apply_s(
        const HSOperator<T, Device> &op,
        T *psi_in,
        T *spsi_out,
        int ncol) const;

    // Inner product <x|y> (real part only).
    Real gamma_dot(const T *x, const T *y) const;

    // Gram matrix: out[i, j] = <a_i | b_j>.
    void gram(
        const T *mat_a,
        const T *mat_b,
        int ncol_a,
        int ncol_b,
        std::vector<T> &out,
        int ld_out) const;

    // Gather / scatter columns.
    void copy_cols(const T *src, const std::vector<int> &cols, Workspace &dst) const;
    void scatter_cols(T *dst, const std::vector<int> &cols, const Workspace &src) const;

    // Project x onto vectors orthogonal to the S-orthonormal basis.
    void project_against(
        const T *basis,
        const T *sbasis,
        const T *hbasis,
        const std::vector<int> &basis_cols,
        Workspace &x,
        Workspace &sx,
        Workspace *hx,
        const std::vector<int> &x_cols) const;

    void project_search_directions(
        const HSOperator<T, Device> &op,
        const T *psi,
        const std::vector<int> &all_cols,
        const std::vector<int> &active_cols);

    // x[c] /= max(prec, eps)  for each active column c.
    void divide_by_preconditioner(const std::vector<int> &active_cols, const Real *prec, Workspace &x) const;

    // -------------------------------------------------------------------------
    // Block-subspace strategy helpers (File 1 style)
    // -------------------------------------------------------------------------
    struct SmallSubspace
    {
        std::vector<T> k;       // K matrix (projected H)
        std::vector<T> m;       // M matrix (projected S)
        std::vector<Real> eval; // eigenvalues
        Workspace psi_l;
        Workspace spsi_l;
        Workspace hpsi_l;
        Workspace w_l;
        Workspace sw_l;
        Workspace hw_l;
        Workspace p_l;
        Workspace sp_l;
        Workspace hp_l;
        Workspace basis;
        Workspace hbasis;
        Workspace sbasis;
        std::vector<T> coeff_state;
        std::vector<T> coeff_p;
        Workspace psi_new;
        Workspace spsi_new;
        Workspace hpsi_new;
        Workspace p_new;
        Workspace sp_new;
        Workspace hp_new;
    };

    struct OperatorBlock
    {
        Workspace vectors;
        Workspace hvectors;
        Workspace svectors;
    };

    void scale_to_unit_snorm(
        Workspace &x,
        Workspace &sx,
        Workspace &hx,
        int ncols) const;

    void hermitize_projected(std::vector<T> &matrix, int dim) const;

    void insert_startup_gram_block(
        const T *left,
        const T *right,
        int row_count,
        int column_count,
        int row_offset,
        int column_offset,
        int matrix_dim,
        std::vector<T> &matrix,
        std::vector<T> &workspace) const;

    void build_startup_global_subspace(
        const T *psi,
        Workspace &w_active,
        Workspace &sw_active,
        Workspace &hw_active,
        int active_count,
        SmallSubspace &subspace);

    void update_startup_global_subspace(
        T *psi,
        const Workspace &w_active,
        const Workspace &sw_active,
        const Workspace &hw_active,
        int active_count,
        const SmallSubspace &subspace);

    void lock_epairs(const Workspace &residual, const std::vector<double> &ethr_band, std::vector<int> &active_cols);

    void build_small_subspace(
        const T *psi,
        const std::vector<int> &cols,
        int nblk,
        SmallSubspace &subspace) const;

    // Returns the number of valid Ritz vectors, including discarded states.
    // An identity fallback reports zero rather than pretending it is a solve.
    int solve_small_generalized(const int dim, const int nstates, SmallSubspace &subspace) const;

    void collect_discarded_ritz(
        const SmallSubspace &subspace,
        int kept_count,
        int valid_count,
        OperatorBlock &candidates) const;

    void project_candidate_complement(const HSOperator<T, Device> &op, const T *psi, OperatorBlock &candidates) const;

    bool recover_competing_ritz(
        const HSOperator<T, Device> &op,
        T *psi,
        Real *eigenvalue,
        bool ordered_ritz,
        std::vector<int> &active_cols,
        const std::vector<double> &ethr_band,
        OperatorBlock &candidates);

    void update_one_block(
        const HSOperator<T, Device> &op,
        T *psi,
        const std::vector<int> &cols,
        int l,
        int nblk,
        SmallSubspace &subspace);

    bool rayleigh_ritz(
        T *psi,
        Real *eigenvalue,
        std::vector<int> &active_cols,
        const std::vector<double> &ethr_band);

    bool refine_stalled_residual(
        const HSOperator<T, Device> &op,
        T *psi,
        Real *eigenvalue,
        const std::vector<double> &ethr_band,
        const std::vector<int> &all_cols,
        int band_hint);

    void compute_eigenpair_residual(const Real *eigenvalue);

    bool orthonormalize(T *psi);

    // Between Ritz solves the columns need not be individual eigenvectors.
    // Use the block gradient, and request a Ritz solve before declaring them
    // converged or changing the locked prefix.
    bool update_block_gradient(const T *psi, Real *eigenvalue, const std::vector<double> &ethr_band);
};

} // namespace hsolver

#endif // DIAGO_PPCG_H

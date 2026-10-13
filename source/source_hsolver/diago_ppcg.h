#ifndef DIAGO_PPCG_H
#define DIAGO_PPCG_H

#include "source_base/module_device/types.h"

#include <complex>
#include <functional>
#include <type_traits>
#include <vector>

namespace hsolver {

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

template <typename T, typename Device>
class DiagoPPCG
{
public:
    // -------------------------------------------------------------------------
    // Type aliases
    // -------------------------------------------------------------------------
    using Real = typename std::conditional<
        std::is_same<T, std::complex<double>>::value, double,
        float>::type;
    using HPsiFunc = std::function<void(T*, T*, int, int)>;
    using SPsiFunc = std::function<void(T*, T*, int, int)>;

    // -------------------------------------------------------------------------
    // Constructor
    // rr_step is the explicit Ritz interval; PW INPUT defaults are supplied by
    // the caller. The internal H/S cache-refresh interval is independent of it.
    // -------------------------------------------------------------------------
    DiagoPPCG(const Real& diag_thr,
              const int& diag_iter_max,
              const int& sbsize,
              const int& rr_step);

    // -------------------------------------------------------------------------
    // Main entry point
    //
    // Returns average number of subspace iterations per band.
    // Optional ABACUS_PPCG_RESIDUAL_TRACE CSV: initial_rr, rayleigh_ritz and
    // final rows contain ordered Ritz values after a successful Ritz solve.
    // block_sweep rows instead contain the current column Rayleigh quotients;
    // exclude those rows when plotting ordered eigenvalue errors.
    // -------------------------------------------------------------------------
    double diag(const HPsiFunc& hpsi_func,
                const SPsiFunc& spsi_func,
                int ld_psi,
                int nband,
                int dim,
                T* psi_in,
                Real* eigenvalue_in,
                const std::vector<double>& ethr_band,
                const Real* prec);

    bool converged() const
    {
        return converged_;
    }

    int active_band_count() const
    {
        return active_band_count_;
    }

    bool iteration_limit_reached() const
    {
        return iteration_limit_reached_;
    }

    Real max_residual() const
    {
        return max_residual_;
    }

    int first_active_band() const
    {
        return first_active_band_;
    }

    Real first_active_threshold() const
    {
        return first_active_threshold_;
    }

    int rayleigh_ritz_count() const
    {
        return rayleigh_ritz_count_;
    }

    int root_recovery_count() const
    {
        return root_recovery_count_;
    }

private:
    friend class DiagoPPCGProjectionTest;
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

    // Cached S-operator (null if identity).
    SPsiFunc spsi_func_;

    // Working storage (column-major: ld_psi_ rows, n_band_ columns).
    std::vector<T> hpsi_;
    std::vector<T> spsi_;
    std::vector<T> w_;       // residual / preconditioned residual
    std::vector<T> sw_;      // S * w
    std::vector<T> hw_;      // H * w
    std::vector<T> p_;       // previous search direction (LOBPCG)
    std::vector<T> sp_;      // S * p
    std::vector<T> hp_;      // H * p
    std::vector<T> rr_hsub_;
    std::vector<T> rr_ssub_;
    std::vector<Real> rr_eval_;

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------
    static inline int idx(int row, int col, int ld)
    {
        return row + col * ld;
    }

    void validate_input(const HPsiFunc& hpsi_func,
                        const T* psi_in, const Real* eigenvalue_in,
                        const std::vector<double>& ethr_band,
                        const Real* prec) const;

    // S-application (identity fallback if spsi_func is null).
    void apply_h(const HPsiFunc& hpsi_func, T* psi_in, T* hpsi_out,
                 int ncol) const;
    void apply_s(const SPsiFunc& spsi_func, T* psi_in, T* spsi_out,
                 int ncol) const;
    void apply_s_current(T* psi_in, T* spsi_out, int ncol) const;

    // Inner product <x|y> (real part only).
    Real gamma_dot(const T* x, const T* y) const;

    // Gram matrix: out[i, j] = <a_i | b_j>.
    void gram(const T* mat_a, const T* mat_b,
              int ncol_a, int ncol_b,
              std::vector<T>& out, int ld_out) const;

    // Gather / scatter columns.
    void copy_cols(const T* src, const std::vector<int>& cols,
                   std::vector<T>& dst) const;
    void scatter_cols(T* dst, const std::vector<int>& cols,
                      const std::vector<T>& src) const;

    // Project x onto vectors orthogonal to the S-orthonormal basis.
    void project_against(const T* basis, const T* sbasis, const T* hbasis,
                         const std::vector<int>& basis_cols,
                         std::vector<T>& x, std::vector<T>& sx, std::vector<T>* hx,
                         const std::vector<int>& x_cols) const;

    void project_search_directions(const HPsiFunc& hpsi_func, const T* psi,
                                    const std::vector<int>& all_cols,
                                    const std::vector<int>& active_cols);

    // x[c] /= max(prec, eps)  for each active column c.
    void divide_by_preconditioner(const std::vector<int>& active_cols,
                                  const Real* prec,
                                  std::vector<T>& x) const;

    // -------------------------------------------------------------------------
    // Block-subspace strategy helpers (File 1 style)
    // -------------------------------------------------------------------------
    struct SmallSubspace
    {
        std::vector<T> k;        // K matrix (projected H)
        std::vector<T> m;        // M matrix (projected S)
        std::vector<Real> eval;  // eigenvalues
        std::vector<T> psi_l;
        std::vector<T> spsi_l;
        std::vector<T> hpsi_l;
        std::vector<T> w_l;
        std::vector<T> sw_l;
        std::vector<T> hw_l;
        std::vector<T> p_l;
        std::vector<T> sp_l;
        std::vector<T> hp_l;
        std::vector<T> basis;
        std::vector<T> hbasis;
        std::vector<T> sbasis;
        std::vector<T> coeff_state;
        std::vector<T> coeff_p;
        std::vector<T> psi_new;
        std::vector<T> spsi_new;
        std::vector<T> hpsi_new;
        std::vector<T> p_new;
        std::vector<T> sp_new;
        std::vector<T> hp_new;
    };

    struct OperatorBlock
    {
        std::vector<T> vectors;
        std::vector<T> hvectors;
        std::vector<T> svectors;
    };

    void scale_to_unit_snorm(std::vector<T>& x,
                             std::vector<T>& sx,
                             std::vector<T>& hx,
                             int ncols) const;

    void hermitize_projected(std::vector<T>& matrix, int dim) const;

    void insert_startup_gram_block(const T* left,
                                   const T* right,
                                   int row_count,
                                   int column_count,
                                   int row_offset,
                                   int column_offset,
                                   int matrix_dim,
                                   std::vector<T>& matrix,
                                   std::vector<T>& workspace) const;

    void build_startup_global_subspace(const T* psi,
                                       std::vector<T>& w_active,
                                       std::vector<T>& sw_active,
                                       std::vector<T>& hw_active,
                                       int active_count,
                                       SmallSubspace& subspace);

    void update_startup_global_subspace(T* psi,
                                        const std::vector<T>& w_active,
                                        const std::vector<T>& sw_active,
                                        const std::vector<T>& hw_active,
                                        int active_count,
                                        const SmallSubspace& subspace);

    void lock_epairs(const std::vector<T>& residual,
                     const std::vector<double>& ethr_band,
                     std::vector<int>& active_cols);

    void build_small_subspace(const T* psi,
                              const std::vector<int>& cols,
                              int nblk,
                              SmallSubspace& subspace) const;

    // Returns the number of valid Ritz vectors, including discarded states.
    // An identity fallback reports zero rather than pretending it is a solve.
    int solve_small_generalized(const int dim,
                                const int nstates,
                                SmallSubspace& subspace) const;

    void collect_discarded_ritz(const SmallSubspace& subspace,
                                int kept_count, int valid_count,
                                OperatorBlock& candidates) const;

    void project_candidate_complement(const HPsiFunc& hpsi_func,
                                      const T* psi, OperatorBlock& candidates) const;

    bool recover_competing_ritz(const HPsiFunc& hpsi_func,
                                T* psi, Real* eigenvalue,
                                bool ordered_ritz,
                                std::vector<int>& active_cols,
                                const std::vector<double>& ethr_band,
                                OperatorBlock& candidates);

    void update_one_block(T* psi,
                          const std::vector<int>& cols,
                          int l,
                          int nblk,
                          SmallSubspace& subspace);

    bool rayleigh_ritz(T* psi, Real* eigenvalue,
                       std::vector<int>& active_cols,
                       const std::vector<double>& ethr_band);

    void compute_eigenpair_residual(const Real* eigenvalue);

    bool orthonormalize(T* psi);

    // Between Ritz solves the columns need not be individual eigenvectors.
    // Use the block gradient, and request a Ritz solve before declaring them
    // converged or changing the locked prefix.
    bool update_block_gradient(const T* psi, Real* eigenvalue,
                               const std::vector<double>& ethr_band);
};

} // namespace hsolver

#endif // DIAGO_PPCG_H

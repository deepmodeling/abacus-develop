#ifndef DIAGO_PPCG_H
#define DIAGO_PPCG_H

#include "source_base/module_device/types.h"

#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <numeric>
#include <complex>
#include <type_traits>

namespace hsolver {

// -----------------------------------------------------------------------------
// DiagoPPCG: Projection Preconditioned Conjugate Gradient solver
// -----------------------------------------------------------------------------
//
// Supports two algorithmic strategies:
//   CONJUGATE_GRADIENT — band-by-band Polak-Ribiere CG with line minimization
//     (File 2 approach).
//   BLOCK_SUBSPACE — block subspace diagonalization (File 1 approach).
//
// BLOCK_SUBSPACE is the production path used by ks_solver=ppcg.
// CONJUGATE_GRADIENT is kept as an explicit fallback strategy.
// -----------------------------------------------------------------------------

enum class PpcgStrategy { BLOCK_SUBSPACE, CONJUGATE_GRADIENT };

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
    // -------------------------------------------------------------------------
    DiagoPPCG(const Real& diag_thr,
              const int& diag_iter_max,
              const int& sbsize,
              const int& rr_step,
              const bool gamma_g0_real,
              const PpcgStrategy strategy = PpcgStrategy::BLOCK_SUBSPACE);

    // -------------------------------------------------------------------------
    // Main entry point
    //
    // Returns average number of subspace iterations per band.
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

private:
    // -------------------------------------------------------------------------
    // Data members
    // -------------------------------------------------------------------------
    int maxiter_;
    int sbsize_;
    int rr_step_;
    Real diag_thr_;
    bool gamma_g0_real_;
    PpcgStrategy strategy_;

    // Problem dimensions (set in diag())
    int ld_psi_ = 0;
    int n_band_ = 0;
    int n_dim_ = 0;

    // Cached S-operator (null if identity).
    SPsiFunc spsi_func_;

    // Working storage (column-major: ld_psi_ rows, n_band_ columns).
    std::vector<T> hpsi_;
    std::vector<T> spsi_;
    std::vector<T> w_;       // residual / preconditioned residual
    std::vector<T> sw_;      // S * w
    std::vector<T> hw_;      // H * w
    std::vector<T> rr_psi_;  // Rayleigh-Ritz rotation workspace
    std::vector<T> rr_spsi_;
    std::vector<T> rr_hpsi_;
    std::vector<T> rr_hsub_;
    std::vector<T> rr_ssub_;
    std::vector<Real> rr_eval_;

    // Polak-Ribiere state (CONJUGATE_GRADIENT strategy)
    std::vector<T> z_old_;      // previous preconditioned residual
    std::vector<Real> beta_denom_;

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

    void force_g0_real(T* x, int ncol) const;

    // S-application (identity fallback if spsi_func is null).
    void apply_h(const HPsiFunc& hpsi_func, T* psi_in, T* hpsi_out,
                 int ncol) const;
    void apply_s(const SPsiFunc& spsi_func, T* psi_in, T* spsi_out,
                 int ncol) const;
    void apply_s_current(T* psi_in, T* spsi_out, int ncol) const;

    // Inner product <x|y> (real part only).
    Real gamma_dot(const T* x, const T* y) const;
    T complex_dot(const T* x, const T* y) const;

    // Gram matrix: out[i, j] = <a_i | b_j>.
    void gram(const T* a, const T* b,
              int ncol_a, int ncol_b,
              std::vector<T>& out, int ld_out) const;

    // Gather / scatter columns.
    void copy_cols(const T* src, const std::vector<int>& cols,
                   std::vector<T>& dst) const;
    void scatter_cols(T* dst, const std::vector<int>& cols,
                      const std::vector<T>& src) const;

    // Project x onto vectors orthogonal to the S-orthonormal basis.
    void project_against(const T* basis, const T* sbasis,
                         const std::vector<int>& basis_cols,
                         std::vector<T>& x, std::vector<T>& sx,
                         const std::vector<int>& x_cols) const;

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
        std::vector<T> basis;
        std::vector<T> hbasis;
        std::vector<T> sbasis;
        std::vector<T> coeff_state;
        std::vector<T> psi_new;
        std::vector<T> spsi_new;
        std::vector<T> hpsi_new;
    };

    void lock_epairs(const std::vector<T>& residual,
                     const std::vector<double>& ethr_band,
                     std::vector<int>& active_cols) const;

    void build_small_subspace(const T* psi,
                              const std::vector<int>& cols,
                              SmallSubspace& subspace) const;

    void solve_small_generalized(int dim, SmallSubspace& subspace) const;

    void update_one_block(T* psi,
                          const std::vector<int>& cols,
                          int l,
                          SmallSubspace& subspace);

    bool is_s_orthonormal(const T* psi, const T* spsi, int ncol) const;

    void s_gram_schmidt(T* psi, T* hpsi, T* spsi, int ncol) const;

    void rayleigh_ritz(T* psi, Real* eigenvalue,
                       std::vector<int>& active_cols,
                       const std::vector<double>& ethr_band);

    // -------------------------------------------------------------------------
    // Conjugate-gradient strategy helpers (File 2 style)
    // -------------------------------------------------------------------------
    void calc_gradient(const Real* prec,
                       const T* hpsi,
                       const T* spsi,
                       const T* psi,
                       const Real* eigenvalue,
                       std::vector<T>& grad) const;

    void orth_gradient(const T* psi, const T* spsi,
                       std::vector<T>& grad) const;

    void update_polak_ribiere(const std::vector<T>& grad,
                              std::vector<T>& p,
                              std::vector<T>& z_old,
                              std::vector<Real>& beta_denom,
                              const Real* prec) const;

    void line_minimize(T* psi, T* hpsi, T* spsi,
                       const T* p, const T* hp, const T* sp,
                       int ncol) const;

    void orth_cholesky(T* psi, T* hpsi, T* spsi, int ncol) const;
};

} // namespace hsolver

#endif // DIAGO_PPCG_H

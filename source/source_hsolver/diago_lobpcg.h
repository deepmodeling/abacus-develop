#ifndef DIAGO_LOBPCG_H_
#define DIAGO_LOBPCG_H_

#include <complex>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/para_gemm.h"
#include "source_hamilt/hamilt.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_hsolver/para_linear_transform.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <source_base/macros.h>

namespace hsolver {

/**
 * @class DiagoLobpcg
 * @brief Locally Optimal Block Preconditioned Conjugate Gradient eigensolver.
 *
 * LOBPCG maintains a block Rayleigh-Ritz subspace:
 *   - First iteration:  W = [X, Z]           (2-block, no valid P yet)
 *   - Subsequent:       W = [X, Z, P]        (3-block)
 * where X = current eigenvectors, Z = preconditioned residual, P = search directions.
 *
 * @note Currently supports CPU only.
 *       GPU support is planned for subsequent phases.
 *
 * @tparam T      Complex floating-point type.
 * @tparam Device Must be base_device::DEVICE_CPU (GPU not yet supported).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoLobpcg
{
  private:
    static_assert(std::is_same<T, std::complex<double>>::value,
                  "DiagoLobpcg currently supports std::complex<double> only.");
    static_assert(std::is_same<Device, base_device::DEVICE_CPU>::value,
                  "DiagoLobpcg currently supports CPU only.");

    using Real = typename GetTypeReal<T>::type;

  public:
    /// @brief H * psi -> hpsi.
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /// @brief S * psi -> spsi.
    using SPsiFunc = std::function<void(const T*, T*, const int, const int)>;

    /// Constructor: stores host preconditioner pointer.
    explicit DiagoLobpcg(const Real* precondition);

    ~DiagoLobpcg();

    /// Allocate workspace and bind h_prec TensorMap (after n_basis is known).
    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    /// Set max inner iterations per SCF step (default 4).
    void set_nline(const int n) { this->nline = n; }

    /// Set hard iteration limit from pw_diag_nmax. Non-positive keeps legacy default.
    void set_max_iter(const int n) { this->max_iter = n; }

    /// Set allowed unconverged bands after max_iter. Negative means report only.
    void set_notconv_max(const int n) { this->notconv_max = n; }

    /// Keep strict failure for NSCF/tests; SCF should report and let outer mixing continue.
    void set_throw_on_notconv_exceed(const bool enabled) { this->throw_on_notconv_exceed = enabled; }

    void set_diag_context(const std::string& context) { this->diag_context = context; }

    /// Generalized diagonalization for non-identity overlap S.
    int diag(const HPsiFunc& hpsi_func,
             const SPsiFunc& spsi_func,
             T* psi_in,
             Real* eigenvalue_in,
             const std::vector<double>& ethr_band);

  private:
    // ---- dimensions ----
    int n_band   = 0;   ///< total bands (global)
    int n_band_l = 0;   ///< local bands
    int n_basis  = 0;   ///< basis functions (lda of psi)
    int n_dim    = 0;   ///< valid dimension (= current_ngk)
    int nline    = 4;   ///< max inner iterations per SCF step
    int max_iter = 0;   ///< hard iteration limit; <=0 uses nline-based default
    int notconv_max = -1; ///< allowed unconverged bands; <0 reports only
    bool throw_on_notconv_exceed = true; ///< throw when notconv exceeds notconv_max after max_iter
    int nsub     = 0;   ///< physical leading dim of hsub (= 3*n_band)
    bool has_pdir = false; ///< true when P block holds valid directions
    std::string diag_context;

    struct State
    {
        std::vector<T> psi;
        std::vector<T> hpsi;
        std::vector<T> spsi;
        std::vector<Real> eigen;
    };
    struct StateQuality
    {
        Real residual = std::numeric_limits<Real>::max();
        int notconv = std::numeric_limits<int>::max();
        bool valid = false;
    };

    // ---- parallel ops ----
    ModuleBase::PGemmCN<T, Device> pmmcn;
    PLinearTransform<T, Device> plintrans;
    PLinearTransform<T, Device> plintrans_batch2;
    PLinearTransform<T, Device> plintrans_batch3;

    // ---- type traits ----
    ct::DataType r_type   = ct::DataType::DT_INVALID;
    ct::DataType t_type   = ct::DataType::DT_INVALID;
    ct::DeviceType dev_type = ct::DeviceType::UnKnown;

    // ---- preconditioner ----
    ct::Tensor prec        = {};  ///< device copy [n_basis]
    const Real* h_prec_ptr = nullptr; ///< host pointer (saved in ctor)
    ct::Tensor h_prec      = {};  ///< host TensorMap [n_basis] (bound in init_iter)

    // ---- eigenvalues & convergence ----
    ct::Tensor eigen     = {}; ///< output eigenvalues [n_band]
    ct::Tensor sub_eigen = {}; ///< subspace eigenvalues [nsub] = [3*n_band]
    ct::Tensor err_st    = {}; ///< residual norm per band [n_band_l]

    // ---- core blocks ----
    // Layout for {n_band_l, n_basis} tensors:
    //   data[ib * n_basis + ig]: band-major contiguous.
    //   BLAS view: n_basis rows x n_band_l cols, column-major.
    ct::Tensor psi   = {};   ///< X (TensorMap to psi_in, no ownership)
    ct::Tensor hpsi  = {};   ///< HX
    ct::Tensor spsi  = {};   ///< SX
    ct::Tensor grad  = {};   ///< Z = T(R)
    ct::Tensor hgrad = {};   ///< HZ
    ct::Tensor sgrad = {};   ///< SZ
    ct::Tensor pdir  = {};   ///< P
    ct::Tensor hpdir = {};   ///< HP
    ct::Tensor spdir = {};   ///< SP

    // ---- subspace matrices ----
    ct::Tensor hsub = {};    ///< H_sub [nsub x nsub]
    ct::Tensor ssub = {};    ///< S_sub [nsub x nsub]

    // ---- workspace ----
    ct::Tensor work    = {}; ///< [n_band_l, n_basis]
    ct::Tensor hwork   = {}; ///< [n_band_l, n_basis]
    ct::Tensor swork   = {}; ///< [n_band_l, n_basis]
    ct::Tensor pwork   = {}; ///< [n_band_l, n_basis] (P update)
    ct::Tensor hpwork  = {}; ///< [n_band_l, n_basis] (HP update)
    ct::Tensor spwork  = {}; ///< [n_band_l, n_basis] (SP update)
    ct::Tensor tmp_hsub = {}; ///< scratch [n_band, n_band]
    ct::Tensor tmp_ssub = {}; ///< scratch [n_band, n_band]
    std::vector<T> plintrans_batch_in;
    std::vector<T> plintrans_batch_out;

    // ---- GEMM constants (following BPCG pattern) ----
    const T one_     = static_cast<T>(1.0);
    const T zero_    = static_cast<T>(0.0);
    const T neg_one_ = static_cast<T>(-1.0);
    const T* one     = nullptr;
    const T* zero    = nullptr;

    // ---- helpers ----

    void calc_prec();

    void calc_hpsi_with_block(const HPsiFunc& hpsi_func,
                              T* psi_in,
                              ct::Tensor& hpsi_out);

    void calc_spsi_with_block(const SPsiFunc& spsi_func,
                              const T* psi_in,
                              ct::Tensor& spsi_out);

    void repair_initial_subspace_s(const HPsiFunc& hpsi_func,
                                   const SPsiFunc& spsi_func);

    /// Generalized R-R.
    void generalized_rayleigh_ritz(ct::Tensor& psi_inout,
                                   ct::Tensor& hpsi_inout,
                                   ct::Tensor& spsi_inout,
                                   ct::Tensor& eigen_out);

    /// Distributed generalized R-R for band-parallel S != I.
    void generalized_rayleigh_ritz_parallel(ct::Tensor& psi_inout,
                                            ct::Tensor& hpsi_inout,
                                            ct::Tensor& spsi_inout,
                                            ct::Tensor& eigen_out);

    /// Generalized residual.
    void compute_residual_s(const ct::Tensor& psi_in,
                            const ct::Tensor& hpsi_in,
                            const ct::Tensor& spsi_in,
                            const ct::Tensor& eigen_in,
                            const ct::Tensor& prec_in,
                            ct::Tensor& grad_out,
                            ct::Tensor& err_out);

    /// S-orthogonalize.
    void orth_projection_s(const ct::Tensor& psi_in,
                           const ct::Tensor& spsi_in,
                           ct::Tensor& hsub_work,
                           ct::Tensor& sgrad_out,
                           ct::Tensor& grad_out);

    void orth_projection_s_with_h(const ct::Tensor& psi_in,
                                  const ct::Tensor& hpsi_in,
                                  const ct::Tensor& spsi_in,
                                  ct::Tensor& hsub_work,
                                  ct::Tensor& hpdir_out,
                                  ct::Tensor& spdir_out,
                                  ct::Tensor& pdir_out);

    /// Generalized subspace update for S != I.
    void lobpcg_update_s(ct::Tensor& psi,
                         ct::Tensor& hpsi,
                         ct::Tensor& spsi,
                         ct::Tensor& grad,
                         ct::Tensor& hgrad,
                         ct::Tensor& sgrad,
                         ct::Tensor& pdir,
                         ct::Tensor& hpdir,
                         ct::Tensor& spdir,
                         ct::Tensor& eigen);

    /// Distributed generalized subspace update for band-parallel S != I.
    void lobpcg_update_s_parallel(ct::Tensor& psi,
                                  ct::Tensor& hpsi,
                                  ct::Tensor& spsi,
                                  ct::Tensor& grad,
                                  ct::Tensor& hgrad,
                                  ct::Tensor& sgrad,
                                  ct::Tensor& pdir,
                                  ct::Tensor& hpdir,
                                  ct::Tensor& spdir,
                                  ct::Tensor& eigen,
                                  bool force_compressed);

    /// psi = psi * U  (via plintrans).
    void rotate_wf(const ct::Tensor& hsub_in,
                   ct::Tensor& psi_out,
                   ct::Tensor& workspace_in);

    /// S-Cholesky orthonormalization.
    void orth_cholesky_s(ct::Tensor& workspace_in,
                         ct::Tensor& psi_out,
                         ct::Tensor& hpsi_out,
                         ct::Tensor& spsi_out,
                         ct::Tensor& hsub_out);

    Real max_error(const ct::Tensor& err_in) const;

    int count_not_converged(const ct::Tensor& err_in,
                            const std::vector<double>& ethr_band) const;

    void validate_ethr_band(const std::vector<double>& ethr_band) const;

    void diag_log(const std::string& context,
                  const std::string& line1,
                  const std::string& line2,
                  const std::string& line3 = std::string()) const;
    std::string error_with_context(const std::string& message) const;

    void report_not_converged(const char* problem_type,
                              const int used_iter,
                              const int max_iter,
                              const std::vector<double>& ethr_band,
                              const std::string& stop_reason) const;

    void pgemm_multiply(const T alpha,
                        const T* a,
                        const T* b,
                        const T beta,
                        T* c);
    void plintrans_act(const T alpha,
                       const T* a,
                       const T* u,
                       const T beta,
                       T* b);
    void plintrans_batched_act(const T alpha,
                               const T* const* a_blocks,
                               const int block_count,
                               const T* u,
                               const T beta,
                               T* const* b_blocks);

    int local_band_start() const;
    void clear_search_directions();
    void save_state(std::vector<T>& psi_out,
                    std::vector<T>& hpsi_out,
                    std::vector<Real>& eigen_out);
    void restore_state(const std::vector<T>& psi_in,
                       const std::vector<T>& hpsi_in,
                       const std::vector<Real>& eigen_in,
                       const bool reset_search);
    void save_state(State& state);
    void restore_state(const State& state,
                       const bool reset_search);
    std::vector<char> make_soft_lock_mask(const ct::Tensor& err_in,
                                          const std::vector<double>& ethr_band,
                                          const int notconv) const;
    bool restore_soft_locked_bands(const State& state,
                                   const std::vector<char>& soft_lock_mask,
                                   const ct::Tensor& err_in,
                                   const std::vector<double>& ethr_band);
    bool update_best_state(State& state,
                           StateQuality& quality,
                           const int candidate_notconv,
                           const Real candidate_residual);
    template <typename UpdateFunc, typename RepairFunc>
    void update_subspace_with_fallback(const char* first_failure,
                                       const char* retry_failure,
                                       State& rollback_state,
                                       const int used_iter,
                                       const UpdateFunc& update_func,
                                       const RepairFunc& repair_func);
    template <typename RecomputeFunc>
    void restore_best_state_if_needed(const char* diag_name,
                                      State& best_state,
                                      const StateQuality& best_quality,
                                      const int used_iter,
                                      Real& final_residual,
                                      int& final_notconv,
                                      const RecomputeFunc& recompute_final_quality);
    int run_lobpcg_loop(const char* problem_type,
                        const char* diag_name,
                        const HPsiFunc& hpsi_func,
                        const SPsiFunc& spsi_func,
                        const std::vector<double>& effective_ethr_band,
                        const int scf_iter);
    bool handle_generalized_rejected_update(const State& rollback_state,
                                            State& best_state,
                                            StateQuality& best_quality,
                                            bool& update_rejected,
                                            const int used_iter,
                                            const int notconv_before_update,
                                            const Real residual_before_update,
                                            const Real residual_after_update,
                                            const Real residual_limit,
                                            const Real residual_growth_limit,
                                            const std::vector<double>& effective_ethr_band);
    void ensure_generalized_update_finite(const ct::Tensor& psi,
                                          const ct::Tensor& hpsi,
                                          const ct::Tensor& spsi,
                                          const ct::Tensor& eigen,
                                          const int nbs,
                                          const int n,
                                          const char* log_context,
                                          const char* error_message) const;
    void update_generalized_subspace(const bool force_compressed);
    void repair_generalized_subspace(const HPsiFunc& hpsi_func,
                                     const SPsiFunc& spsi_func);
    void project_generalized_search_direction(const HPsiFunc& hpsi_func,
                                              const SPsiFunc& spsi_func);

    // ---- memory-op aliases ----
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_var_op    = ct::kernels::set_memory<Real, ct_Device>;
    using syncmem_var_h2d_op = ct::kernels::synchronize_memory<Real, ct_Device, ct::DEVICE_CPU>;
    using syncmem_var_d2h_op = ct::kernels::synchronize_memory<Real, ct::DEVICE_CPU, ct_Device>;

    using setmem_complex_op    = ct::kernels::set_memory<T, ct_Device>;
    using syncmem_complex_op   = ct::kernels::synchronize_memory<T, ct_Device, ct_Device>;
    using syncmem_complex_h2d_op = ct::kernels::synchronize_memory<T, ct_Device, ct::DEVICE_CPU>;
    using syncmem_complex_d2h_op = ct::kernels::synchronize_memory<T, ct::DEVICE_CPU, ct_Device>;

};

} // namespace hsolver
#endif // DIAGO_LOBPCG_H_

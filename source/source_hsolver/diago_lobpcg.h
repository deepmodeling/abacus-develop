#ifndef DIAGO_LOBPCG_H_
#define DIAGO_LOBPCG_H_

#include <complex>
#include <functional>
#include <type_traits>
#include <vector>

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/para_gemm.h"
#include "source_hamilt/hamilt.h"
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
 * @note Currently supports NC (S=I) + CPU only.
 *       GPU and USPP support are planned for subsequent phases.
 *
 * @tparam T      Complex floating-point type.
 * @tparam Device Must be base_device::DEVICE_CPU (GPU not yet supported).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoLobpcg
{
  private:
    static_assert(std::is_same<Device, base_device::DEVICE_CPU>::value,
                  "DiagoLobpcg currently supports CPU only.");

    using Real = typename GetTypeReal<T>::type;

  public:
    /// @brief H * psi -> hpsi.
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /// @brief S * psi -> spsi.
    using SPsiFunc = std::function<void(const T*, T*, const int, const int)>;

    /// Constructor — stores host preconditioner pointer.
    explicit DiagoLobpcg(const Real* precondition);

    ~DiagoLobpcg();

    /// Allocate workspace and bind h_prec TensorMap (after n_basis is known).
    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    /// Set max inner iterations per SCF step (default 4).
    void set_nline(const int n) { this->nline = n; }

    /// Generalized diagonalization interface. Currently supports S = I only.
    void diag(const HPsiFunc& hpsi_func,
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
    int nsub     = 0;   ///< physical leading dim of hsub (= 3*n_band)
    bool has_pdir = false; ///< true when P block holds valid directions

    // ---- parallel ops ----
    ModuleBase::PGemmCN<T, Device> pmmcn;
    PLinearTransform<T, Device> plintrans;

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
    //   data[ib * n_basis + ig] — band-major contiguous.
    //   BLAS view: n_basis rows × n_band_l cols, column-major.
    ct::Tensor psi   = {};   ///< X (TensorMap → psi_in, no ownership)
    ct::Tensor hpsi  = {};   ///< HX
    ct::Tensor spsi  = {};   ///< SX  (Phase 2)
    ct::Tensor grad  = {};   ///< Z = T(R)
    ct::Tensor hgrad = {};   ///< HZ
    ct::Tensor pdir  = {};   ///< P
    ct::Tensor hpdir = {};   ///< HP

    // ---- subspace matrices ----
    ct::Tensor hsub = {};    ///< H_sub [nsub × nsub]

    // ---- workspace ----
    ct::Tensor work    = {}; ///< [n_band_l, n_basis]
    ct::Tensor hwork   = {}; ///< [n_band_l, n_basis]
    ct::Tensor pwork   = {}; ///< [n_band_l, n_basis] (P update)
    ct::Tensor hpwork  = {}; ///< [n_band_l, n_basis] (HP update)
    ct::Tensor tmp_hsub = {}; ///< scratch [n_band, n_band]

    // ---- GEMM constants (following BPCG pattern) ----
    Device* ctx = {};
    const T one_     = static_cast<T>(1.0);
    const T zero_    = static_cast<T>(0.0);
    const T neg_one_ = static_cast<T>(-1.0);
    const T* one     = nullptr;
    const T* zero    = nullptr;
    const T* neg_one = nullptr;

    // ---- helpers ----

    void calc_prec();

    void calc_hpsi_with_block(const HPsiFunc& hpsi_func,
                              T* psi_in,
                              ct::Tensor& hpsi_out);

    void calc_spsi_with_block(const SPsiFunc& spsi_func,
                              const T* psi_in,
                              ct::Tensor& spsi_out);

    /// Standard R-R: H_sub = psi^H * hpsi → heevd → rotate.
    /// multiply(hpsi, psi) = psi^H * hpsi.
    void rayleigh_ritz(ct::Tensor& psi_inout,
                       ct::Tensor& hpsi_inout,
                       ct::Tensor& eigen_out);

    /// Generalized R-R. NOT IMPLEMENTED — aborts.
    void generalized_rayleigh_ritz(ct::Tensor& psi_inout,
                                   ct::Tensor& hpsi_inout,
                                   ct::Tensor& spsi_inout,
                                   ct::Tensor& eigen_out);

    /// NC residual: R = HX - X*Lambda, Z = R ./ prec.
    /// CPU-only: direct loops.
    void compute_residual(const ct::Tensor& psi_in,
                          const ct::Tensor& hpsi_in,
                          const ct::Tensor& eigen_in,
                          const ct::Tensor& prec_in,
                          ct::Tensor& grad_out,
                          ct::Tensor& err_out);

    /// USPP residual. NOT IMPLEMENTED — aborts.
    void compute_residual_s(const ct::Tensor& psi_in,
                            const ct::Tensor& hpsi_in,
                            const ct::Tensor& spsi_in,
                            const ct::Tensor& eigen_in,
                            const ct::Tensor& prec_in,
                            ct::Tensor& grad_out,
                            ct::Tensor& err_out);

    /// grad -= psi * (psi^H * grad).  multiply(grad, psi) = psi^H * grad.
    void orth_projection(const ct::Tensor& psi_in,
                         ct::Tensor& hsub_work,
                         ct::Tensor& grad_out);

    /// S-orthogonalize. NOT IMPLEMENTED — aborts.
    void orth_projection_s(const ct::Tensor& psi_in,
                           const ct::Tensor& spsi_in,
                           ct::Tensor& hsub_work,
                           ct::Tensor& grad_out);

    /// Core subspace update (2-block first, then 3-block).
    /// Orthonormalizes W = [X, Z, P] before Rayleigh-Ritz for stability.
    void lobpcg_update(ct::Tensor& psi,
                       ct::Tensor& hpsi,
                       ct::Tensor& grad,
                       ct::Tensor& hgrad,
                       ct::Tensor& pdir,
                       ct::Tensor& hpdir,
                       ct::Tensor& eigen);

    /// psi = psi * U  (via plintrans).
    void rotate_wf(const ct::Tensor& hsub_in,
                   ct::Tensor& psi_out,
                   ct::Tensor& workspace_in);

    /// S=I Cholesky: psi^H*psi → potrf(U) → trtri → psi *= U^{-1}, hpsi *= U^{-1}.
    void orth_cholesky(ct::Tensor& workspace_in,
                       ct::Tensor& psi_out,
                       ct::Tensor& hpsi_out,
                       ct::Tensor& hsub_out);

    /// S-Cholesky. NOT IMPLEMENTED — aborts.
    void orth_cholesky_s(ct::Tensor& workspace_in,
                         ct::Tensor& psi_out,
                         ct::Tensor& hpsi_out,
                         ct::Tensor& spsi_out,
                         ct::Tensor& hsub_out);

    bool test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band);

    // ---- memory-op aliases ----
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_var_op    = ct::kernels::set_memory<Real, ct_Device>;
    using resmem_var_op    = ct::kernels::resize_memory<Real, ct_Device>;
    using delmem_var_op    = ct::kernels::delete_memory<Real, ct_Device>;
    using syncmem_var_h2d_op = ct::kernels::synchronize_memory<Real, ct_Device, ct::DEVICE_CPU>;
    using syncmem_var_d2h_op = ct::kernels::synchronize_memory<Real, ct::DEVICE_CPU, ct_Device>;

    using setmem_complex_op    = ct::kernels::set_memory<T, ct_Device>;
    using delmem_complex_op    = ct::kernels::delete_memory<T, ct_Device>;
    using resmem_complex_op    = ct::kernels::resize_memory<T, ct_Device>;
    using syncmem_complex_op   = ct::kernels::synchronize_memory<T, ct_Device, ct_Device>;
    using syncmem_complex_h2d_op = ct::kernels::synchronize_memory<T, ct_Device, ct::DEVICE_CPU>;
    using syncmem_complex_d2h_op = ct::kernels::synchronize_memory<T, ct::DEVICE_CPU, ct_Device>;
};

} // namespace hsolver
#endif // DIAGO_LOBPCG_H_

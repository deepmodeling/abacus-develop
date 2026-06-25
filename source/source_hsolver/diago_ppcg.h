#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/para_gemm.h"
#include "source_hsolver/para_linear_transform.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/lapack.h>
#include <source_base/macros.h>

#include <functional>
#include <vector>

namespace hsolver {

/**
 * @class DiagoPPCG
 * @brief A class for diagonalization using the Projected-PCG method.
 *
 * Key differences from BPCG:
 * 1. Polak-Ribiere beta formula for conjugate direction mixing (BPCG uses
 *    Fletcher-Reeves), giving better convergence for non-quadratic problems.
 * 2. Soft-locking of converged bands: once a band's residual falls below its
 *    convergence threshold, its search direction is zeroed in subsequent
 *    iterations, reducing the effective problem size.
 * 3. Conjugate history (z_old, grad_old) is consistently rotated alongside
 *    the wavefunctions during Cholesky orthonormalization, preserving the
 *    conjugacy relationship across basis rotations.
 *
 * @tparam T The floating-point type used for calculations.
 * @tparam Device The device used for calculations (cpu only).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoPPCG
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    explicit DiagoPPCG(const Real* precondition);
    ~DiagoPPCG();

    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    void diag(const HPsiFunc& hpsi_func,
              T* psi_in,
              Real* eigenvalue_in,
              const std::vector<double>& ethr_band);

  private:
    int n_band = 0;
    int n_band_l = 0;
    int n_basis = 0;
    int n_dim = 0;
    int nline = 4;

    ModuleBase::PGemmCN<T, Device> pmmcn;
    PLinearTransform<T, Device> plintrans;

    ct::DataType r_type  = ct::DataType::DT_INVALID;
    ct::DataType t_type  = ct::DataType::DT_INVALID;
    ct::DeviceType device_type = ct::DeviceType::UnKnown;

    ct::Tensor prec = {}, h_prec = {};

    /// Polak-Ribiere beta for conjugate direction mixing
    ct::Tensor beta = {};
    /// Error state per band
    ct::Tensor err_st = {};
    /// Eigenvalues
    ct::Tensor eigen = {};

    /// Wavefunction and H|psi>
    ct::Tensor psi = {}, hpsi = {};

    /// Subspace Hamiltonian
    ct::Tensor hsub = {};

    /// Search direction d (mixed preconditioned gradient)
    ct::Tensor grad = {}, hgrad = {}, grad_old = {};

    /// Previous preconditioned gradient z_old (before mixing),
    /// needed for the Polak-Ribiere beta formula.
    ct::Tensor z_old = {};

    /// Work array
    ct::Tensor work = {};

    /// Per-band convergence mask: 1 = converged (locked), 0 = active.
    std::vector<char> conv_mask;

    /// Update the precondition array from host to device.
    void calc_prec();

    /// Apply H to psi: hpsi_out = H |psi_in>
    void calc_hpsi_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& hpsi_out);

    /// Subspace diagonalization: solve H_sub * c = eps * c
    void diag_hsub(
        const ct::Tensor& psi_in,
        const ct::Tensor& hpsi_in,
        ct::Tensor& hsub_out,
        ct::Tensor& eigenvalue_out);

    /// Rotate wavefunctions: psi_out = psi_out * hsub_in
    void rotate_wf(
        const ct::Tensor& hsub_in,
        ct::Tensor& psi_out,
        ct::Tensor& workspace_in);

    /// Compute the PPCG search direction with Polak-Ribiere mixing.
    /// Converged bands (conv_mask[ib] == 1) are skipped and their
    /// search direction is zeroed out (soft-locking).
    void calc_grad_with_block(
        const ct::Tensor& prec_in,
        ct::Tensor& err_out,
        ct::Tensor& beta_out,
        ct::Tensor& psi_in,
        ct::Tensor& hpsi_in,
        ct::Tensor& grad_out,
        ct::Tensor& grad_old_out);

    /// Project grad to be orthogonal to all columns of psi:
    /// grad = grad - psi * (psi^T * grad)
    void orth_projection(
        const ct::Tensor& psi_in,
        ct::Tensor& hsub_in,
        ct::Tensor& grad_out);

    /// Line minimization: optimize psi along the conjugate direction.
    void line_minimize(
        ct::Tensor& grad_in,
        ct::Tensor& hgrad_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out);

    /// Cholesky orthonormalization of psi, hpsi, and conjugate history.
    /// Also rotates z_old and grad_old by the Cholesky factor so that
    /// conjugate history remains consistent with the current basis.
    void orth_cholesky(
        ct::Tensor& workspace_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out);

    /// Initial subspace setup: H|psi>, diag, rotate.
    void calc_hsub_with_block(
        const HPsiFunc& hpsi_func,
        T *psi_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out,
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out);

    /// Final subspace cleanup: diag and rotate psi only.
    void calc_hsub_with_block_exit(
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out,
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out);

    /// Check convergence: returns true if any band still unconverged.
    bool test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band);

    /// Update conv_mask from current err_st vs per-band thresholds.
    void update_convergence(const std::vector<double>& ethr_band);

    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_var_op = ct::kernels::set_memory<Real, ct_Device>;
    using syncmem_var_h2d_op = ct::kernels::synchronize_memory<Real, ct_Device, ct::DEVICE_CPU>;
    using syncmem_var_d2h_op = ct::kernels::synchronize_memory<Real, ct::DEVICE_CPU, ct_Device>;

    using setmem_complex_op = ct::kernels::set_memory<T, ct_Device>;
    using syncmem_complex_op = ct::kernels::synchronize_memory<T, ct_Device, ct_Device>;
};

} // namespace hsolver
#endif // DIAGO_PPCG_H_

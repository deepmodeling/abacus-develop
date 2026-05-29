#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/para_gemm.h"
#include "source_hsolver/kernels/hegvd_op.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <source_base/macros.h>

namespace hsolver {

// Projected Preconditioned Conjugate Gradient (block) eigensolver.
// This implementation follows an LOBPCG-style subspace projection:
//   V = [X, W, P] (or [X, W] for the first iter)
//   solve (V^H H V) c = (V^H V) c Λ
//   update X <- V c(:,1:nband)
//   update P <- [W, P] c_{W,P}(:,1:nband)
// with W from preconditioned residual projected to the complement of X (and P).
//
// Notes:
// - Designed to match the existing diag interface used by BPCG.
// - Preconditioner is treated as a diagonal Real vector of length n_basis.

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

    ct::DataType r_type = ct::DataType::DT_INVALID;
    ct::DataType t_type = ct::DataType::DT_INVALID;
    ct::DeviceType device_type = ct::DeviceType::UnKnown;

    // Host pointer mapped preconditioner + device copy
    ct::Tensor h_prec = {};
    ct::Tensor prec = {};

    // Work vectors (column-major, lda = n_basis): each tensor stores a matrix (n_basis x ncols)
    // as a contiguous array; we use Tensor with shape {ncols, n_basis} for contiguous column blocks.
    ct::Tensor X = {};   // mapped to psi_in
    ct::Tensor HX = {};
    ct::Tensor R = {};
    ct::Tensor W = {};
    ct::Tensor HW = {};
    ct::Tensor P = {};
    ct::Tensor HP = {};

    ct::Tensor V = {};   // basis [X,W,P] packed
    ct::Tensor HV = {};  // H*V

    ct::Tensor hcc = {}; // V^H H V
    ct::Tensor scc = {}; // V^H V
    ct::Tensor vcc = {}; // eigenvectors of projected problem
    ct::Tensor eval = {}; // eigenvalues of projected problem

    ct::Tensor work = {}; // generic workspace (ncols x n_basis)

    // Parallel matmul helper (A^H B)
    ModuleBase::PGemmCN<T, Device> pmmcn;

    // Device memory helpers
    Device* ctx = {};
    const T one_ = static_cast<T>(1.0);
    const T zero_ = static_cast<T>(0.0);
    const T neg_one_ = static_cast<T>(-1.0);
    const T* one = &one_;
    const T* zero = &zero_;
    const T* neg_one = &neg_one_;

    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_var_op = base_device::memory::synchronize_memory_op<Real, Device, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;

    void calc_prec();

    void apply_h(const HPsiFunc& hpsi_func, const ct::Tensor& in_vecs, ct::Tensor& out_vecs, const int nvec);

    void pack_basis(const int ncols, const bool has_p);

    void compute_projected_mats(const int ncols);

    void solve_projected(const int ncols);

    void update_from_projected(const int ncols, const bool has_p);

    void compute_residual_and_precond(const std::vector<double>& ethr_band, bool& not_conv);

    void orthonormalize_block(ct::Tensor& A, ct::Tensor* HA, const int ncols);

    void project_out(const ct::Tensor& basis, const int ncols_basis, ct::Tensor& vecs, const int ncols_vecs);

    bool check_convergence(const ct::Tensor& residual, const std::vector<double>& ethr_band);
};

} // namespace hsolver

#endif // DIAGO_PPCG_H_

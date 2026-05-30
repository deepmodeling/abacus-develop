#ifndef MODULE_HSOLVER_DIAGO_CG_H_
#define MODULE_HSOLVER_DIAGO_CG_H_

#include <functional>
#include <vector>

#include <source_base/macros.h>
#include <source_base/kernels/math_kernel_op.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

namespace hsolver {

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagoCG final
{
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;

  public:
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;
    using SPsiFunc = std::function<void(T*, T*, const int, const int)>;
    using SubspaceFunc = std::function<void(T*, T*, const int, const int, const bool)>;

    DiagoCG(const std::string& basis_type, const std::string& calculation);
    DiagoCG(const std::string& basis_type,
            const std::string& calculation,
            const bool& need_subspace,
            const SubspaceFunc& subspace_func,
            const Real& pw_diag_thr,
            const int& pw_diag_nmax,
            const int& nproc_in_pool);

    ~DiagoCG();

    double diag(const HPsiFunc& hpsi_func,
                const SPsiFunc& spsi_func,
                const int ld_psi,
                const int nband,
                const int dim,
                T* psi_in,
                Real* eigenvalue_in,
                const std::vector<double>& ethr_band,
                const Real* prec = nullptr);

  private:
    Device * ctx_ = {};
    int notconv_ = 0;
    int n_band_ = 0;
    int n_basis_ = 0;
    double avg_iter_ = 0;
    std::vector<int> iter_band;
    Real pw_diag_thr_ = 1e-5;
    int pw_diag_nmax_ = 0;
    int nproc_in_pool_ = 0;
    std::string basis_type_ = {};
    std::string calculation_ = {};

    bool need_subspace_ = false;
    HPsiFunc hpsi_func_ = nullptr;
    SPsiFunc spsi_func_ = nullptr;
    SubspaceFunc subspace_func_ = nullptr;

    // P1: 类成员 workspace，避免 diag_once 中重复分配
    ct::Tensor ws_phi_m_;
    ct::Tensor ws_hphi_;
    ct::Tensor ws_sphi_;
    ct::Tensor ws_pphi_;
    ct::Tensor ws_cg_;
    ct::Tensor ws_scg_;
    ct::Tensor ws_grad_;
    ct::Tensor ws_g0_;
    ct::Tensor ws_lagrange_;
    ct::Tensor ws_temp_;  // orth_grad 中的投影临时向量

    void ensure_workspace(int n_basis, int n_band);

    void calc_grad(const ct::Tensor& prec,
                   ct::Tensor& grad,
                   ct::Tensor& hphi,
                   ct::Tensor& sphi,
                   ct::Tensor& pphi);

    void orth_grad(const ct::Tensor& psi,
                   const int& m,
                   ct::Tensor& grad,
                   ct::Tensor& scg,
                   ct::Tensor& lagrange);

    void calc_gamma_cg(const int& iter,
                       const Real& cg_norm,
                       const Real& theta,
                       const ct::Tensor& prec,
                       const ct::Tensor& scg,
                       const ct::Tensor& grad,
                       const ct::Tensor& phi_m,
                       Real& gg_last,
                       ct::Tensor& g0,
                       ct::Tensor& cg);

    bool update_psi(const ct::Tensor& pphi,
                    const ct::Tensor& cg,
                    const ct::Tensor& scg,
                    const double& ethreshold,
                    Real &cg_norm,
                    Real &theta,
                    Real &eigen,
                    ct::Tensor& phi_m,
                    ct::Tensor& sphi,
                    ct::Tensor& hphi);

    void schmit_orth(const int& m, const ct::Tensor& psi, const ct::Tensor& sphi, ct::Tensor& phi_m);

    void diag_once(const ct::Tensor& prec,
                   ct::Tensor& psi,
                   ct::Tensor& eigen,
                   const std::vector<double>& ethr_band);

    bool test_exit_cond(const int& ntry, const int& notconv) const;

    using dot_real_op = ModuleBase::dot_real_op<T, Device>;
    const T * one_ = nullptr, * zero_ = nullptr, * neg_one_ = nullptr;
};

} // namespace hsolver

#endif
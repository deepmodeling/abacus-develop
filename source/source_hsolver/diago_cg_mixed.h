#ifndef SOURCE_HSOLVER_DIAGO_CG_MIXED_H_
#define SOURCE_HSOLVER_DIAGO_CG_MIXED_H_

#include <functional>
#include <source_base/macros.h>
#include <source_base/kernels/math_kernel_op.h>
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

namespace hsolver {

template <typename T> struct GetFloatType { using type = T; };
template <> struct GetFloatType<std::complex<double>> { using type = std::complex<float>; };
template <> struct GetFloatType<double> { using type = float; };

template <typename T> struct GetFloatRealType { using type = typename GetTypeReal<T>::type; };
template <> struct GetFloatRealType<std::complex<double>> { using type = float; };
template <> struct GetFloatRealType<double> { using type = float; };

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagoCGMixed final
{
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using T_float = typename GetFloatType<T>::type;
    using Real_float = typename GetFloatRealType<T>::type;

  public:
    using Func = std::function<void(const ct::Tensor&, ct::Tensor&)>;

    DiagoCGMixed(const std::string& basis_type, const std::string& calculation);
    DiagoCGMixed(const std::string& basis_type, const std::string& calculation,
                 const bool& need_subspace, const Func& subspace_func,
                 const Real& pw_diag_thr, const int& pw_diag_nmax, const int& nproc_in_pool);
    ~DiagoCGMixed() = default;

    void diag(const Func& hpsi_func, const Func& spsi_func,
              ct::Tensor& psi, ct::Tensor& eigen,
              const std::vector<double>& ethr_band, const ct::Tensor& prec = {});

  private:
    int notconv_ = 0, n_band_ = 0, n_basis_ = 0, avg_iter_ = 0, pw_diag_nmax_ = 0, nproc_in_pool_ = 0;
    Real pw_diag_thr_ = 1e-5;
    std::string basis_type_, calculation_;
    bool need_subspace_ = false;
    Func hpsi_func_, spsi_func_, subspace_func_;

    void convert_d2f(const ct::Tensor& d, ct::Tensor& f);
    void convert_f2d(const ct::Tensor& f, ct::Tensor& d);

    void calc_grad(const ct::Tensor& prec_f, ct::Tensor& grad, ct::Tensor& hphi,
                   ct::Tensor& sphi, ct::Tensor& pphi);
    void orth_grad(const ct::Tensor& psi, const int& m, ct::Tensor& grad,
                   ct::Tensor& scg, ct::Tensor& lagrange);
    void calc_gamma_cg(const int& iter, const Real& cg_norm, const Real& theta,
                       const ct::Tensor& prec_f, const ct::Tensor& scg,
                       const ct::Tensor& grad, const ct::Tensor& phi_m,
                       Real& gg_last, ct::Tensor& g0, ct::Tensor& cg);
    bool update_psi(const ct::Tensor& pphi, const ct::Tensor& cg, const ct::Tensor& scg,
                    const double& ethreshold, Real& cg_norm, Real& theta, Real& eigen,
                    ct::Tensor& phi_m, ct::Tensor& sphi, ct::Tensor& hphi);
    void schmit_orth(const int& m, const ct::Tensor& psi, const ct::Tensor& sphi, ct::Tensor& phi_m);
    void diag_once(const ct::Tensor& prec, ct::Tensor& psi, ct::Tensor& eigen,
                   const std::vector<double>& ethr_band);
    bool test_exit_cond(const int& ntry, const int& notconv) const;
};

} // namespace hsolver
#endif

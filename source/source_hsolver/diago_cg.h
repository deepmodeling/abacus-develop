#ifndef MODULE_HSOLVER_DIAGO_CG_H_
#define MODULE_HSOLVER_DIAGO_CG_H_

#include <functional>

#include <module_base/macros.h>
#include <module_base/kernels/math_kernel_op.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

namespace hsolver {

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagoCG final
{
    // private: accessibility within class is private by default
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;
  public:
    using Func = std::function<void(const ct::Tensor&, ct::Tensor&)>;
    // Constructor need:
    // 1. temporary mock of Hamiltonian "Hamilt_PW"
    // 2. precondition pointer should point to place of precondition array.
    DiagoCG(const std::string& basis_type, const std::string& calculation);
    DiagoCG(
        const std::string& basis_type,
        const std::string& calculation,
        const bool& need_subspace,
        const Func& subspace_func,
        const Real& pw_diag_thr,
        const int& pw_diag_nmax,
        const int& nproc_in_pool);
    
    ~DiagoCG();

    // virtual void init(){};
    // refactor hpsi_info
    // this is the diag() function for CG method
    void diag(const Func& hpsi_func,
              const Func& spsi_func,
              ct::Tensor& psi,
              ct::Tensor& eigen,
              const std::vector<double>& ethr_band,
              const ct::Tensor& prec = {});

  private:
    Device * ctx_ = {};
    /// static variables, used for passing control variables
    /// record for how many bands not have convergence eigenvalues
    int notconv_ = 0;
    /// inside variables and vectors, used by inside functions.
    /// row size for input psi matrix
    int n_band_ = 0;
    /// col size for input psi matrix
    int n_basis_ = 0;
    /// average iteration steps for cg diagonalization
    int avg_iter_ = 0;
    /// threshold for cg diagonalization
    Real pw_diag_thr_ = 1e-5;
    /// maximum iteration steps for cg diagonalization
    int pw_diag_nmax_ = 0;
    /// number of processors in a node
    int nproc_in_pool_ = 0;
    /// basis_type of psi
    std::string basis_type_ = {};
    /// calculation type of ABACUS
    std::string calculation_ = {};

    bool need_subspace_ = false;
    /// A function object that performs the hPsi calculation.
    std::function<void(const ct::Tensor&, ct::Tensor&)> hpsi_func_ = nullptr;
    /// A function object that performs the sPsi calculation.
    std::function<void(const ct::Tensor&, ct::Tensor&)> spsi_func_ = nullptr;
    /// A function object that performs the subspace calculation.
    std::function<void(const ct::Tensor&, ct::Tensor&)> subspace_func_ = nullptr;

    void calc_grad(
        const ct::Tensor& prec,
        ct::Tensor& grad,
        ct::Tensor& hphi,
        ct::Tensor& sphi,
        ct::Tensor& pphi);

    void orth_grad(
        const ct::Tensor& psi, 
        const int& m, 
        ct::Tensor& grad, 
        ct::Tensor& scg,
        ct::Tensor& lagrange);

    void calc_gamma_cg(
        const int& iter,
        const Real& cg_norm, 
        const Real& theta,
        const ct::Tensor& prec,
        const ct::Tensor& scg,
        const ct::Tensor& grad,
        const ct::Tensor& phi_m,
        Real& gg_last,
        ct::Tensor& g0,
        ct::Tensor& cg);

    bool update_psi(
        const ct::Tensor& pphi,
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

    // used in diag() for template replace Hamilt with Hamilt_PW
    void diag_mock(const ct::Tensor& prec,
                   ct::Tensor& psi,
                   ct::Tensor& eigen,
                   const std::vector<double>& ethr_band);

    bool test_exit_cond(const int& ntry, const int& notconv) const;

    using dot_real_op = ModuleBase::dot_real_op<T, Device>;
    const T * one_ = nullptr, * zero_ = nullptr, * neg_one_ = nullptr;
};

} // namespace hsolver

#endif // MODULE_HSOLVER_DIAGO_CG_H_
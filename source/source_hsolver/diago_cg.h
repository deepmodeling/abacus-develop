#ifndef MODULE_HSOLVER_DIAGO_CG_H_
#define MODULE_HSOLVER_DIAGO_CG_H_

#include <vector>

#include <source_base/macros.h>
#include <source_base/kernels/math_kernel_op.h>

#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/hs_operator.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

namespace hsolver {

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagoCG final
{
    // private: accessibility within class is private by default
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, std::complex<double>)
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;
  public:
    // Constructor need:
    // 1. basis type and calculation type of ABACUS
    // 2. diag_comm: the pool over which the subspace matrices are reduced
    DiagoCG(const std::string& basis_type, const std::string& calculation);
    DiagoCG(
        const std::string& basis_type,
        const std::string& calculation,
        const bool& need_subspace,
        const diag_comm_info& diag_comm,
        const Real& pw_diag_thr,
        const int& pw_diag_nmax);

    ~DiagoCG();

    // this is the diag() function for CG method
    // returns avg_iter
    double diag(const HSOperator<T, Device>& op,
                const int ld_psi,
                const int nband,
                const int dim,
                T* psi_in,
                Real* eigenvalue_in,
                const std::vector<double>& ethr_band,
                const Real* prec = nullptr);

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
    double avg_iter_ = 0;
    /// std::vector for iter count of each band
    std::vector<int> iter_band;
    /// threshold for cg diagonalization
    Real pw_diag_thr_ = 1e-5;
    /// maximum iteration steps for cg diagonalization
    int pw_diag_nmax_ = 0;
    /// communicator of the pool sharing the plane waves
    const diag_comm_info diag_comm_;
    /// basis_type of psi
    std::string basis_type_ = {};
    /// calculation type of ABACUS
    std::string calculation_ = {};

    bool need_subspace_ = false;
    /// The H and S operator being diagonalized, set for the duration of diag().
    const HSOperator<T, Device>* op_ = nullptr;

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
    void diag_once(const ct::Tensor& prec,
                   ct::Tensor& psi,
                   ct::Tensor& eigen,
                   const std::vector<double>& ethr_band);

    bool test_exit_cond(const int& ntry, const int& notconv) const;

    /// subspace rotation of the current nband vectors (packed, leading dimension dim)
    void diag_subspace(const T* psi_in, T* psi_out, const int dim, const int nband, const bool S_orth);

    using dot_real_op = ModuleBase::dot_real_op<T, Device>;
    const T * one_ = nullptr, * zero_ = nullptr, * neg_one_ = nullptr;
};

} // namespace hsolver

#endif // MODULE_HSOLVER_DIAGO_CG_H_

#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_base/macros.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"

namespace hsolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolverPW
{
  private:
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in,
              wavefunc* pwf_in,
              const bool initialed_psi_in);

    /// @brief solve function for pw
    /// @param pHamilt interface to hamilt
    /// @param psi reference to psi
    /// @param pes interface to elecstate
    /// @param method_in dav or cg
    /// @param skip_charge
    void solve(hamilt::Hamilt<T, Device>* pHamilt,
               psi::Psi<T, Device>& psi,
               elecstate::ElecState* pes,
               double* out_eigenvalues,
               const std::vector<bool>& is_occupied_in,
               const std::string method_in,
               const std::string calculation_type_in,
               const std::string basis_type_in,
               const bool use_paw_in,
               const bool use_uspp_in,
               const int rank_in_pool_in,
               const int nproc_in_pool_in,
               const int scf_iter_in,
               const bool need_subspace_in,
               const int diag_iter_max_in,
               const double iter_diag_thr_in,
               const bool skip_charge);
    
  protected:
    // diago caller
    void hamiltSolvePsiK(hamilt::Hamilt<T, Device>* hm,
                         psi::Psi<T, Device>& psi,
                         std::vector<Real>& pre_condition,
                         Real* eigenvalue);

    // psi initializer && change k point in psi
    void updatePsiK(hamilt::Hamilt<T, Device>* pHamilt, psi::Psi<T, Device>& psi, const int ik);

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<Real>& h_diag, const int ik, const int npw);

    void output_iterInfo();

    bool initialed_psi = false;

    ModulePW::PW_Basis_K* wfc_basis = nullptr;

    wavefunc* pwf = nullptr;

    int scf_iter = 1; // Start from 1
    bool need_subspace = false;
    int diag_iter_max = 50;
    double iter_diag_thr = 1.0e-2;  // threshold for diagonalization

    std::string method = "none";

  private:
    Device* ctx = {};

    std::string calculation_type = "scf";
    std::string basis_type = "pw";

    bool use_paw = false;
    bool use_uspp = false;

    int rank_in_pool = 0;
    int nproc_in_pool = 1;

    int nspin = 1;

#ifdef USE_PAW
    void paw_func_in_kloop(const int ik);

    void call_paw_cell_set_currentk(const int ik);

    void paw_func_after_kloop(psi::Psi<T, Device>& psi, elecstate::ElecState* pes);
#endif
    
};


} // namespace hsolver

#endif
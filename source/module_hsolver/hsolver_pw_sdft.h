#ifndef HSOLVERPW_SDFT_H
#define HSOLVERPW_SDFT_H
#include "hsolver_pw.h"
#include "module_hamilt_pw/hamilt_stodft/hamilt_sdft_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
namespace hsolver
{
template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolverPW_SDFT : public HSolverPW<T, Device>
{
  public:
    HSolverPW_SDFT(K_Vectors* pkv,
                   ModulePW::PW_Basis_K* wfc_basis_in,
                   wavefunc* pwf_in,
                   Stochastic_WF<T, Device>& stowf,
                   StoChe<double>& stoche,
                   hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto,
                   const std::string calculation_type_in,
                   const std::string basis_type_in,
                   const std::string method_in,
                   const bool use_paw_in,
                   const bool use_uspp_in,
                   const int nspin_in,
                   const int scf_iter_in,
                   const int diag_iter_max_in,
                   const double diag_thr_in,
                   const bool need_subspace_in,
                   const bool initialed_psi_in)
        : HSolverPW<T, Device>(wfc_basis_in,
                               pwf_in,
                               calculation_type_in,
                               basis_type_in,
                               method_in,
                               use_paw_in,
                               use_uspp_in,
                               nspin_in,
                               scf_iter_in,
                               diag_iter_max_in,
                               diag_thr_in,
                               need_subspace_in,
                               initialed_psi_in)
    {
        stoiter.init(pkv, wfc_basis_in, stowf, stoche, p_hamilt_sto);
    }

    void solve(hamilt::Hamilt<T, Device>* pHamilt,
               psi::Psi<T, Device>& psi,
               elecstate::ElecState* pes,
               ModulePW::PW_Basis_K* wfc_basis,
               Stochastic_WF<T, Device>& stowf,
               const int istep,
               const int iter,
               const bool skip_charge);

    Stochastic_Iter<T, Device> stoiter;
};
} // namespace hsolver
#endif
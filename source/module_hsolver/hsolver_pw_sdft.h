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
  protected:
    using Real = typename GetTypeReal<T>::type;

  public:
    HSolverPW_SDFT(K_Vectors* pkv,
                   ModulePW::PW_Basis_K* wfc_basis_in,
                   Stochastic_WF<T, Device>& stowf,
                   StoChe<Real, Device>& stoche,
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
                   const bool need_subspace_in)
        : HSolverPW<T, Device>(wfc_basis_in,
                               calculation_type_in,
                               basis_type_in,
                               method_in,
                               use_paw_in,
                               use_uspp_in,
                               nspin_in,
                               scf_iter_in,
                               diag_iter_max_in,
                               diag_thr_in,
                               need_subspace_in)
    {
        stoiter.init(pkv, wfc_basis_in, stowf, stoche, p_hamilt_sto);
    }

    void solve(hamilt::Hamilt<T, Device>* pHamilt,
               psi::Psi<T, Device>& psi,
               psi::Psi<T>& psi_cpu,
               elecstate::ElecState* pes,
               ModulePW::PW_Basis_K* wfc_basis,
               Stochastic_WF<T, Device>& stowf,
               const int istep,
               const int iter,
               const bool skip_charge);

    Stochastic_Iter<T, Device> stoiter;

  protected:
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
};
} // namespace hsolver
#endif
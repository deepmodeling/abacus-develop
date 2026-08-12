#ifndef HSOLVERPW_SDFT_H
#define HSOLVERPW_SDFT_H
#include "hsolver_pw.h"
#include "source_pw/module_stodft/hamilt_sdft_pw.h"
#include "source_pw/module_stodft/sto_iter.h"
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
                   const bool use_uspp_in,
                   const int nspin_in,
                   const int scf_iter_in,
                   const int diag_iter_max_in,
                   const double diag_thr_in,
                   const bool need_subspace_in,
                   const int nbands_in,
                   const bool diago_smooth_ethr_in,
                   const int pw_diag_ndim_in,
                   const int diag_subspace_in,
                   const int nb2d_in,
                   const bool ks_run_in,
                   const bool all_ks_run_in,
                   const int bndpar_in)
        : HSolverPW<T, Device>(wfc_basis_in,
                               calculation_type_in,
                               basis_type_in,
                               method_in,
                               use_uspp_in,
                               nspin_in,
                               scf_iter_in,
                               diag_iter_max_in,
                               diag_thr_in,
                               need_subspace_in,
                               nbands_in,
                               diago_smooth_ethr_in,
                               pw_diag_ndim_in,
                               diag_subspace_in,
                               nb2d_in),
          ks_run(ks_run_in), all_ks_run(all_ks_run_in), bndpar(bndpar_in)
    {
        stoiter.init(pkv, wfc_basis_in, stowf, stoche, p_hamilt_sto);
    }

    void solve(const UnitCell& ucell,
               hamilt::Hamilt<T, Device>* pHamilt,
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
    const bool ks_run;     // true if the current process runs the KS part of the SDFT calculation
    const bool all_ks_run; // true if every process runs the KS part
    const int bndpar;      // number of band-parallel groups

    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
};
} // namespace hsolver
#endif
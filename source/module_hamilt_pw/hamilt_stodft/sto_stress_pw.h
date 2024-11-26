#ifndef STO_STRESS_PW_H
#define STO_STRESS_PW_H

#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "sto_wf.h"
// qianrui create 2021-6-4

template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Sto_Stress_PW : public Stress_Func<FPTYPE, Device>
{
  public:
    Sto_Stress_PW(){};
    ~Sto_Stress_PW(){};

    // calculate the stress in PW basis
    void cal_stress(ModuleBase::matrix& sigmatot,
                    const elecstate::ElecState& elec,
                    ModulePW::PW_Basis* rho_basis,
                    ModuleSymmetry::Symmetry* p_symm,
                    Structure_Factor* p_sf,
                    K_Vectors* p_kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
                    const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf,
                    const Charge* const chr,
                    pseudopot_cell_vnl* nlpp_in,
                    const UnitCell& ucell_in);

  private:
    void sto_stress_kin(ModuleBase::matrix& sigma,
                        const ModuleBase::matrix& wg,
                        ModuleSymmetry::Symmetry* p_symm,
                        K_Vectors* p_kv,
                        ModulePW::PW_Basis_K* wfc_basis,
                        const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
                        const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf);

    void sto_stress_nl(ModuleBase::matrix& sigma,
                       const ModuleBase::matrix& wg,
                       Structure_Factor* p_sf,
                       ModuleSymmetry::Symmetry* p_symm,
                       K_Vectors* p_kv,
                       ModulePW::PW_Basis_K* wfc_basis,
                       const pseudopot_cell_vnl& nlpp,
                       const UnitCell& ucell,
                       const psi::Psi<std::complex<FPTYPE>, Device>& psi,
                       const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf);
  private:
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using setmem_var_op = base_device::memory::set_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;
};
#endif

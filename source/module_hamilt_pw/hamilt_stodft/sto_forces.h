#ifndef STO_FORCES_H
#define STO_FORCES_H

#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_psi/psi.h"
#include "sto_wf.h"

template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Sto_Forces : public Forces<FPTYPE, Device>
{
  public:
    /* This routine is a driver routine which compute the forces
     * acting on the atoms, the complete forces in plane waves
     * is computed from 4 main parts
     * (1) cal_force_loc: contribution due to local potential.
     * (2) cal_foce_ew: contribution due to ewald potential.
     * (3) cal_force_cc: contributino due to NLCC.
     * (4) cal_nl: contribution due to the non-local pseudopotential.
     * (4) cal_scc: contributino due to incomplete SCF calculation.
     */
    Sto_Forces(const int nat_in) : Forces<FPTYPE, Device>(nat_in){};
    ~Sto_Forces(){};

    void cal_stoforce(ModuleBase::matrix& force,
                      const elecstate::ElecState& elec,
                      ModulePW::PW_Basis* rho_basis,
                      ModuleSymmetry::Symmetry* p_symm,
                      const Structure_Factor* p_sf,
                      K_Vectors* pkv,
                      ModulePW::PW_Basis_K* wfc_basis,
                      const pseudopot_cell_vnl& nlpp,
                      UnitCell& ucell,
                      const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
                      const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf);

  private:
    void cal_sto_force_nl(ModuleBase::matrix& forcenl,
                          const ModuleBase::matrix& wg,
                          K_Vectors* p_kv,
                          ModulePW::PW_Basis_K* wfc_basis,
                          const Structure_Factor* p_sf,
                          const pseudopot_cell_vnl& nlpp,
                          const UnitCell& ucell,
                          const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
                          const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf);

    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;
};

#endif

#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_STRESS_PW_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_STRESS_PW_H

#include "module_elecstate/elecstate.h"
#include "stress_func.h"

template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Stress_PW : public Stress_Func<FPTYPE, Device>
{
  public:
    Stress_PW(const elecstate::ElecState* pelec_in) : pelec(pelec_in){};

    // calculate the stress in PW basis
    void cal_stress(ModuleBase::matrix& smearing_sigmatot,
                    UnitCell& ucell,
                    pseudopot_cell_vnl* nlpp,
                    ModulePW::PW_Basis* rho_basis,
                    ModuleSymmetry::Symmetry* p_symm,
                    Structure_Factor* p_sf,
                    K_Vectors* p_kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const psi::Psi<complex<FPTYPE>, Device>* d_psi_in = nullptr);

  protected:
    // call the vdw stress
    void stress_vdw(ModuleBase::matrix& smearing_sigma,
                    UnitCell& ucell); // force and stress calculated in vdw together.

    // the stress from the non-local pseudopotentials in uspp
    // which is due to the dependence of the Q function on the atomic position
    void stress_us(ModuleBase::matrix& sigma,
                   ModulePW::PW_Basis* rho_basis,
                   pseudopot_cell_vnl* ppcell_in,
                   const UnitCell& ucell); // nonlocal part of uspp in PW basis

    const elecstate::ElecState* pelec = nullptr;
};
#endif

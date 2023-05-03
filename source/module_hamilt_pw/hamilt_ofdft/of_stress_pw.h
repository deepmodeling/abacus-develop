#ifndef OF_STRESS_PW_H
#define OF_STRESS_PW_H

#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"

class OF_Stress_PW : public Stress_Func<double>
{
  public:
    OF_Stress_PW(const elecstate::ElecState* pelec_in, ModulePW::PW_Basis* rhopw_in)
        : pelec(pelec_in), rhopw(rhopw_in){};

    // calculate the stress in OFDFT
    void cal_stress(ModuleBase::matrix& sigmatot,
                    ModuleBase::matrix& kinetic_stress,
                    UnitCell& ucell,
                    ModuleSymmetry::Symmetry& symm,
                    Structure_Factor& sf,
                    K_Vectors& kv,
                    ModulePW::PW_Basis_K* wfc_basis = nullptr,
                    const psi::Psi<complex<double>>* psi_in = nullptr);

  protected:
    // call the vdw stress
    void stress_vdw(ModuleBase::matrix& smearing_sigma,
                    UnitCell& ucell); // force and stress calculated in vdw together.

    const elecstate::ElecState* pelec = nullptr;
    ModulePW::PW_Basis* rhopw = nullptr;
};
#endif

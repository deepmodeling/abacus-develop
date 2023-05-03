#ifndef STRESS_PW_H
#define STRESS_PW_H

#include "module_elecstate/elecstate.h"
#include "stress_func.h"

template <typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Stress_PW : public Stress_Func<FPTYPE, Device>
{
  public:
    Stress_PW(const elecstate::ElecState* pelec_in) : pelec(pelec_in){};

    // calculate the stress in PW basis
    void cal_stress(ModuleBase::matrix& smearing_sigmatot,
                    UnitCell& ucell,
                    ModulePW::PW_Basis* rho_basis,
                    ModuleSymmetry::Symmetry& symm,
                    Structure_Factor& sf,
                    K_Vectors& kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const psi::Psi<complex<FPTYPE>>* psi_in = nullptr,
                    const psi::Psi<complex<FPTYPE>, Device>* d_psi_in = nullptr);

  protected:
    // call the vdw stress
    void stress_vdw(ModuleBase::matrix& smearing_sigma,
                    UnitCell& ucell); // force and stress calculated in vdw together.

    const elecstate::ElecState* pelec = nullptr;
};
#endif

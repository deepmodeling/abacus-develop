#ifndef STO_STRESS_PW_H
#define STO_STRESS_PW_H

#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "sto_wf.h"
// qianrui create 2021-6-4

class Sto_Stress_PW : public Stress_Func<double>
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
                    const psi::Psi<complex<double>>* psi_in,
                    Stochastic_WF& stowf,
                    const Charge* const chr);

  private:
    void sto_stress_kin(ModuleBase::matrix& sigma,
                        const ModuleBase::matrix& wg,
                        ModuleSymmetry::Symmetry* p_symm,
                        K_Vectors* p_kv,
                        ModulePW::PW_Basis_K* wfc_basis,
                        const psi::Psi<complex<double>>* psi_in,
                        Stochastic_WF& stowf);

    void sto_stress_nl(ModuleBase::matrix& sigma,
                       const ModuleBase::matrix& wg,
                       Structure_Factor* p_sf,
                       ModuleSymmetry::Symmetry* p_symm,
                       K_Vectors* p_kv,
                       ModulePW::PW_Basis_K* wfc_basis,
                       const psi::Psi<complex<double>>* psi_in,
                       Stochastic_WF& stowf);
};
#endif

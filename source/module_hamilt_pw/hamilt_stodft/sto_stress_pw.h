#ifndef STO_STRESS_PW_H
#define STO_STRESS_PW_H

#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "sto_wf.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/module_charge/charge.h"
//qianrui create 2021-6-4

class Sto_Stress_PW:public Stress_Func<double>
{
	public :
	
	Sto_Stress_PW (){};
	~Sto_Stress_PW (){};

	//calculate the stress in PW basis
    void cal_stress(ModuleBase::matrix& sigmatot,
                    const ModuleBase::matrix& wg,
                    ModulePW::PW_Basis* rho_basis,
                    ModuleSymmetry::Symmetry& symm,
                    Structure_Factor& sf,
                    K_Vectors& kv,
                    ModulePW::PW_Basis_K* wfc_basis,
                    const psi::Psi<complex<double>>* psi_in,
                    Stochastic_WF& stowf,
                    const Charge* const chr);

  private:
    void sto_stress_kin(ModuleBase::matrix& sigma,
                        const ModuleBase::matrix& wg,
                        ModuleSymmetry::Symmetry& symm,
                        K_Vectors& kv,
                        ModulePW::PW_Basis_K* wfc_basis,
                        const psi::Psi<complex<double>>* psi_in,
                        Stochastic_WF& stowf);

    void sto_stress_nl(ModuleBase::matrix& sigma,
                       const ModuleBase::matrix& wg,
                       ModuleSymmetry::Symmetry& symm,
                       K_Vectors& kv,
                       ModulePW::PW_Basis_K* wfc_basis,
                       const psi::Psi<complex<double>>* psi_in,
                       Stochastic_WF& stowf);
};
#endif

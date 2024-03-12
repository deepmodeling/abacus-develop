#ifndef STO_FORCES_H
#define STO_FORCES_H

#include "module_hamilt_pw/hamilt_pwdft/forces.h"
#include "module_psi/psi.h"
#include "sto_wf.h"

class Sto_Forces : public Forces<double>
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
    Sto_Forces(const int nat_in):Forces(nat_in){};
    ~Sto_Forces(){};

    void cal_stoforce(ModuleBase::matrix& force,
                      const elecstate::ElecState& elec,
                      ModulePW::PW_Basis* rho_basis,
                      ModuleSymmetry::Symmetry* p_symm,
                      Structure_Factor* p_sf,
                      K_Vectors* pkv,
                      ModulePW::PW_Basis_K* wfc_basis,
                      const psi::Psi<std::complex<double>>* psi_in,
                      Stochastic_WF& stowf);

  private:
    void cal_sto_force_nl(ModuleBase::matrix& forcenl,
                          const ModuleBase::matrix& wg,
                          K_Vectors* p_kv,
                          ModulePW::PW_Basis_K* wfc_basis,
                          const psi::Psi<complex<double>>* psi_in,
                          Stochastic_WF& stowf);
};

#endif

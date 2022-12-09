#ifndef STO_FORCES_H
#define STO_FORCES_H

#include "forces.h"
#include "./sto_wf.h"
#include "module_psi/psi.h"

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
    Sto_Forces(){};
    ~Sto_Forces(){};

    void init(ModuleBase::matrix& matrix, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>* psi_in,Stochastic_WF& stowf, const Charge* const chr);

private:
    void cal_sto_force_nl(ModuleBase::matrix& forcenl, const ModuleBase::matrix& wg, const psi::Psi<complex<double>>* psi_in, Stochastic_WF& stowf);

};

#endif

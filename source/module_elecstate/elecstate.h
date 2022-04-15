#ifndef ELECSTATE_H
#define ELECSTATE_H

#include "module_psi/psi.h"
#include "module_hamilt/matrixblock.h"
#include "src_pw/charge.h"

namespace ModuleElecS
{

class ElecState
{
    public:
    virtual void init(Charge* chg_in
    /*const Basis &basis, const Cell &cell*/) = 0;
    
    //return current electronic density rho, as a input for constructing Hamiltonian
    virtual const ModuleHamilt::MatrixBlock<double> getRho()const = 0;
    
    //calculate electronic charge density on grid points or density matrix in real space
    //the consequence charge density rho saved into rho_out, preparing for charge mixing. 
    virtual void updateRhoK(const ModulePsi::Psi<std::complex<double>> &psi)=0;
    virtual void updateRhoK(const ModulePsi::Psi<double> &psi){return;}
    
    //update charge density for next scf step
    virtual void getNewRho() = 0;

};

}
#endif
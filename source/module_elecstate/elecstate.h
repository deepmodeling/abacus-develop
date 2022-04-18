#ifndef ELECSTATE_H
#define ELECSTATE_H

#include "module_hamilt/matrixblock.h"
#include "module_psi/psi.h"
#include "src_pw/charge.h"

namespace elecstate
{

class ElecState
{
  public:
    virtual void init(Charge *chg_in
                      /*const Basis &basis, const Cell &cell*/)
        = 0;

    // return current electronic density rho, as a input for constructing Hamiltonian
    virtual const hamilt::MatrixBlock<double> getRho() const = 0;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    virtual void updateRhoK(const psi::Psi<std::complex<double>> &psi) = 0;
    virtual void updateRhoK(const psi::Psi<double> &psi)
    {
        return;
    }

    // update charge density for next scf step
    // in this function, 1. input rho for construct Hamilt and 2. calculated rho from Psi will mix to 3. new charge
    // density rho among these rho,
    // 1. input rho would be store to file for restart
    // 2. calculated rho should be near with input rho when convergence has achieved
    // 3. new rho should be input rho for next scf step.
    virtual void getNewRho() = 0;
};

} // namespace elecstate
#endif
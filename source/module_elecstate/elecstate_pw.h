#ifndef ELECSTATEPW_H
#define ELECSTATEPW_H

#include "elecstate.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"

namespace elecstate
{

class ElecStatePW : public ElecState
{
  public:
    void init(Charge* chg_in
              /*const Basis &basis, const Cell &cell*/) override;

    // return current electronic density rho, as a input for constructing Hamiltonian
    const hamilt::MatrixBlock<double> getRho() const override;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    void updateRhoK(const psi::Psi<std::complex<double>>& psi) override;

    // update charge density for next scf step
    void getNewRho() override;

  private:
    Charge* pchg;
};

} // namespace elecstate

#endif
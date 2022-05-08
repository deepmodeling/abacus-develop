#ifndef ELECSTATEPW_H
#define ELECSTATEPW_H

#include "elecstate.h"
#include "src_pw/pw_basis.h"

namespace elecstate
{

class ElecStatePW : public ElecState
{
  public:
    ElecStatePW(const PW_Basis* basis_in, Charge* chg_in, int nbands_in) : basis(basis_in)
    {
        init(chg_in, basis_in->Klist->nks, nbands_in);
    }
    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override;
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

  private:
    const PW_Basis* basis;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    void updateRhoK(const psi::Psi<std::complex<double>>& psi); // override;
    // sum over all pools for rho and ebands
    void parallelK();
    // calcualte rho for each k
    void rhoBandK(const psi::Psi<std::complex<double>>& psi);
};

} // namespace elecstate

#endif
#ifndef ELECSTATELCAO_H
#define ELECSTATELCAO_H

#include "elecstate.h"
#include "src_lcao/LCAO_hamilt.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"

namespace elecstate
{

class ElecStateLCAO : public ElecState
{
  public:
    ElecStateLCAO(Charge* chg_in,
                  int nks_in,
                  int nbands_in,
                  Local_Orbital_Charge* loc_in,
                  LCAO_Hamilt* uhm_in,
                  Local_Orbital_wfc* lowf_in)
    {
        init(chg_in, nks_in, nbands_in);
        this->loc = loc_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
    }
    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override;
    virtual void psiToRho(const psi::Psi<double>& psi) override;
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

    static int out_wfc_lcao;

  private:
    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    // void updateRhoK(const psi::Psi<std::complex<double>>& psi) ;//override;
    // sum over all pools for rho and ebands
    // void parallelK();
    // calcualte rho for each k
    // void rhoBandK(const psi::Psi<std::complex<double>>& psi);

    template <typename T> void cal_dm(const ModuleBase::matrix& wg, const psi::Psi<T>& wfc, psi::Psi<T>& dm);

    Local_Orbital_Charge* loc = nullptr;
    LCAO_Hamilt* uhm = nullptr;
    Local_Orbital_wfc* lowf = nullptr;
};

} // namespace elecstate

#endif
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
    ElecStateLCAO(Charge* chg_in = nullptr,
                        const K_Vectors* klist_in = nullptr,
                        int nks_in = 1,
                        int nbands_in = 1,
                        Local_Orbital_Charge* loc_in = nullptr,
                        LCAO_Hamilt* uhm_in = nullptr,
                        Local_Orbital_wfc* lowf_in = nullptr)
    {
        init(chg_in, klist_in, nks_in, nbands_in);
        this->loc = loc_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->classname = "ElecStateLCAO";
    }
    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override;
    virtual void psiToRho(const psi::Psi<double>& psi) override;
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

    virtual void print_psi(const psi::Psi<double>& psi_in)override;
    virtual void print_psi(const psi::Psi<std::complex<double>>& psi_in)override;

    static int out_wfc_lcao;
    static bool need_psi_grid;

  protected:
    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    // void updateRhoK(const psi::Psi<std::complex<double>>& psi) ;//override;
    // sum over all pools for rho and ebands
    // void parallelK();
    // calcualte rho for each k
    // void rhoBandK(const psi::Psi<std::complex<double>>& psi);

    Local_Orbital_Charge* loc = nullptr;
    LCAO_Hamilt* uhm = nullptr;
    Local_Orbital_wfc* lowf = nullptr;
};

} // namespace elecstate

#endif
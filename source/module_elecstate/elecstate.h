#ifndef ELECSTATE_H
#define ELECSTATE_H

#include "module_psi/psi.h"
#include "src_pw/charge.h"
#include "src_pw/klist.h"

namespace elecstate
{

class ElecState
{
  public:
    ElecState(){};
    virtual ~ElecState(){};
    virtual void init(Charge *chg_in, // pointer for class Charge
                      const K_Vectors *klist_in,
                      int nk_in, // number of k points
                      int nb_in) // number of bands
    {
        this->charge = chg_in;
        this->klist = klist_in;
        this->ekb.create(nk_in, nb_in);
        this->wg.create(nk_in, nb_in);
    }

    // return current electronic density rho, as a input for constructing Hamiltonian
    virtual const double *getRho(int spin) const;

    // allocate memory for charge
    virtual void allocateRho(const int &nspin_in, const int &nrxx_in, const int &ngmc_in);

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    virtual void psiToRho(const psi::Psi<std::complex<double>> &psi)
    {
        return;
    }
    virtual void psiToRho(const psi::Psi<double> &psi)
    {
        return;
    }
    // virtual void updateRhoK(const psi::Psi<std::complex<double>> &psi) = 0;
    // virtual void updateRhoK(const psi::Psi<double> &psi)=0

    // update charge density for next scf step
    // in this function, 1. input rho for construct Hamilt and 2. calculated rho from Psi will mix to 3. new charge
    // density rho among these rho,
    // 1. input rho would be store to file for restart
    // 2. calculated rho should be near with input rho when convergence has achieved
    // 3. new rho should be input rho for next scf step.
    virtual void getNewRho()
    {
        return;
    }

    // calculate wg from ekb
    virtual void calculate_weights();

    virtual void print_psi(const psi::Psi<double>& psi_in)
    {
        return;
    }
    virtual void print_psi(const psi::Psi<std::complex<double>>& psi_in)
    {
        return;
    }

    // pointer to charge density
    Charge *charge = nullptr;
    // pointer to k points lists
    const K_Vectors* klist = nullptr;
    // energy for sum of electrons
    double eband = 0.0;
    // Fermi energy
    double ef = 0.0;
    // correction energy for metals
    double demet = 0.0;

    // band energy at each k point, each band.
    ModuleBase::matrix ekb;
    // occupation weight for each k-point and band
    ModuleBase::matrix wg;

    std::string classname = "none";

  protected:
    // calculate ebands for all k points and all occupied bands
    void calEBand();

    // print and check for band energy and occupations
    void print_band(const int &ik, const int &printe, const int &iter);
};

} // namespace elecstate
#endif

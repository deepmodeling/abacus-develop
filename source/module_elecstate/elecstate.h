#ifndef ELECSTATE_H
#define ELECSTATE_H

#include "module_psi/psi.h"
#include "src_pw/charge.h"
#include "src_pw/klist.h"
#include "potentials/potential_new.h"

namespace elecstate
{

class ElecState
{
  public:
    ElecState(){};
    ElecState(Charge* charge_in){this->charge = charge_in;};
    virtual ~ElecState()
    {
        if(this->pot != nullptr) 
        {
            delete this->pot;
            this->pot = nullptr;
        }
    };
    void init_ks(Charge *chg_in, // pointer for class Charge
                      const K_Vectors *klist_in,
                      int nk_in); // number of k points

    // return current electronic density rho, as a input for constructing Hamiltonian
    virtual const double *getRho(int spin) const;

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
    // use occupied weights from INPUT and skip calculate_weights
    void fixed_weights(const double * const ocp_kb);
    // if nupdown is not 0(TWO_EFERMI case), 
    // nelec_spin will be fixed and weights will be constrained 
    void init_nelec_spin();
    //used to record number of electrons per spin index
    //for NSPIN=2, it will record number of spin up and number of spin down
    //for NSPIN=4, it will record total number, magnetization for x, y, z direction  
    std::vector<double> nelec_spin;

    //calculate nbands and 
    void cal_nbands();

    virtual void print_psi(const psi::Psi<double>& psi_in)
    {
        return;
    }
    virtual void print_psi(const psi::Psi<std::complex<double>>& psi_in)
    {
        return;
    }

    void init_scf(const int istep, const ModuleBase::ComplexMatrix& strucfac);

    // pointer to potential
    Potential* pot = nullptr;
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

    // print and check for band energy and occupations
    void print_band(const int &ik, const int &printe, const int &iter);

  protected:
    // calculate ebands for all k points and all occupied bands
    void calEBand();

  private:
    bool skip_weights = false;
};

} // namespace elecstate
#endif

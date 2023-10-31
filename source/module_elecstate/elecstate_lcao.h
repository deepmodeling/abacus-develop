#ifndef ELECSTATELCAO_H
#define ELECSTATELCAO_H

#include "elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_elecstate/module_dm/density_matrix.h"

namespace elecstate
{
template <typename TK>
class ElecStateLCAO : public ElecState
{
  public:
    ElecStateLCAO(){} // will be called by ElecStateLCAO_TDDFT
    ElecStateLCAO(Charge* chg_in ,
                  const K_Vectors* klist_in ,
                  int nks_in,
                  Local_Orbital_Charge* loc_in ,
                  LCAO_Hamilt* uhm_in ,
                  Local_Orbital_wfc* lowf_in ,
                  ModulePW::PW_Basis* rhopw_in ,
                  ModulePW::PW_Basis_Big* bigpw_in )
    {
        init_ks(chg_in, klist_in, nks_in, rhopw_in, bigpw_in);
        this->loc = loc_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->classname = "ElecStateLCAO";
    }

    virtual ~ElecStateLCAO()
    {
        if (this->DM != nullptr)
            delete this->DM;
    }

    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<TK>& psi) override;
    //virtual void psiToRho(const psi::Psi<double>& psi) override;
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

    virtual void print_psi(const psi::Psi<TK>& psi_in, const int istep = -1) override;
    //virtual void print_psi(const psi::Psi<std::complex<double>>& psi_in, const int istep = -1) override;

    // initial density matrix
    void init_DM(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin);
    DensityMatrix<TK,double>* get_DM() const { return const_cast<DensityMatrix<TK,double>*>(this->DM); } 
    static int out_wfc_lcao;
    static int out_wfc_flag;
    static bool need_psi_grid;

    double get_spin_constrain_energy() override;

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
    DensityMatrix<TK,double>* DM = nullptr;
    
};

} // namespace elecstate

#endif
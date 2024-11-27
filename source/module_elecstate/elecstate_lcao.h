#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ELECSTATE_ELECSTATE_LCAO_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ELECSTATE_ELECSTATE_LCAO_H

#include "elecstate.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

#include <vector>

namespace elecstate
{
template <typename TK>
class ElecStateLCAO : public ElecState
{
  public:
    ElecStateLCAO()
    {
    } // will be called by ElecStateLCAO_TDDFT
    ElecStateLCAO(Charge* chg_in,
                  const K_Vectors* klist_in,
                  int nks_in,
                  Gint_Gamma* gint_gamma_in, // mohan add 2024-04-01
                  Gint_k* gint_k_in,         // mohan add 2024-04-01
                  ModulePW::PW_Basis* rhopw_in,
                  ModulePW::PW_Basis_Big* bigpw_in)
    {
        init_ks(chg_in, klist_in, nks_in, rhopw_in, bigpw_in);
        this->gint_gamma = gint_gamma_in; // mohan add 2024-04-01
        this->gint_k = gint_k_in;         // mohan add 2024-04-01
        this->classname = "ElecStateLCAO";
    }

    virtual ~ElecStateLCAO()
    {
        if (this->DM != nullptr)
        {
            delete this->DM;
        }
    }

    // void init(Charge* chg_in):charge(chg_in){} override;

    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<TK>& psi) override;
    // virtual void psiToRho(const psi::Psi<double>& psi) override;
    //  return current electronic density rho, as a input for constructing Hamiltonian
    //  const double* getRho(int spin) const override;
    virtual void cal_tau(const psi::Psi<TK>& psi) override;

    // update charge density for next scf step
    // void getNewRho() override;

    // initial density matrix
    void init_DM(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin);
    DensityMatrix<TK, double>* get_DM() const
    {
        return const_cast<DensityMatrix<TK, double>*>(this->DM);
    }
    static int out_wfc_lcao;
    static bool need_psi_grid;

    double get_spin_constrain_energy() override;

#ifdef __PEXSI
    // use for pexsi

    /**
     * @brief calculate electronic charge density from pointers of density matrix calculated by pexsi
     * @param pexsi_DM: pointers of density matrix (DMK) calculated by pexsi
     * @param pexsi_EDM: pointers of energy-weighed density matrix (EDMK) calculated by pexsi, needed by MD, will be
     * stored in DensityMatrix::pexsi_EDM
     */
    void dmToRho(std::vector<TK*> pexsi_DM, std::vector<TK*> pexsi_EDM);
#endif

    DensityMatrix<TK, double>* DM = nullptr;

  protected:
    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    // void updateRhoK(const psi::Psi<std::complex<double>>& psi) ;//override;
    // sum over all pools for rho and ebands
    // void parallelK();
    // calcualte rho for each k
    // void rhoBandK(const psi::Psi<std::complex<double>>& psi);

    Gint_Gamma* gint_gamma = nullptr; // mohan add 2024-04-01
    Gint_k* gint_k = nullptr;         // mohan add 2024-04-01
};

template <typename TK>
int ElecStateLCAO<TK>::out_wfc_lcao = 0;

template <typename TK>
bool ElecStateLCAO<TK>::need_psi_grid = true;

} // namespace elecstate

#endif

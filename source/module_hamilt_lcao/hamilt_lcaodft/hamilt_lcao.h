#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/module_dm/density_matrix.h"

namespace hamilt
{

// template first for type of k space H matrix elements
// template second for type of temporary matrix, gamma_only fix-gamma-matrix + S-gamma, multi-k fix-Real + S-Real
template <typename TK, typename TR>
class HamiltLCAO : public Hamilt<TK>
{
  public:
    /**
     * @brief Constructor of Hamiltonian for LCAO base
     * HR and SR will be allocated with Operators
    */
    HamiltLCAO(Gint_Gamma* GG_in,
               Gint_k* GK_in,
               LCAO_gen_fixedH* genH_in,
               LCAO_Matrix* LM_in,
               Local_Orbital_Charge* loc_in,
               elecstate::Potential* pot_in,
               const K_Vectors& kv_in,
               elecstate::DensityMatrix<TK,double>* DM_in);
    /**
     * @brief Constructor of vacuum Operators, only HR and SR will be initialed as empty HContainer
    */
    HamiltLCAO(LCAO_Matrix* LM_in, const K_Vectors& kv_in);

    ~HamiltLCAO()
    {
        if (this->ops != nullptr)
        {
            delete this->ops;
        }
        delete this->hR;
        delete this->sR;
    };

    /// get pointer of Operator<TK> ops
    Operator<TK>*& getOperator();
    /// get hk-pointer of std::vector<TK>, the return will be LM->Hloc or LM->Hloc2
    std::vector<TK>& getHk(LCAO_Matrix* LM);
    /// get sk-pointer of std::vector<TK>, the return will be this->sk
    std::vector<TK>& getSk(LCAO_Matrix* LM);
    /// get HR pointer of *this->hR, which is a HContainer<TR> and contains H(R)
    HContainer<TR>*& getHR();
    /// get SR pointer of *this->sR, which is a HContainer<TR> and contains S(R)
    HContainer<TR>*& getSR();
    /// refresh the status of HR
    void refresh() override;

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik) override;

    /**
     * @brief special for LCAO, update SK only
     * 
     * @param ik target K point
     * @param hk_type 0: SK is row-major, 1: SK is collumn-major
     * @return void
    */
    void updateSk(const int ik, LCAO_Matrix* LM_in, const int hk_type=0);

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    // void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<TK>& hk_in, MatrixBlock<TK>& sk_in) override;

  private:
    const K_Vectors *kv = nullptr;

    // Real space Hamiltonian
    HContainer<TR>* hR = nullptr;
    HContainer<TR>* sR = nullptr;

    // special case for NSPIN=2 , data of HR should be separated into two parts
    // save them in this->hRS2;
    std::vector<TR> hRS2;
    int refresh_times = 0;

    // sk and hk will be refactored to HamiltLCAO later
    //std::vector<TK> sk;
    //std::vector<TK> hk;
};

} // namespace hamilt

#endif
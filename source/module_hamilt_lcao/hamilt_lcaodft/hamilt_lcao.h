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

namespace hamilt
{

// template first for type of k space H matrix elements
// template second for type of temporary matrix, gamma_only fix-gamma-matrix + S-gamma, multi-k fix-Real + S-Real
template <typename T>
class HamiltLCAO : public Hamilt<std::complex<double>>
{
  public:
    HamiltLCAO(Gint_Gamma* GG_in,
               LCAO_gen_fixedH* genH_in,
               LCAO_Matrix* LM_in,
               Local_Orbital_Charge* loc_in,
               elecstate::Potential* pot_in,
               const K_Vectors& kv_in);

    HamiltLCAO(Gint_k* GK_in,
               LCAO_gen_fixedH* genH_in,
               LCAO_Matrix* LM_in,
               Local_Orbital_Charge* loc_in,
               elecstate::Potential* pot_in,
               const K_Vectors& kv_in);

    ~HamiltLCAO()
    {
        if (this->ops != nullptr)
        {
            delete this->ops;
        }
        if (this->opsd != nullptr)
        {
            delete this->opsd;
        }
    };

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik) override;

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    // void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<T>& hk_in, MatrixBlock<T>& sk_in) override;

  private:
    K_Vectors kv;
};

} // namespace hamilt

#endif
#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "hamilt.h"
#include "src_lcao/LCAO_gen_fixedH.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_gint/gint_gamma.h"
#include "module_gint/gint_k.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"

namespace hamilt
{

// template first for type of k space H matrix elements
// template second for type of temporary matrix, gamma_only fix-gamma-matrix + S-gamma, multi-k fix-Real + S-Real
template <typename T> class HamiltLCAO : public Hamilt
{
  public:
    HamiltLCAO(
      Gint_Gamma* GG_in, 
      LCAO_gen_fixedH* genH_in, 
      LCAO_Matrix* LM_in, 
      Local_Orbital_Charge* loc_in);

    HamiltLCAO(
      Gint_k* GK_in, 
      LCAO_gen_fixedH* genH_in, 
      LCAO_Matrix* LM_in, 
      Local_Orbital_Charge* loc_in);

    ~HamiltLCAO(){
      if(this->ops!= nullptr)
      {
        delete this->ops;
      }
      if(this->opsd!= nullptr)
      {
        delete this->opsd;
      }
    };

    // for target K point, update consequence of hPsi() and matrix()
    virtual void updateHk(const int ik);

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    //void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in) override;
};

} // namespace hamilt

#endif
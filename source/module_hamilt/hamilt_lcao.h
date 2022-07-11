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
      LCAO_Hamilt* uhm_in,
      Local_Orbital_wfc* lowf_in,
      Local_Orbital_Charge* loc_in)
    {
        this->GG = GG_in;
        this->genH = genH_in;
        this->LM = LM_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->loc = loc_in;
        this->classname = "HamiltLCAO";
    }
    HamiltLCAO(
      Gint_k* GK_in, 
      LCAO_gen_fixedH* genH_in, 
      LCAO_Matrix* LM_in, 
      LCAO_Hamilt* uhm_in,
      Local_Orbital_wfc* lowf_in,
      Local_Orbital_Charge* loc_in)
    {
        this->GK = GK_in;
        this->genH = genH_in;
        this->LM = LM_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->loc = loc_in;
        this->classname = "HamiltLCAO";
    }
    ~HamiltLCAO(){
      if(allocated_smatrix) delete[] smatrix_k;
    };

    void constructHamilt() override;

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    //void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in) override;

  private:
    // specific code for matrix()
    void getMatrix(MatrixBlock<T>& hk_in, MatrixBlock<T>& sk_in);

    // there are H and S matrix for each k point in reciprocal space
    // type double for gamma_only case, type complex<double> for multi-k-points case
    T* hmatrix_k = nullptr;
    T* smatrix_k = nullptr;

    // temporary class members
    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH* genH = nullptr;

    LCAO_Matrix* LM = nullptr;

    LCAO_Hamilt* uhm = nullptr;

    Local_Orbital_wfc* lowf;

    Local_Orbital_Charge* loc;

    //only used for Gamma_only case
    bool allocated_smatrix = false;
};

} // namespace hamilt

#endif
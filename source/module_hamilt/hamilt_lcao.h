#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "hamilt.h"
#include "src_lcao/LCAO_gen_fixedH.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/gint_gamma.h"
#include "src_lcao/gint_k.h"

namespace hamilt
{

// memory for storing Hamiltonian matrix and overlap matrix for one k point
template <typename T> class LocalMatrix
{
  public:
    std::vector<T> hloc;
    std::vector<T> sloc;
    size_t size = 1;
    void resize(size_t size_in, T value)
    {
        size = size_in;
        hloc.resize(size, value);
        sloc.resize(size, value);
    }
    T* getH()
    {
        return hloc.data();
    }
    T* getS()
    {
        return sloc.data();
    }
};

// template first for type of k space H matrix elements
// template second for type of temporary matrix, gamma_only fix-gamma-matrix + S-gamma, multi-k fix-Real + S-Real
template <typename T, typename T1> class HamiltLCAO : public Hamilt
{
  public:
    HamiltLCAO(Gint_Gamma* GG_in, Gint_k* GK_in, LCAO_gen_fixedH* genH_in, LCAO_Matrix* LM_in)
    {
        this->GG = GG_in;
        this->GK = GK_in;
        this->genH = genH_in;
        this->LM = LM_in;
    }
    //~HamiltLCAO();

    // construct Hamiltonian matrix with inputed electonic density
    void constructHamilt(const int iter, const MatrixBlock<double> rho) override;

    // for target K point, update consequence of hPsi() and matrix()
    void updateHk(const int ik) override;

    // core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const psi::Psi<std::complex<double>>& psi, psi::Psi<std::complex<double>>& hpsi) const override
    {
        // should be updated for iterative diagonalization method
        return;
    };

    // core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    // not used in PW base
    //void matrix(MatrixBlock<std::complex<double>> &hk_in, MatrixBlock<std::complex<double>> &sk_in) override;
    void matrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in) override;

  private:
    void ch_mock();
    void hk_fixed_mock(const int ik);
    void hk_update_mock(const int ik);

    // specific code for matrix()
    void getMatrix(MatrixBlock<T>& hk_in, MatrixBlock<T>& sk_in);

    // there are H and S matrix for each k point in reciprocal space
    // type double for gamma_only case, type complex<double> for multi-k-points case
    LocalMatrix<T> kM;
    // there are H and S matrix for fixed T+VNL and overlap terms in real space
    // type double for nspin<4, type complex<double> for nspin==4
    LocalMatrix<T1> fixedRealM;

    // control if fixed Matrix should be construct
    static bool isFixedDone;

    void constructFixedReal();
    void constructUpdateReal();

    // temporary class members
    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH* genH = nullptr;

    LCAO_Matrix* LM = nullptr;
};
template <typename T, typename T1> bool HamiltLCAO<T, T1>::isFixedDone = false;

} // namespace hamilt

#endif
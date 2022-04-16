#ifndef DIAGCG_H
#define DIAGCG_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"

#if ((defined __CUDA) || (defined __ROCM))

#ifdef __CUDA
#include "src_pw/hamilt_pw.cuh"
#else
#include "src_pw/hamilt_pw_hip.h"
#endif

#else
#include "src_pw/hamilt_pw.h"
#endif
#include "src_pw/pw_basis.h"

namespace hsolver
{

class DiagoCG : public DiagH
{
  public:
    // Constructor need:
    // 1. temporary mock of Hamiltonian "Hamilt_PW"
    // 2. precondition pointer should point to place of precondition array.
    DiagoCG(Hamilt_PW *hpw_in, const double *precondition_in);
    ~DiagoCG();

    // virtual void init(){};

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in) override;

  private:
    /// static variables, used for passing control variables
    /// if eigenvalue and eigenvectors should be reordered after diagonalization, it is always be true.
    bool reorder = false;
    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    /// temp operator pointer
    Hamilt_PW *hpw = nullptr;

    int test_cg = 0;

    /// inside variables and vectors, used by inside functions.
    /// row size for input psi matrix
    int n_band = 0;
    /// col size for input psi matrix
    int dmx = 0;
    /// non-zero col size for inputted psi matrix
    int dim = 0;
    /// precondition for cg diag
    const double *precondition = nullptr;
    /// eigenvalue results
    double *eigenvalue = nullptr;

    /// temp vector for S|psi> for one band, size dim
    std::vector<std::complex<double>> sphi;
    /// temp vector for , size dim
    std::vector<std::complex<double>> scg;
    /// temp vector for H|psi> for one band, size dim
    std::vector<std::complex<double>> hphi;
    /// temp vector for , size dim
    std::vector<std::complex<double>> gradient;
    /// temp vector for , size dim
    std::vector<std::complex<double>> cg;
    /// temp vector for , size dim
    std::vector<std::complex<double>> g0;
    /// temp vector for store psi in sorting with eigenvalues, size dim
    std::vector<std::complex<double>> pphi;
    /// temp vector for matrix eigenvector * vector S|psi> , size m_band
    std::vector<std::complex<double>> lagrange;
    /// temp vector for new psi for one band, size dim
    std::vector<std::complex<double>> phi_m;

    void calculate_gradient();

    void orthogonal_gradient(const psi::Psi<std::complex<double>> &eigenfunction, const int m);

    void calculate_gamma_cg(const int iter, double &gg_last, const double &cg0, const double &theta);

    bool update_psi(double &cg_norm, double &theta, double &eigenvalue);

    void schmit_orth(const int &m, const psi::Psi<std::complex<double>> &psi);

    // used in diag() for template replace Hamilt with Hamilt_PW
    void diag_mock(psi::Psi<std::complex<double>> &phi, double *eigenvalue_in);
};

} // namespace hsolver
#endif

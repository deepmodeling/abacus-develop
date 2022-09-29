#ifndef DIAGCG_H
#define DIAGCG_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "src_pw/structure_factor.h"

#include "module_psi/include/types.h"
#include "module_psi/include/device.h"
#include "module_psi/include/memory.h"

namespace hsolver
{

template<typename FPTYPE = double, typename Device = psi::DEVICE_CPU>
class DiagoCG : public DiagH
{
  public:
    // Constructor need:
    // 1. temporary mock of Hamiltonian "Hamilt_PW"
    // 2. precondition pointer should point to place of precondition array.
    DiagoCG(const FPTYPE *precondition_in);
    ~DiagoCG();

    // virtual void init(){};
    // refactor hpsi_info
    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &psi, FPTYPE *eigenvalue_in);

  private:
    /// static variables, used for passing control variables
    /// if eigenvalue and eigenvectors should be reordered after diagonalization, it is always be true.
    bool reorder = false;
    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    int test_cg = 0;

    /// inside variables and vectors, used by inside functions.
    /// row size for input psi matrix
    int n_band = 0;
    /// col size for input psi matrix
    int dmx = 0;
    /// non-zero col size for inputted psi matrix
    int dim = 0;
    /// precondition for cg diag
    const FPTYPE *precondition = nullptr;
    /// eigenvalue results
    FPTYPE *eigenvalue = nullptr;


    /// temp vector for new psi for one band, size dim
    psi::Psi<std::complex<FPTYPE>, Device>* phi_m = nullptr;
    // psi::Psi<std::complex<FPTYPE>, Device>* phi_m = nullptr;
    /// temp vector for S|psi> for one band, size dim
    std::vector<std::complex<FPTYPE>> sphi;
    /// temp vector for H|psi> for one band, size dim
    std::vector<std::complex<FPTYPE>> hphi;

    /// temp vector for , size dim
    psi::Psi<std::complex<FPTYPE>, Device>* cg = nullptr;
    // psi::Psi<std::complex<FPTYPE>, Device>* cg = nullptr;
    /// temp vector for , size dim
    std::vector<std::complex<FPTYPE>> scg;
    /// temp vector for store psi in sorting with eigenvalues, size dim
    std::vector<std::complex<FPTYPE>> pphi;

    /// temp vector for , size dim
    std::vector<std::complex<FPTYPE>> gradient;
    /// temp vector for , size dim
    std::vector<std::complex<FPTYPE>> g0;
    /// temp vector for matrix eigenvector * vector S|psi> , size m_band
    std::vector<std::complex<FPTYPE>> lagrange;

    /// device type of psi
    psi::AbacusDevice_t device = {};
    Device * ctx = {};

    void calculate_gradient();

    void orthogonal_gradient(hamilt::Hamilt *phm_in, const psi::Psi<std::complex<FPTYPE>> &eigenfunction, const int m);

    void calculate_gamma_cg(const int iter, FPTYPE &gg_last, const FPTYPE &cg0, const FPTYPE &theta);

    bool update_psi(FPTYPE &cg_norm, FPTYPE &theta, FPTYPE &eigenvalue);

    void schmit_orth(const int &m, const psi::Psi<std::complex<FPTYPE>> &psi);

    // used in diag() for template replace Hamilt with Hamilt_PW
    void diag_mock(hamilt::Hamilt *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &phi, FPTYPE *eigenvalue_in);
};

} // namespace hsolver
#endif

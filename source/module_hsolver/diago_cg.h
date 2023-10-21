#ifndef DIAGCG_H
#define DIAGCG_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

#include "module_psi/kernels/types.h"
#include "module_psi/kernels/device.h"
#include "module_psi/kernels/memory_op.h"

#include "module_hsolver/kernels/math_kernel_op.h"
#include <module_base/macros.h>

namespace hsolver {

template<typename T = std::complex<double>, typename Device = psi::DEVICE_CPU>
class DiagoCG : public DiagH<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  public:
    // Constructor need:
    // 1. temporary mock of Hamiltonian "Hamilt_PW"
    // 2. precondition pointer should point to place of precondition array.
    DiagoCG(const Real *precondition_in);
    ~DiagoCG();

    // virtual void init(){};
    // refactor hpsi_info
    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<T, Device> &psi, Real *eigenvalue_in) override;

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
    const Real *precondition = nullptr;
    /// eigenvalue results
    Real *eigenvalue = nullptr;
    Real *d_precondition = nullptr;

    /// temp vector for new psi for one band, size dim
    psi::Psi<T, Device>* phi_m = nullptr;
    // psi::Psi<T, Device>* phi_m = nullptr;
    /// temp vector for S|psi> for one band, size dim
    T* sphi = nullptr;
    /// temp vector for H|psi> for one band, size dim
    T* hphi = nullptr;

    /// temp vector for , size dim
    psi::Psi<T, Device>* cg = nullptr;
    // psi::Psi<T, Device>* cg = nullptr;
    /// temp vector for , size dim
    T* scg = nullptr;
    /// temp vector for store psi in sorting with eigenvalues, size dim
    T* pphi = nullptr;

    /// temp vector for , size dim
    T* gradient = nullptr;
    /// temp vector for , size dim
    T* g0 = nullptr;
    /// temp vector for matrix eigenvector * vector S|psi> , size m_band
    T* lagrange = nullptr;

    /// device type of psi
    psi::AbacusDevice_t device = {};
    Device * ctx = {};
    psi::DEVICE_CPU *cpu_ctx = {};

    void calculate_gradient();

    void orthogonal_gradient(hamilt::Hamilt<T, Device> *phm_in, const psi::Psi<T, Device> &eigenfunction, const int m);

    void calculate_gamma_cg(const int iter, Real &gg_last, const Real &cg0, const Real &theta);

    bool update_psi(Real &cg_norm, Real &theta, Real &eigenvalue);

    void schmit_orth(const int &m, const psi::Psi<T, Device> &psi);

    // used in diag() for template replace Hamilt with Hamilt_PW
    void diag_mock(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<T, Device> &phi, Real *eigenvalue_in);

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
    using dot_real_op = hsolver::dot_real_op<T, Device>;

    using setmem_complex_op = psi::memory::set_memory_op<T, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<T, Device>;
    using resmem_complex_op = psi::memory::resize_memory_op<T, Device>;
    using syncmem_complex_op = psi::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<T, psi::DEVICE_CPU, Device>;
    using castmem_complex_d2h_op = psi::memory::cast_memory_op<T, T, psi::DEVICE_CPU, Device>;

    using resmem_var_op = psi::memory::resize_memory_op<Real, Device>;
    using setmem_var_h_op = psi::memory::set_memory_op<Real, psi::DEVICE_CPU>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<Real, Device, psi::DEVICE_CPU>;

    const T * one = nullptr, * zero = nullptr, * neg_one = nullptr;
};

} // namespace hsolver
#endif
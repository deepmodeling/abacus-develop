#ifndef DIAGO_NEW_DAV_H
#define DIAGO_NEW_DAV_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_base/macros.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/kernels/device.h"

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = psi::DEVICE_CPU>
class Diago_NewDav : public DiagH<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    Diago_NewDav(const Real* precondition_in);
    ~Diago_NewDav();

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt<T, Device>* phm_in,
              psi::Psi<T, Device>& phi,
              Real* eigenvalue_in,
              std::vector<bool>& is_occupied);

    static int PW_DIAG_NDIM;

  private:
    bool is_subspace = false;

    int test_david = 0;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    /// row size for input psi matrix
    int n_band = 0;
    /// col size for input psi matrix
    int dmx = 0;
    /// non-zero col size for inputted psi matrix
    int dim = 0;
    // maximum dimension of the reduced basis set
    int nbase_x = 0;
    /// precondition for cg diag
    const Real* precondition = nullptr;
    Real* d_precondition = nullptr;

    /// eigenvalue results
    Real* eigenvalue = nullptr;

    T* hphi = nullptr; // the product of H and psi in the reduced basis set

    T* hcc = nullptr; // Hamiltonian on the reduced basis

    T* scc = nullptr; // Overlap on the reduced basis

    T* vcc = nullptr; // Eigenvectors of hc

    T* lagrange_matrix = nullptr;

    /// device type of psi
    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    psi::AbacusDevice_t device = {};

    void cal_grad(hamilt::Hamilt<T, Device>* phm_in,
                  const int& dim,
                  const int& nbase,
                  const int& notconv,
                  psi::Psi<T, Device>& basis,
                  T* hphi,
                  T* vcc,
                  const int* unconv,
                  Real* eigenvalue);

    void cal_elem(const int& dim,
                  int& nbase,
                  const int& notconv,
                  const psi::Psi<T, Device>& basis,
                  const T* hphi,
                  T* hcc,
                  T* scc,
                  bool init);

    void refresh(const int& dim,
                 const int& nband,
                 int& nbase,
                 const Real* eigenvalue,
                 const psi::Psi<T, Device>& psi,
                 psi::Psi<T, Device>& basis,
                 T* hphi,
                 T* hcc,
                 T* scc,
                 T* vcc);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     T* hcc,
                     T* scc,
                     const int& nbase_x,
                     Real* eigenvalue,
                     T* vcc,
                     bool init);

    void diag_mock(hamilt::Hamilt<T, Device>* phm_in,
                   psi::Psi<T, Device>& psi,
                   Real* eigenvalue_in,
                   std::vector<bool>& is_occupied);

    using resmem_complex_op = psi::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<T, Device>;
    using setmem_complex_op = psi::memory::set_memory_op<T, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<Real, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<Real, Device>;
    using setmem_var_op = psi::memory::set_memory_op<Real, Device>;

    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<Real, Device, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<Real, psi::DEVICE_CPU, Device>;
    using syncmem_complex_op = psi::memory::synchronize_memory_op<T, Device, Device>;
    using castmem_complex_op = psi::memory::cast_memory_op<std::complex<double>, T, Device, Device>;
    using syncmem_h2d_op = psi::memory::synchronize_memory_op<T, Device, psi::DEVICE_CPU>;
    using syncmem_d2h_op = psi::memory::synchronize_memory_op<T, psi::DEVICE_CPU, Device>;

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    consts<T> cs;
    const T *one = nullptr, *zero = nullptr, *neg_one = nullptr;
};

template <typename Real, typename Device>
int Diago_NewDav<Real, Device>::PW_DIAG_NDIM = 4;

} // namespace hsolver

#endif
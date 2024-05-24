#ifndef DIAGO_NEW_DAV_H
#define DIAGO_NEW_DAV_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_base/macros.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_base/module_device/device.h"

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class Diago_DavSubspace : public DiagH<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    Diago_DavSubspace(const Real* precondition_in);
    ~Diago_DavSubspace();

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt<T, Device>* phm_in,
              psi::Psi<T, Device>& phi,
              Real* eigenvalue_in,
              std::vector<bool>& is_occupied);

    static int PW_DIAG_NDIM;

    static double dav_large_thr;

  private:
    bool is_subspace = false;

    int test_david = 0;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    /// row size for input psi matrix
    int n_band = 0;
    /// non-zero col size for inputted psi matrix
    int dim = 0;
    // maximum dimension of the reduced basis set
    int nbase_x = 0;
    /// precondition for cg diag
    const Real* precondition = nullptr;
    Real* d_precondition = nullptr;

    /// eigenvalue results
    Real* eigenvalue_in_dav = nullptr;

    T* hphi = nullptr; // the product of H and psi in the reduced basis set

    T* hcc = nullptr; // Hamiltonian on the reduced basis

    T* scc = nullptr; // Overlap on the reduced basis

    T* vcc = nullptr; // Eigenvectors on the reduced basis

    /// device type of psi
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

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
                  T* scc);

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
                     bool init,
                     bool is_subspace);

    void diag_once(hamilt::Hamilt<T, Device>* phm_in,
                   psi::Psi<T, Device>& psi,
                   Real* eigenvalue_in,
                   std::vector<bool>& is_occupied);

    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;

    using resmem_real_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_real_op = base_device::memory::delete_memory_op<Real, Device>;
    using setmem_real_op = base_device::memory::set_memory_op<Real, Device>;

    using resmem_real_h_op = base_device::memory::resize_memory_op<Real, base_device::DEVICE_CPU>;
    using delmem_real_h_op = base_device::memory::delete_memory_op<Real, base_device::DEVICE_CPU>;
    using setmem_real_h_op = base_device::memory::set_memory_op<Real, base_device::DEVICE_CPU>;

    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using castmem_complex_op = base_device::memory::cast_memory_op<std::complex<double>, T, Device, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    consts<T> cs;
    const T *one = nullptr, *zero = nullptr, *neg_one = nullptr;
};

template <typename Real, typename Device>
int Diago_DavSubspace<Real, Device>::PW_DIAG_NDIM = 4;

template <typename Real, typename Device>
double Diago_DavSubspace<Real, Device>::dav_large_thr = 1e-5;

} // namespace hsolver

#endif
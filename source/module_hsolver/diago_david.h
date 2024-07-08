#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "diagh.h"
#include "module_hsolver/diag_comm_info.h"
#include "module_hsolver/diag_const_nums.h"

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoDavid : public DiagH<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  
  public:

    DiagoDavid(const Real* precondition_in, 
               const int david_ndim_in,
               const bool use_paw_in,
               const diag_comm_info& diag_comm_in);

    virtual ~DiagoDavid() override;

    int diag(hamilt::Hamilt<T, Device>* phm_in,
                      psi::Psi<T, Device>& psi,
                      Real* eigenvalue_in,
                      const Real david_diag_thr,
                      const int david_maxiter,
                      const int ntry_max = 5,
                      const int notconv_max = 0);

  private:
    int david_ndim = 4;
    bool use_paw = false;
    int test_david = 0;

    diag_comm_info diag_comm;

    /// row size for input psi matrix
    int n_band = 0;
    /// non-zero col size for inputted psi matrix
    int dim = 0;
    int dmx = 0;
    // maximum dimension of the reduced basis set
    int nbase_x = 0;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    /// precondition for cg diag
    const Real* precondition = nullptr;
    Real* d_precondition = nullptr;

    /// eigenvalue results
    Real* eigenvalue = nullptr;

    T* hphi = nullptr; // the product of H and psi in the reduced basis set

    T* sphi = nullptr; // the Product of S and psi in the reduced basis set

    T* hcc = nullptr; // Hamiltonian on the reduced basis

    T* scc = nullptr; // Overlap on the reduced basis

    T* vcc = nullptr; // Eigenvectors of hc

    T* lagrange_matrix = nullptr;

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
                  T* sphi,
                  const T* vcc,
                  const int* unconv,
                  const Real* eigenvalue);

    void cal_elem(const int& dim,
                  int& nbase,
                  const int& notconv,
                  const psi::Psi<T, Device>& basis,
                  const T* hphi,
                  const T* sphi,
                  T* hcc,
                  T* scc);

    void refresh(const int& dim,
                 const int& nband,
                 int& nbase,
                 const Real* eigenvalue,
                 const psi::Psi<T, Device>& psi,
                 psi::Psi<T, Device>& basis,
                 T* hphi,
                 T* sphi,
                 T* hcc,
                 T* scc,
                 T* vcc);

    void SchmitOrth(const int& dim,
                    const int nband,
                    const int m,
                    psi::Psi<T, Device>& basis,
                    const T* sphi,
                    T* lagrange_m,
                    const int mm_size,
                    const int mv_size);

    void planSchmitOrth(const int nband, int* pre_matrix_mm_m, int* pre_matrix_mv_m);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     const T* hcc,
                     const T* scc,
                     const int& nbase_x,
                     Real* eigenvalue,
                     T* vcc);

    int diag_mock(hamilt::Hamilt<T, Device>* phm_in,
                   psi::Psi<T, Device>& psi,
                   Real* eigenvalue_in,
                   const Real david_diag_thr,
                   const int david_maxiter);

    bool check_block_conv(const int &ntry, const int &notconv, const int &ntry_max, const int &notconv_max);

    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_var_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<Real, Device>;
    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;

    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using castmem_complex_op = base_device::memory::cast_memory_op<std::complex<double>, T, Device, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    const_nums<T> cs;
    const T* one = nullptr, * zero = nullptr, * neg_one = nullptr;
};
} // namespace hsolver

#endif

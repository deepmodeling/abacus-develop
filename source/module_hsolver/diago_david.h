//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_psi/kernels/device.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

namespace hsolver
{

template <typename FPTYPE = double, typename Device = psi::DEVICE_CPU> class DiagoDavid : public DiagH<FPTYPE, Device>
{
  public:
    DiagoDavid(const FPTYPE* precondition_in);
    ~DiagoDavid();

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt<FPTYPE, Device>* phm_in,
              psi::Psi<std::complex<FPTYPE>, Device>& phi,
              FPTYPE* eigenvalue_in);

    static int PW_DIAG_NDIM;

  private:
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
    const FPTYPE* precondition = nullptr;
    FPTYPE* d_precondition = nullptr;

    /// eigenvalue results
    FPTYPE* eigenvalue = nullptr;

    std::complex<FPTYPE>* hphi = nullptr; // the product of H and psi in the reduced basis set

    std::complex<FPTYPE>* sphi = nullptr; // the Product of S and psi in the reduced basis set

    std::complex<FPTYPE>* hcc = nullptr; // Hamiltonian on the reduced basis

    std::complex<FPTYPE>* scc = nullptr; // Overlap on the reduced basis

    std::complex<FPTYPE>* vcc = nullptr; // Eigenvectors of hc

    std::complex<FPTYPE>* lagrange_matrix = nullptr;

    /// device type of psi
    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    psi::AbacusDevice_t device = {};

    void cal_grad(hamilt::Hamilt<FPTYPE, Device>* phm_in,
                  const int& dim,
                  const int& nbase,
                  const int& notconv,
                  psi::Psi<std::complex<FPTYPE>, Device>& basis,
                  std::complex<FPTYPE>* hphi,
                  std::complex<FPTYPE>* sphi,
                  const std::complex<FPTYPE>* vcc,
                  const int* unconv,
                  const FPTYPE* eigenvalue);

    void cal_elem(const int& dim,
                  int& nbase,
                  const int& notconv,
                  const psi::Psi<std::complex<FPTYPE>, Device>& basis,
                  const std::complex<FPTYPE>* hphi,
                  const std::complex<FPTYPE>* sphi,
                  std::complex<FPTYPE>* hcc,
                  std::complex<FPTYPE>* scc);

    void refresh(const int& dim,
                 const int& nband,
                 int& nbase,
                 const FPTYPE* eigenvalue,
                 const psi::Psi<std::complex<FPTYPE>, Device>& psi,
                 psi::Psi<std::complex<FPTYPE>, Device>& basis,
                 std::complex<FPTYPE>* hphi,
                 std::complex<FPTYPE>* sphi,
                 std::complex<FPTYPE>* hcc,
                 std::complex<FPTYPE>* scc,
                 std::complex<FPTYPE>* vcc);

    void SchmitOrth(const int& dim,
                    const int nband,
                    const int m,
                    psi::Psi<std::complex<FPTYPE>, Device>& basis,
                    const std::complex<FPTYPE>* sphi,
                    std::complex<FPTYPE>* lagrange_m,
                    const int mm_size,
                    const int mv_size);

    void planSchmitOrth(const int nband, int* pre_matrix_mm_m, int* pre_matrix_mv_m);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     const std::complex<FPTYPE>* hcc,
                     const std::complex<FPTYPE>* scc,
                     const int& nbase_x,
                     FPTYPE* eigenvalue,
                     std::complex<FPTYPE>* vcc);

    void diag_mock(hamilt::Hamilt<FPTYPE, Device>* phm_in,
                   psi::Psi<std::complex<FPTYPE>, Device>& psi,
                   FPTYPE* eigenvalue_in);

    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using setmem_complex_op = psi::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using setmem_var_op = psi::memory::set_memory_op<FPTYPE, Device>;

    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, Device>;
    using syncmem_complex_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, Device>;
    using castmem_complex_op = psi::memory::cast_memory_op<std::complex<double>, std::complex<FPTYPE>, Device, Device>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;

    using hpsi_info = typename hamilt::Operator<std::complex<FPTYPE>, Device>::hpsi_info;

    const std::complex<FPTYPE> * one = nullptr, * zero = nullptr, * neg_one = nullptr;
};

} // namespace hsolver

#endif

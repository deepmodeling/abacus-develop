#ifndef DIAGO_NEW_DAV_H
#define DIAGO_NEW_DAV_H

#include "module_base/macros.h"   // GetRealType
#include "module_base/module_device/device.h"   // base_device
#include "module_base/module_device/memory_op.h"// base_device::memory"

#include "module_hsolver/diag_comm_info.h"
#include "module_hsolver/diag_const_nums.h"

#include <vector>
#include <functional>

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class Diago_DavSubspace
{
  private:
    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    Diago_DavSubspace(const std::vector<Real>& precondition_in,
                      const int& nband_in,
                      const int& nbasis_in,
                      const int& david_ndim_in,
                      const double& diag_thr_in,
                      const int& diag_nmax_in,
                      const bool& need_subspace_in,
                      const diag_comm_info& diag_comm_in);

    ~Diago_DavSubspace();

    // See diago_david.h for information on the HPsiFunc function type
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    int diag(const HPsiFunc& hpsi_func,
             T* psi_in,
             const int psi_in_dmax,
             Real* eigenvalue_in,
             const double* ethr_band,
             const bool& scf_type);

  private:
    /// for MPI communication
    const diag_comm_info diag_comm;

    /// the threshold for this electronic iteration
    const double diag_thr;

    /// maximal iteration number
    const int iter_nmax;

    /// is diagH_subspace needed?
    const bool is_subspace;

    /// the first dimension of the matrix to be diagonalized
    const int n_band = 0;

    /// the second dimension of the matrix to be diagonalized
    const int dim = 0;

    /// the maximum dimension of the reduced basis set
    const int nbase_x = 0;

    /// precondition for diag
    const std::vector<Real>& precondition;
    Real* d_precondition = nullptr;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    T* psi_in_iter = nullptr;

    /// the product of H and psi in the reduced basis set
    T* hphi = nullptr;

    /// Hamiltonian on the reduced basis
    T* hcc = nullptr;

    /// Overlap on the reduced basis
    T* scc = nullptr;

    /// Eigenvectors on the reduced basis
    T* vcc = nullptr;

    /// device type of psi
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

    void cal_grad(const HPsiFunc& hpsi_func,
                  const int& dim,
                  const int& nbase,
                  const int& notconv,
                  T* psi_iter,
                  T* hphi,
                  T* vcc,
                  const int* unconv,
                  std::vector<Real>* eigenvalue_iter);

    void cal_elem(const int& dim, int& nbase, const int& notconv, const T* psi_iter, const T* hphi, T* hcc, T* scc);

    void refresh(const int& dim,
                 const int& nband,
                 int& nbase,
                 const Real* eigenvalue,
                 T* psi_iter,
                 T* hphi,
                 T* hcc,
                 T* scc,
                 T* vcc);

    // void diagH_LAPACK(const int nstart,
    //                   const int nbands,
    //                   const T* hcc,
    //                   const T* sc,
    //                   const int ldh, // nstart
    //                   Real* e,
    //                   T* vcc);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     T* hcc,
                     T* scc,
                     const int& nbase_x,
                     std::vector<Real>* eigenvalue_iter,
                     T* vcc);

    int diag_once(const HPsiFunc& hpsi_func,
                  T* psi_in,
                  const int psi_in_dmax,
                  Real* eigenvalue_in,
                  const double* ethr_band);

    bool test_exit_cond(const int& ntry, const int& notconv, const bool& scf);

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

    const T *one = nullptr, *zero = nullptr, *neg_one = nullptr;
    const T one_ = static_cast<T>(1.0), zero_ = static_cast<T>(0.0), neg_one_ = static_cast<T>(-1.0);
};

} // namespace hsolver

#endif
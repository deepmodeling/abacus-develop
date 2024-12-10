#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "module_base/macros.h"   // GetRealType
#include "module_base/module_device/device.h"   // base_device
#include "module_base/module_device/memory_op.h"// base_device::memory

#include "module_hsolver/diag_comm_info.h"

#include <vector>
#include <functional>

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoDavid
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  
  public:

    DiagoDavid(const Real* precondition_in,
               const int nband_in,
               const int dim_in,
               const int david_ndim_in,
               const bool use_paw_in,
               const diag_comm_info& diag_comm_in);

     ~DiagoDavid();


    // declare type of matrix-blockvector functions.
    // the function type is defined as a std::function object.
    /**
     * @brief A function type representing the HX function.
     *
     * This function type is used to define a matrix-blockvector operator H.
     * For eigenvalue problem HX = λX or generalized eigenvalue problem HX = λSX,
     * this function computes the product of the Hamiltonian matrix H and a blockvector X.
     * 
     * Called as follows:
     * hpsi(X, HX, ld, nvec) where X and HX are (ld, nvec)-shaped blockvectors.
     * Result HX = H * X is stored in HX.
     *
     * @param[out] X      Head address of input blockvector of type `T*`.
     * @param[in]  HX     Head address of output blockvector of type `T*`.
     * @param[in]  ld     Leading dimension of blockvector.
     * @param[in]  nvec   Number of vectors in a block.
     * 
     * @warning X and HX are the exact address to read input X and store output H*X,
     * @warning both of size ld * nvec.
     */
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /**
     * @brief A function type representing the SX function.
     * 
     * nrow is leading dimension of spsi, npw is leading dimension of psi, nbands is number of vecs
     *
     * This function type is used to define a matrix-blockvector operator S.
     * For generalized eigenvalue problem HX = λSX,
     * this function computes the product of the overlap matrix S and a blockvector X.
     *
     * @param[in]   X       Pointer to the input blockvector.
     * @param[out] SX       Pointer to the output blockvector.
     * @param[in] ld_psi    Leading dimension of psi and spsi. Dimension of X&SX: ld * nvec.
     * @param[in] nvec      Number of vectors.
     */
    using SPsiFunc = std::function<void(T*, T*, const int, const int)>;

    int diag(
      const HPsiFunc& hpsi_func,  // function void hpsi(T*, T*, const int, const int) 
      const SPsiFunc& spsi_func,  // function void spsi(T*, T*, const int, const int, const int) 
      const int ld_psi,           // Leading dimension of the psi input
      T *psi_in,                  // Pointer to eigenvectors
      Real* eigenvalue_in,        // Pointer to store the resulting eigenvalues
      const std::vector<double>& ethr_band, // Convergence threshold for the Davidson iteration
      const int david_maxiter,    // Maximum allowed iterations for the Davidson method
      const int ntry_max = 5,     // Maximum number of diagonalization attempts (5 by default)
      const int notconv_max = 0); // Maximum number of allowed non-converged eigenvectors

  private:
    bool use_paw = false;
    int test_david = 0;

    diag_comm_info diag_comm;

    /// number of required eigenpairs
    const int nband;
    /// dimension of the input matrix to be diagonalized
    const int dim;
    /// maximum dimension of the reduced basis set
    const int nbase_x;
    /// dimension of the subspace allowed in Davidson
    const int david_ndim = 4;
    /// number of unconverged eigenvalues
    int notconv = 0;

    /// precondition for diag, diagonal approximation of matrix A(i.e. Hamilt)
    const Real* precondition = nullptr;
    Real* d_precondition = nullptr;

    /// eigenvalue results
    Real* eigenvalue = nullptr;

    T *basis = nullptr;  /// pointer to basis set(dim, nbase_x), leading dimension = dim

    T* hpsi = nullptr;    /// the product of H and psi in the reduced basis set

    T* spsi = nullptr;    /// the Product of S and psi in the reduced basis set

    T* hcc = nullptr;     /// Hamiltonian on the reduced basis

    T* scc = nullptr;     /// overlap on the reduced basis

    T* vcc = nullptr;     /// eigenvectors of hc

    T* lagrange_matrix = nullptr;

    /// device type of psi
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

    int diag_once(const HPsiFunc& hpsi_func,
                  const SPsiFunc& spsi_func,
                  const int dim,
                  const int nband,
                  const int ld_psi,
                  T *psi_in,
                  Real* eigenvalue_in,
                  const std::vector<double>& ethr_band,
                  const int david_maxiter);

    void cal_grad(const HPsiFunc& hpsi_func,
                  const SPsiFunc& spsi_func,
                  const int& dim,
                  const int& nbase,
                  const int nbase_x,
                  const int& notconv,
                  T* hpsi,
                  T* spsi,
                  const T* vcc,
                  const int* unconv,
                  const Real* eigenvalue);

    void cal_elem(const int& dim,
                  int& nbase,
                  const int nbase_x,
                  const int& notconv,
                  const T* hpsi,
                  const T* spsi,
                  T* hcc,
                  T* scc);

    void refresh(const int& dim,
                 const int& nband,
                 int& nbase,
                 const int nbase_x,
                 const Real* eigenvalue,
                 const T *psi_in,
                 const int ld_psi,
                 T* hpsi,
                 T* spsi,
                 T* hcc,
                 T* scc,
                 T* vcc);

    void SchmidtOrth(const int& dim,
                     const int nband,
                     const int m,
                     const T* spsi,
                     T* lagrange_m,
                     const int mm_size,
                     const int mv_size);

    void planSchmidtOrth(const int nband, std::vector<int>& pre_matrix_mm_m, std::vector<int>& pre_matrix_mv_m);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     const T* hcc,
                     const T* scc,
                     const int& nbase_x,
                     Real* eigenvalue,
                     T* vcc);

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

    // using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info; // Dependence of hpsi removed

    const T *one = nullptr, *zero = nullptr, *neg_one = nullptr;
    const T one_ = static_cast<T>(1.0), zero_ = static_cast<T>(0.0), neg_one_ = static_cast<T>(-1.0);
};
} // namespace hsolver

#endif

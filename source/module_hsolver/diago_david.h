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
/**
 * @class DiagoDavid
 * @brief A class that implements the block-Davidson algorithm for solving generalized eigenvalue problems.
 *
 * The DiagoDavid class provides methods for performing iterative diagonalization using the Davidson algorithm.
 * It supports both real and complex data types and can be executed on different devices (CPU or GPU).
 *
 * @tparam T The data type of the matrices and arrays (e.g., float, double, std::complex<float>, std::complex<double>).
 * @tparam Device The device type (e.g., base_device::DEVICE_CPU or DEVICE_GPU).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoDavid
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  
  public:

    /**
     * @brief Constructor for the DiagoDavid class.
     * 
     * @param[in] precondition_in Pointer to the preconditioning matrix.
     * @param[in] nband_in Number of eigenpairs required(i.e. bands).
     * @param[in] dim_in Dimension of the matrix.
     * @param[in] david_ndim_in Dimension of the reduced basis set of Davidson.
     *                      `david_ndim_in` * `nband_in` is the maximum allowed size of
     *                      the reduced basis set before \b restart of Davidson.
     * @param[in] use_paw_in Flag indicating whether to use PAW.
     * @param[in] diag_comm_in Communication information for diagonalization.
     * 
     * @tparam T The data type of the matrices and arrays.
     * @tparam Device The device type (base_device::DEVICE_CPU or DEVICE_GPU).
     * 
     * @note Auxiliary memory is allocated in the constructor and deallocated in the destructor.
     */
    DiagoDavid(const Real* precondition_in,
               const int nband_in,
               const int dim_in,
               const int david_ndim_in,
               const bool use_paw_in,
               const diag_comm_info& diag_comm_in);

    /**
     * @brief Destructor for the DiagoDavid class.
     * 
     * This destructor releases the dynamically allocated memory used by the class members.
     * It deletes the basis, hpsi, spsi, hcc, scc, vcc, lagrange_matrix, and eigenvalue arrays.
     * 
     */
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

    /**
     * @brief Performs iterative diagonalization using the David algorithm.
     * 
     * @warning Please see docs of `HPsiFunc` for more information about the hpsi mat-vec interface.
     * 
     * @tparam T The type of the elements in the matrix.
     * @tparam Device The device type (CPU or GPU).
     * @param hpsi_func The function object that computes the matrix-blockvector product H * psi.
     * @param spsi_func The function object that computes the matrix-blockvector product overlap S * psi.
     * @param ld_psi The leading dimension of the psi_in array.
     * @param psi_in The input wavefunction.
     * @param eigenvalue_in The array to store the eigenvalues.
     * @param david_diag_thr The convergence threshold for the diagonalization.
     * @param david_maxiter The maximum number of iterations for the diagonalization.
     * @param ntry_max The maximum number of attempts for the diagonalization restart.
     * @param notconv_max The maximum number of bands unconverged allowed.
     * @return The total number of iterations performed during the diagonalization.
     * 
     * @note ntry_max is an empirical parameter that should be specified in external routine, default 5
     *       notconv_max is determined by the accuracy required for the calculation, default 0
     */
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

    /**
     * Calculates the preconditioned gradient of the eigenvectors in Davidson method.
     *
     * @param hpsi_func The function to calculate the matrix-blockvector product H * psi.
     * @param spsi_func The function to calculate the matrix-blockvector product overlap S * psi.
     * @param dim The dimension of the blockvector.
     * @param nbase The current dimension of the reduced basis.
     * @param nbase_x The maximum dimension of the reduced basis set.
     * @param notconv The number of unconverged eigenpairs.
     * @param hpsi The output array for the Hamiltonian H times blockvector psi.
     * @param spsi The output array for the overlap matrix S times blockvector psi.
     * @param vcc The input array for the eigenvector coefficients.
     * @param unconv The array of indices for the unconverged eigenpairs.
     * @param eigenvalue The array of eigenvalues.
     */
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

    /**
     * Calculates the elements of the diagonalization matrix for the DiagoDavid class.
     * 
     * @param dim The dimension of the problem.
     * @param nbase The current dimension of the reduced basis.
     * @param nbase_x The maximum dimension of the reduced basis set.
     * @param notconv The number of newly added basis vectors.
     * @param hpsi The output array for the Hamiltonian H times blockvector psi.
     * @param spsi The output array for the overlap matrix S times blockvector psi.
     * @param hcc Pointer to the array where the calculated Hamiltonian matrix elements will be stored.
     * @param scc Pointer to the array where the calculated overlap matrix elements will be stored.
     */
    void cal_elem(const int& dim,
                  int& nbase,
                  const int nbase_x,
                  const int& notconv,
                  const T* hpsi,
                  const T* spsi,
                  T* hcc,
                  T* scc);

    /**
     * Refreshes the diagonalization solver by updating the basis and the reduced Hamiltonian.
     * 
     * @param dim The dimension of the problem.
     * @param nband The number of bands.
     * @param nbase The number of basis states.
     * @param nbase_x The maximum dimension of the reduced basis set.
     * @param eigenvalue_in Pointer to the array of eigenvalues.
     * @param psi_in Pointer to the array of wavefunctions.
     * @param ld_psi The leading dimension of the wavefunction array.
     * @param hpsi Pointer to the output array for the updated basis set.
     * @param spsi Pointer to the output array for the updated basis set (nband-th column).
     * @param hcc Pointer to the output array for the updated reduced Hamiltonian.
     * @param scc Pointer to the output array for the updated overlap matrix.
     * @param vcc Pointer to the output array for the updated eigenvector matrix.
     * 
     */
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

    /**
     * SchmidtOrth function performs orthogonalization of the starting eigenfunction to those already calculated.
     * It takes the dimension of the basis, number of bands, index of the current band, starting eigenfunction psi_m,
     * lagrange_m array, mm_size, and mv_size as input parameters.
     *
     * @param dim The dimension of the basis.
     * @param nband The number of bands.
     * @param m The index of the current band.
     * @param spsi Pointer to the starting eigenfunction psi_m.
     * @param lagrange_m Pointer to the lagrange_m array.
     * @param mm_size The size of the square matrix for future lagranges.
     * @param mv_size The size of the lagrange_m array.
     */
    void SchmidtOrth(const int& dim,
                     const int nband,
                     const int m,
                     const T* spsi,
                     T* lagrange_m,
                     const int mm_size,
                     const int mv_size);

    /**
     * @brief Plans the Schmidt orthogonalization for a given number of bands.
     * 
     * @tparam T The type of the elements in the vectors.
     * @tparam Device The device on which the computation will be performed.
     * @param nband The number of bands.
     * @param pre_matrix_mm_m The vector to store the matrix sizes.
     * @param pre_matrix_mv_m The vector to store the number of matrix-vector multiplications.
     */
    void planSchmidtOrth(const int nband, std::vector<int>& pre_matrix_mm_m, std::vector<int>& pre_matrix_mv_m);

    void diag_zhegvx(const int& nbase,
                     const int& nband,
                     const T* hcc,
                     const T* scc,
                     const int& nbase_x,
                     Real* eigenvalue,
                     T* vcc);

    /**
     * @brief Check the convergence of block eigenvectors in the Davidson iteration.
     *
     * This function determines whether the block eigenvectors have reached convergence
     * during the iterative diagonalization process. Convergence is judged based on
     * the number of eigenvectors that have not converged and the maximum allowed
     * number of such eigenvectors.
     *
     * @tparam T The data type for the eigenvalues and eigenvectors (e.g., float, double).
     * @tparam Device The device type (e.g., base_device::DEVICE_CPU).
     * @param ntry The current number of tries for diagonalization.
     * @param notconv The current number of eigenvectors that have not converged.
     * @param ntry_max The maximum allowed number of tries for diagonalization.
     * @param notconv_max The maximum allowed number of eigenvectors that can fail to converge.
     * @return true if the eigenvectors are considered converged or the maximum number
     *         of tries has been reached, false otherwise.
     *
     * @note Exits the diagonalization loop if either the convergence criteria
     *       are met or the maximum number of tries is exceeded.
     */
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

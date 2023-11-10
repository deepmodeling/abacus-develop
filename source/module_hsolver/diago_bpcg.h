#ifndef DIAGO_BPCG_H_
#define DIAGO_BPCG_H_

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

#include "module_psi/kernels/types.h"
#include "module_psi/kernels/device.h"
#include "module_psi/kernels/memory_op.h"

#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include <module_base/macros.h>

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>

namespace hsolver {

/**
 * @class DiagoBPCG
 * @brief A class for diagonalization  using the Blocked-PCG method.
 * @tparam T The floating-point type used for calculations.
 * @tparam Device The device used for calculations (e.g., cpu or gpu).
 */
template<typename T = std::complex<double>, typename Device = psi::DEVICE_CPU>
class DiagoBPCG : public DiagH<T, Device>
{
  private:
    // Note GetTypeReal<T>::type will 
    // return T if T is real type(float, double), 
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;
  // Column major psi in this class
  public:
    /**
     * @brief Constructor for DiagoBPCG class.
     *
     * @param precondition precondition data passed by the "Hamilt_PW" class.
     */
    explicit DiagoBPCG(const Real* precondition);

    /**
     * @brief Destructor for DiagoBPCG class.
     */
    ~DiagoBPCG();

    /**
     * @brief Initialize the class before diagonalization.
     *
     * This function allocates all the related variables, such as hpsi, hsub, before the diag call.
     * It is called by the HsolverPW::initDiagh() function.
     *
     * @param psi_in The input wavefunction psi.
     */
    void init_iter(const psi::Psi<T, Device> &psi_in);

    /**
     * @brief Diagonalize the Hamiltonian using the CG method.
     *
     * This function is an override function for the CG method. It is called by the HsolverPW::solve() function.
     *
     * @param phm_in A pointer to the hamilt::Hamilt object representing the Hamiltonian operator.
     * @param psi The input wavefunction psi matrix with [dim: n_basis x n_band, column major].
     * @param eigenvalue_in Pointer to the eigen array with [dim: n_band, column major].
     */
    void diag(hamilt::Hamilt<T, Device> *phm_in, psi::Psi<T, Device> &psi, Real *eigenvalue_in) override;


  private:
    /// the number of rows of the input psi
    int n_band = 0;
    /// the number of cols of the input psi
    int n_basis = 0;
    /// max iter steps for all-band cg loop
    int nline = 4;
    /// cg convergence thr
    Real all_band_cg_thr = 1E-5;

    ct::DataType r_type  = ct::DataType::DT_INVALID;
    ct::DataType t_type  = ct::DataType::DT_INVALID;
    ct::DeviceType device_type = ct::DeviceType::UnKnown;

    ct::Tensor prec = {}, h_prec = {};

    /// The coefficient for mixing the current and previous step gradients, used in iterative methods.
    ct::Tensor beta = {};
    /// Error state value, if it is smaller than the given threshold, then exit the iteration.
    ct::Tensor err_st = {};
    /// Calculated eigen
    ct::Tensor eigen = {};

    /// Pointer to the input wavefunction.
    /// Note: this pointer does not own memory, instead it ref the psi_in object.
    /// H|psi> matrix.
    ct::Tensor psi = {}, hpsi = {};
    
    ct::Tensor hsub = {};

    /// H|psi> - epsilo * psi, grad of the given problem.
    /// Dim: n_basis * n_band, column major, lda = n_basis_max.
    ct::Tensor grad = {}, hgrad = {}, grad_old = {};

    /// work for some calculations within this class, including rotate_wf call
    ct::Tensor work = {};

    psi::Psi<T, Device>* grad_wrapper;
    /**
     * @brief Update the precondition array.
     *
     * This function updates the precondition array by copying the host precondition
     * to the device in a 'gpu' runtime environment. The address of the precondition
     * array is passed by the constructor function called by hsolver::HSolverPW::initDiagh.
     * The precondition will be updated before the DiagoBPCG<T, Device>::diag call.
     *
     * @note prec[dim: n_band]
     *
     * @param dev Reference to the AbacusDevice_t object, speciy which device used in the calc_prec function.
     * @param prec Pointer to the host precondition array with [dim: n_band, column major]
     * @param h_prec Pointer to the host precondition array with [dim: n_band, column major].
     * @param d_prec Pointer to the device precondition array with [dim: n_band, column major].
     */
    void calc_prec();

    /**
     *
     * @brief Apply the H operator to psi and obtain the hpsi matrix.
     *
     * This function calculates the matrix product of the Hamiltonian operator (H) and the input wavefunction (psi).
     * The resulting matrix is stored in the output array hpsi_out.
     *
     * @note hpsi = H|psi>;
     *
     * psi_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_out[dim: n_basis x n_band, column major, lda = n_basis_max].
     *
     * @param hamilt_in A pointer to the hamilt::Hamilt object representing the Hamiltonian operator.
     * @param psi_in The input wavefunction psi.
     * @param hpsi_out Pointer to the array where the resulting hpsi matrix will be stored.
     */
    void calc_hpsi_with_block(
        hamilt::Hamilt<T, Device>* hamilt_in, 
        const psi::Psi<T, Device>& psi_in,  
        ct::Tensor& hpsi_out);

    /**
     * @brief Diagonalization of the subspace matrix.
     *
     * All the matrix used in this function are stored and used as the column major.
     * psi_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * eigenvalue_out[dim: n_basis_max, column major].
     *
     * @param psi_in Input wavefunction matrix with [dim: n_basis x n_band, column major].
     * @param hpsi_in H|psi> matrix with [dim: n_basis x n_band, column major].
     * @param hsub_out Output Hamiltonian subtracted matrix with [dim: n_band x n_band, column major]
     * @param eigenvalue_out Computed eigen array with [dim: n_band]
     */
    void diag_hsub(
        const ct::Tensor& psi_in, 
        const ct::Tensor& hpsi_in,
        ct::Tensor& hsub_out, 
        ct::Tensor& eigenvalue_out);

    /**
     * @brief Inplace matrix multiplication to obtain the initial guessed wavefunction.
     *
     * hsub_in[dim: n_band x n_band, column major, lda = n_band],
     * workspace_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * psi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     *
     * @param hsub_in Subspace matrix input, dim [n_basis, n_band] with column major.
     * @param workspace_in Workspace matrix, dim [n_basis, n_band] with column major..
     * @param psi_out output wavefunction matrix with dim [n_basis, n_band], column major.
     */
    void rotate_wf(
        const ct::Tensor& hsub_in,
        ct::Tensor& psi_out, 
        ct::Tensor& workspace_in);

    /**
     * @brief Calculate the gradient for all bands used in CG method.
     *
     * prec_in[dim: n_basis_max, column major],
     * err_out[dim: n_band, column major],
     * beta_out[dim: n_band, column major],
     * psi_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_in[dim: n_basis x n_band, column major, lda = n_basis_max],
     * grad_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * grad_old_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     *
     * @param prec_in Input preconditioner.
     * @param err_out Output error state value. If it is smaller than a given threshold, exit the iteration.
     * @param beta_out Output beta coefficient.
     * @param psi_in Input wavefunction matrix.
     * @param hpsi_in Product of psi_in and Hamiltonian.
     * @param grad_out Output gradient matrix.
     * @param grad_old_out Previous gradient matrix.
     * @note The steps involved in optimization are:
     *   1. normalize psi
     *   2. calculate the epsilo
     *   3. calculate the gradient by hpsi - epsilo * psi
     *   4. gradient mix with the previous gradient
     *   5. Do precondition
     */
    void calc_grad_with_block(
        const ct::Tensor& prec_in, 
        ct::Tensor& err_out, 
        ct::Tensor& beta_out,
        ct::Tensor& psi_in, ct::Tensor& hpsi_in,
        ct::Tensor& grad_out, ct::Tensor& grad_old_out);

    /**
     *
     * @brief Apply the Hamiltonian operator to psi and obtain the hpsi matrix.
     *
     * psi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hsub_out[dim: n_band x n_band, column major, lda = n_band],
     * eigenvalue_out[dim: n_basis_max, column major].
     *
     * @param hamilt_in Pointer to the Hamiltonian object.
     * @param psi_in Input wavefunction.
     * @param psi_out Output wavefunction.
     * @param hpsi_out Product of psi_out and Hamiltonian.
     * @param hsub_out Subspace matrix output.
     * @param eigenvalue_out Computed eigen.
     */
    void calc_hsub_with_block(
        hamilt::Hamilt<T, Device>* hamilt_in,
        const psi::Psi<T, Device>& psi_in,
        ct::Tensor& psi_out, ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out, ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out);
    
    /**
     *
     * @brief Apply the Hamiltonian operator to psi and obtain the hpsi matrix.
     *
     * psi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hpsi_out[dim: n_basis x n_band, column major, lda = n_basis_max],
     * hsub_out[dim: n_band x n_band, column major, lda = n_band],
     * eigenvalue_out[dim: n_basis_max, column major].
     *
     * @param hamilt_in Pointer to the Hamiltonian object.
     * @param psi_in Input wavefunction.
     * @param psi_out Output wavefunction.
     * @param hpsi_out Product of psi_out and Hamiltonian.
     * @param hsub_out Subspace matrix output.
     * @param eigenvalue_out Computed eigen.
     */
    void calc_hsub_with_block_exit(
        ct::Tensor& psi_out, 
        ct::Tensor& hpsi_out,
        ct::Tensor& hsub_out, 
        ct::Tensor& workspace_in,
        ct::Tensor& eigenvalue_out);

    /**
     * @brief Orthogonalize column vectors in grad to column vectors in psi.
     *
     * hsub_in and workspace_in are only used to store intermediate variables of the gemm operator.
     *
     * @param psi_in Input wavefunction array, [dim: n_basis x n_band, column major, lda = n_basis_max].
     * @param hsub_in Subspace matrix input, [dim: n_band x n_band, column major, lda = n_band].
     * @param grad_out Input and output gradient array, [dim: n_basis x n_band, column major, lda = n_basis_max]..
     * @note This function is a member of the DiagoBPCG class.
     */
    void orth_projection(
        const ct::Tensor& psi_in,
        ct::Tensor& hsub_in,
        ct::Tensor& grad_out);

    /**
     *
     *@brief Optimize psi as well as the hpsi.
     *
     *@param grad_in Input gradient array, [dim: n_basis x n_band, column major, lda = n_basis_max].
     *@param hgrad_in Product of grad_in and Hamiltonian, [dim: n_basis x n_band, column major, lda = n_basis_max].
     *@param psi_out Input and output wavefunction array, [dim: n_basis x n_band, column major, lda = n_basis_max].
     *@param hpsi_out Product of psi_out and Hamiltonian, [dim: n_basis x n_band, column major, lda = n_basis_max].
     *@note The steps involved in optimization are:
     *  1. Normalize the gradient.
     *  2. Calculate theta.
     *  3. Update psi as well as hpsi.
     */
    void line_minimize(
        ct::Tensor& grad_in,
        ct::Tensor& hgrad_in,
        ct::Tensor& psi_out,
        ct::Tensor& hpsi_out);

    /**
     * @brief Orthogonalize and normalize the column vectors in psi_out using Cholesky decomposition.
     *
     * @param workspace_in Workspace memory, [dim: n_basis x n_band, column major, lda = n_basis_max]..
     * @param psi_out Input and output wavefunction array. [dim: n_basis x n_band, column major, lda = n_basis_max].
     * @param hpsi_out Input and output hpsi array. [dim: n_basis x n_band, column major, lda = n_basis_max].
     * @param hsub_out Input Hamiltonian product array. [dim: n_band x n_band, column major, lda = n_band].
     */
    void orth_cholesky(
        ct::Tensor& workspace_in, 
        ct::Tensor& psi_out, 
        ct::Tensor& hpsi_out, 
        ct::Tensor& hsub_out);

    /**
     * @brief Checks if the error satisfies the given threshold.
     *
     * @param err_in Pointer to the error array.[dim: n_band, column major]
     * @param thr_in The threshold.
     * @return Returns true if all error values are less than or equal to the threshold, false otherwise.
     */
    bool test_error(const ct::Tensor& err_in, Real thr_in);

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_var_op = ct::kernels::set_memory<Real, ct_Device>;
    using resmem_var_op = ct::kernels::resize_memory<Real, ct_Device>;
    using delmem_var_op = ct::kernels::delete_memory<Real, ct_Device>;
    using syncmem_var_h2d_op = ct::kernels::synchronize_memory<Real, ct_Device, ct::DEVICE_CPU>;
    using syncmem_var_d2h_op = ct::kernels::synchronize_memory<Real, ct::DEVICE_CPU, ct_Device>;

    using setmem_complex_op = ct::kernels::set_memory<T, ct_Device>;
    using delmem_complex_op = ct::kernels::delete_memory<T, ct_Device>;
    using resmem_complex_op = ct::kernels::resize_memory<T, ct_Device>;
    using syncmem_complex_op = ct::kernels::synchronize_memory<T, ct_Device, ct_Device>;

    using calc_grad_with_block_op = hsolver::calc_grad_with_block_op<T, Device>;
    using line_minimize_with_block_op = hsolver::line_minimize_with_block_op<T, Device>;

};

} // namespace hsolver
#endif // DIAGO_BPCG_H_
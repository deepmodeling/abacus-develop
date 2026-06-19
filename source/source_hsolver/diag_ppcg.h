#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/para_gemm.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_hsolver/para_linear_transform.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <source_base/macros.h>

namespace hsolver
{

/**
 * @class DiagoPPCG
 * @brief A class for diagonalization using the Projected Preconditioned Conjugate Gradient (PPCG) method.
 *
 * The DiagoPPCG class implements a block LOBPCG-like algorithm for solving generalized eigenvalue problems.
 * It uses an expanded subspace [X, W, P] where X is the current eigenvector approximation,
 * W is the preconditioned residual, and P is the conjugate search direction from previous steps.
 *
 * @tparam T The floating-point type used for calculations.
 * @tparam Device The device used for calculations (e.g., cpu or gpu).
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoPPCG
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    /**
     * @brief Constructor for DiagoPPCG class.
     *
     * @param precondition_in Pointer to the host precondition array.
     */
    explicit DiagoPPCG(const Real* precondition_in);

    /**
     * @brief Destructor for DiagoPPCG class.
     */
    ~DiagoPPCG();

    /**
     * @brief Initialize the class before diagonalization.
     *
     * @param nband The number of bands of all processes.
     * @param nband_l The number of bands of current process.
     * @param nbasis The number of basis functions. Leading dimension of psi.
     * @param ndim The number of valid dimension of psi.
     */
    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    /**
     * @brief Diagonalize the Hamiltonian using the PPCG method.
     *
     * @param hpsi_func A function computing the product of the Hamiltonian matrix H
     * and a wavefunction blockvector X.
     * @param psi_in Pointer to input wavefunction psi matrix with [dim: n_basis x n_band, column major].
     * @param eigenvalue_in Pointer to the eigen array with [dim: n_band].
     * @param ethr_band Convergence threshold for each band.
     */
    void diag(const HPsiFunc& hpsi_func, T* psi_in, Real* eigenvalue_in, const std::vector<double>& ethr_band);

  private:
    /// the number of bands of all processes
    int n_band = 0;
    /// the number of bands of current process
    int n_band_l = 0;
    /// the number of cols of the input psi
    int n_basis = 0;
    /// valid dimension of psi
    int n_dim = 0;
    /// max iter steps for ppcg loop
    int nline = 4;

    /// parallel matrix multiplication
    ModuleBase::PGemmCN<T, Device> pmmcn;
    PLinearTransform<T, Device> plintrans;

    ct::DataType r_type = ct::DataType::DT_INVALID;
    ct::DataType t_type = ct::DataType::DT_INVALID;
    ct::DeviceType device_type = ct::DeviceType::UnKnown;

    ct::Tensor h_prec = {};
    ct::Tensor prec = {};
    ct::Tensor eigen = {};

    /// Number of globally converged (locked) bands
    int nlocked = 0;
    /// Locked eigenvalues on CPU
    std::vector<Real> eigen_locked;
    /// MPI band distribution for error gathering
    std::vector<int> all_n_band_l;
    std::vector<int> band_displs;
    ct::Tensor err_st = {};

    ct::Tensor psi = {}, hpsi = {};
    ct::Tensor w = {}, hw = {};
    ct::Tensor p = {}, hp = {};
    ct::Tensor work = {};

    Device* ctx = {};
    const T *one = nullptr, *zero = nullptr, *neg_one = nullptr;
    const T one_ = static_cast<T>(1.0), zero_ = static_cast<T>(0.0), neg_one_ = static_cast<T>(-1.0);

    /**
     * @brief Update the precondition array from host to device.
     */
    void calc_prec();

    /**
     * @brief Apply the H operator to psi and obtain the hpsi matrix.
     */
    void calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, ct::Tensor& hpsi_out);

    /**
     * @brief Calculate the preconditioned residual (gradient) and error.
     *
     * @param prec_in Input preconditioner.
     * @param err_out Output error for each local band.
     * @param psi_in Input wavefunction.
     * @param hpsi_in H|psi> matrix.
     * @param grad_out Output preconditioned residual.
     */
    void calc_grad(const ct::Tensor& prec_in,
                   ct::Tensor& err_out,
                   ct::Tensor& psi_in,
                   ct::Tensor& hpsi_in,
                   ct::Tensor& grad_out,
                   const int nlocked_in = 0);

    /**
     * @brief Orthogonalize grad to psi using S-inner product (S=I for norm-conserving).
     *
     * @param psi_in Input wavefunction.
     * @param hsub_tmp Workspace for overlap matrix.
     * @param grad_out Input/Output gradient.
     */
    void orth_projection(const ct::Tensor& psi_in, ct::Tensor& hsub_tmp, ct::Tensor& grad_out);

    /**
     * @brief Perform expanded subspace diagonalization and update X, P, HX, HP.
     *
     * @param hpsi_func Hamiltonian application function.
     * @param has_p If true, use 3-block [X, W, P]; otherwise use 2-block [X, W].
     * @param psi_out Input/Output wavefunction.
     * @param hpsi_out Input/Output H|psi>.
     * @param p_out Input/Output search direction.
     * @param hp_out Input/Output H|p>.
     */
    void diag_subspace(const HPsiFunc& hpsi_func,
                       const bool has_p,
                       ct::Tensor& psi_out,
                       ct::Tensor& hpsi_out,
                       ct::Tensor& p_out,
                       ct::Tensor& hp_out,
                       const int nlocked_in = 0);

    /**
     * @brief Orthogonalize and normalize psi using Cholesky decomposition.
     */
    void orth_cholesky(ct::Tensor& workspace_in, ct::Tensor& psi_out, ct::Tensor& hpsi_out, ct::Tensor& hsub_out);

    /**
     * @brief Update locking status: scan errors from current nlocked forward
     *        and lock bands that have converged.
     */
    void update_locking(const ct::Tensor& err_in, const std::vector<double>& ethr_band);

    /**
     * @brief Check if all bands have converged.
     */
    bool test_error(const ct::Tensor& err_in, const std::vector<double>& ethr_band);

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
    using syncmem_complex_h2d_op = ct::kernels::synchronize_memory<T, ct_Device, ct::DEVICE_CPU>;
    using syncmem_complex_d2h_op = ct::kernels::synchronize_memory<T, ct::DEVICE_CPU, ct_Device>;

    using gemm_op = ModuleBase::gemm_op<T, Device>;
};

} // namespace hsolver

#endif // DIAGO_PPCG_H_

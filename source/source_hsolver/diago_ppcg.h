#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_device/types.h"
#include "source_base/parallel_reduce.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <functional>
#include <vector>

namespace hsolver {

/**
 * @class DiagoPPCG
 * @brief Band-by-band projected preconditioned conjugate gradient method
 *        for eigenvalue problems.
 *
 * This implementation solves the generalized eigenvalue problem H*psi = lambda*psi
 * using a band-by-band preconditioned conjugate gradient approach, closely
 * following the algorithm of DiagoCG. Each band is solved sequentially with
 * explicit projection (Schmidt orthogonalization) against previously converged
 * eigenvectors and a 2D Rayleigh-Ritz line minimization in the [psi, cg] subspace.
 *
 * Key algorithmic components:
 * - Preconditioned gradient with Rayleigh-quotient projection
 * - Polak-Ribiere CG direction update with automatic restart
 * - 2D subspace line minimization via atan-based optimal angle
 * - Outer retry loop with subspace diagonalization for improved initial guesses
 *
 * @tparam T The floating-point type (float, double, complex<float>, complex<double>)
 * @tparam Device The device type (CPU or GPU)
 */
template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoPPCG
{
  private:
    using Real = typename GetTypeReal<T>::type;
    using ct_Device = typename ct::PsiToContainer<Device>::type;

  public:
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;
    using SPsiFunc = std::function<void(T*, T*, const int, const int)>;
    using SubspaceFunc = std::function<void(T*, T*, const int, const int, const bool)>;

    /**
     * @brief Constructor.
     * @param basis_type The basis type ("pw" or "lcao")
     * @param calculation The calculation type
     */
    DiagoPPCG(const std::string& basis_type, const std::string& calculation);

    /**
     * @brief Constructor with subspace function support.
     * @param basis_type The basis type ("pw" or "lcao")
     * @param calculation The calculation type
     * @param need_subspace Whether to use subspace diagonalization
     * @param subspace_func Subspace diagonalization function
     * @param pw_diag_thr Convergence threshold
     * @param pw_diag_nmax Max iterations
     * @param nproc_in_pool Number of processors in pool
     */
    DiagoPPCG(const std::string& basis_type,
              const std::string& calculation,
              const bool& need_subspace,
              const SubspaceFunc& subspace_func,
              const Real& pw_diag_thr,
              const int& pw_diag_nmax,
              const int& nproc_in_pool);

    /// Destructor
    ~DiagoPPCG();

    /**
     * @brief Solve eigenvalue problem using PPCG method.
     *
     * @param hpsi_func Function to compute H|psi>
     * @param spsi_func Function to compute S|psi> (identity for norm-conserving)
     * @param ld_psi Leading dimension of psi matrix
     * @param nband Number of bands
     * @param dim Dimension of basis
     * @param psi_in Input/output wavefunction matrix [dim: nband x ld_psi]
     * @param eigenvalue_in Output eigenvalues [dim: nband]
     * @param ethr_band Convergence threshold per band
     * @param prec Preconditioner array [dim: dim]
     * @return Average number of iterations
     */
    double diag(const HPsiFunc& hpsi_func,
                const SPsiFunc& spsi_func,
                const int ld_psi,
                const int nband,
                const int dim,
                T* psi_in,
                Real* eigenvalue_in,
                const std::vector<double>& ethr_band,
                const Real* prec = nullptr);

  private:
    Device* ctx_ = {};
    int n_band_ = 0;
    int n_basis_ = 0;
    int ld_psi_ = 0;

    double avg_iter_ = 0.0;
    int notconv_ = 0;
    int pw_diag_nmax_ = 0;
    Real pw_diag_thr_ = 1e-5;

    std::string basis_type_;
    std::string calculation_;

    HPsiFunc hpsi_func_ = nullptr;
    SPsiFunc spsi_func_ = nullptr;
    SubspaceFunc subspace_func_ = nullptr;

    bool need_subspace_ = false;
    int nproc_in_pool_ = 0;

    /// Constants
    const T* one_ = nullptr;
    const T* zero_ = nullptr;
    const T* neg_one_ = nullptr;

    /**
     * @brief Project gradient to be orthogonal to all lower states.
     */
    void orth_gradient(const ct::Tensor& psi,
                       const int& m,
                       ct::Tensor& grad,
                       ct::Tensor& s_grad,
                       ct::Tensor& lagrange);

    /**
     * @brief Compute gradient: g = H|ψ> - λ*S|ψ>, precondition and project.
     */
    void compute_gradient(const ct::Tensor& prec,
                          const ct::Tensor& hpsi,
                          const ct::Tensor& spsi,
                          const Real& eigenvalue,
                          ct::Tensor& grad,
                          ct::Tensor& ppsi);

    /**
     * @brief Update CG direction using Polak-Ribiere formula.
     * @return gamma coefficient for correction term
     */
    Real update_cg_direction(ct::Tensor& cg,
                             ct::Tensor& scg,
                             const ct::Tensor& grad,
                             const ct::Tensor& s_grad,
                             const ct::Tensor& prec,
                             ct::Tensor& g0,
                             Real& gg_last,
                             const int& iter);

    /**
     * @brief Perform line minimization in [ψ, cg] subspace and update ψ.
     * @return true if converged
     */
    bool line_minimization(const ct::Tensor& ppsi,
                           const ct::Tensor& cg,
                           const ct::Tensor& scg,
                           const double& ethreshold,
                           Real& cg_norm,
                           Real& theta,
                           Real& eigenvalue,
                           ct::Tensor& psi_m,
                           ct::Tensor& spsi_m,
                           ct::Tensor& hpsi_m);

    /**
     * @brief Schmidt orthogonalize ψ_m against previously converged states.
     */
    void schmidt_orth(const int& m,
                      const ct::Tensor& psi,
                      const ct::Tensor& spsi_m,
                      ct::Tensor& psi_m);

    bool test_exit_cond(const int& ntry, const int& notconv) const;
};

} // namespace hsolver

#endif // DIAGO_PPCG_H_

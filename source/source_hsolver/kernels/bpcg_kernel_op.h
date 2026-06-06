#ifndef MODULE_HSOLVER_BPCG_KERNEL_H
#define MODULE_HSOLVER_BPCG_KERNEL_H
#include "source_base/macros.h"
#include "source_base/module_device/types.h"
namespace hsolver
{/**
 * @brief Compute the optimal frequency for subspace diagonalization calls
 *        based on problem size (n_band * n_basis).
 *
 * Large problems benefit from more frequent subspace diagonalization to
 * maintain orthogonality, while small problems can reduce overhead by
 * calling it less frequently.
 *
 * @param n_band Number of bands (eigenvectors).
 * @param n_basis Number of basis functions.
 * @param nline Base iteration count for SCF > 1.
 * @return Optimal frequency (number of CG steps between calc_hsub_with_block calls).
 */
int compute_optimal_freq(const int n_band, const int n_basis, const int nline);

/**
 * @brief Compute the convergence rate from historical error values.
 *
 * The convergence rate is computed as the ratio of current average error
 * to previous average error. A rate close to 1.0 indicates slow convergence
 * (stagnation), while a rate close to 0.0 indicates fast convergence.
 *
 * @param prev_avg_err Average error from the previous measurement point.
 * @param curr_avg_err Average error from the current measurement point.
 * @return Convergence rate in range [0.0, 1.0+]. Values > 1.0 indicate divergence.
 */
double compute_convergence_rate(const double prev_avg_err, const double curr_avg_err);

/**
 * @brief Compute adaptive diagonalization frequency based on convergence rate.
 *
 * Uses the convergence rate to dynamically adjust how often subspace
 * diagonalization is performed. Slow convergence triggers more frequent
 * diagonalization to reset the search direction, while fast convergence
 * reduces the frequency to save computation.
 *
 * @param conv_rate Convergence rate (curr_avg_err / prev_avg_err).
 * @param base_freq Base frequency from problem-size-based computation.
 * @param nline Maximum allowed frequency (CG iteration limit).
 * @return Adaptive frequency in range [1, nline].
 */
int compute_adaptive_freq(const double conv_rate, const int base_freq, const int nline);


template <typename T, typename Device>
struct line_minimize_with_block_op
{
    /// @brief dot_real_op computes the dot product of the given complex
    /// arrays(treated as float arrays). And there's may have MPI communications
    /// while enabling planewave parallization strategy.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param A : input array arr
    /// \param dim : size of A
    /// \param lda : leading dimention of A
    /// \param batch : batch size, the size of the result array res
    ///
    /// \return res : the result vector
    /// T : dot product result
    void operator()(T* grad_out,
                    T* hgrad_out,
                    T* psi_out,
                    T* hpsi_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band);
};

template <typename T, typename Device>
struct calc_grad_with_block_op
{
    /// @brief dot_real_op computes the dot product of the given complex
    /// arrays(treated as float arrays). And there's may have MPI communications
    /// while enabling planewave parallization strategy.
    ///
    /// Input Parameters
    /// \param dev : the type of computing device
    /// \param A : input array arr
    /// \param dim : size of A
    /// \param lda : leading dimention of A
    /// \param batch : batch size, the size of the result array res
    ///
    /// \return res : the result vector
    /// T : dot product result
    using Real = typename GetTypeReal<T>::type;
    void operator()(const Real* prec_in,
                    Real* err_out,
                    Real* beta_out,
                    T* psi_out,
                    T* hpsi_out,
                    T* grad_out,
                    T* grad_old_out,
                    const int& n_basis,
                    const int& n_basis_max,
                    const int& n_band);
};

template <typename T, typename Device>
struct apply_eigenvalues_op
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& nbase,
                    const int& nbase_x,
                    const int& notconv,
                    T* result,
                    const T* vectors,
                    const Real* eigenvalues);
};

template <typename T, typename Device>
struct precondition_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   const Real* precondition,
                   const Real* eigenvalues);
};

template <typename T, typename Device>
struct normalize_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& dim,
                   T* psi_iter,
                   const int& nbase,
                   const int& notconv,
                   Real* psi_norm = nullptr);
};

template <typename T, typename Device> struct refresh_hcc_scc_vcc_op {
    using Real = typename GetTypeReal<T>::type;
  /// @brief refresh hcc scc vcc
  ///
  /// Input Parameters
  /// \param n : first dimension of matrix
  /// \param ldh : leading dimension of hcc, scc, vcc
  /// \param nbase : matrix size
  /// \param eigenvalue : input eigenvalue
  /// \param one : constant one
  ///
  /// Output Parameters
  /// \param hcc : output matrix hcc
  /// \param scc : output matrix scc
  /// \param vcc : output matrix vcc
  void operator()(const int &n,
                  T *hcc,
                  T *scc,
                  T *vcc,
                  const int &ldh,
                  const Real *eigenvalue,
                  const T& one);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

template <typename T>
struct line_minimize_with_block_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(T *grad_out, T *hgrad_out, T *psi_out, T *hpsi_out,
                  const int &n_basis, const int &n_basis_max,
                  const int &n_band);
};

template <typename T>
struct calc_grad_with_block_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const Real *prec_in, Real *err_out, Real *beta_out,
                  T *psi_out, T *hpsi_out, T *grad_out, T *grad_old_out,
                  const int &n_basis, const int &n_basis_max,
                  const int &n_band);
};

template <typename T>
struct apply_eigenvalues_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int& nbase,
                  const int& nbase_x,
                  const int& notconv,
                  T* result,
                  const T* vectors,
                  const Real* eigenvalues);
};

template <typename T>
struct precondition_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int& dim,
                 T* psi_iter,
                 const int& nbase,
                 const int& notconv,
                 const Real* precondition,
                 const Real* eigenvalues);
};

template <typename T>
struct normalize_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int& dim,
                 T* psi_iter,
                 const int& nbase,
                 const int& notconv,
                 Real* psi_norm = nullptr);
};

template <typename T>
struct refresh_hcc_scc_vcc_op<T, base_device::DEVICE_GPU> {
  using Real = typename GetTypeReal<T>::type;
  void operator()(const int &n,
                  T *hcc,
                  T *scc,
                  T *vcc,
                  const int &ldh,
                  const Real *eigenvalue,
                  const T& one);
};
#endif
} // namespace hsolver

#endif
#ifndef HSOLVER_EXP_LANCZOS_H
#define HSOLVER_EXP_LANCZOS_H
#include "source_hsolver/module_linalg/exp_action.h"
#include "source_hsolver/module_linalg/linear_workspace.h"

namespace hsolver
{
template <typename T, typename Device>
class ExpLanczos final : public ExpAction<T, Device>
{
  private:
    const diag_comm_info comm_;
    LinearWorkspace<T, Device> work_;
    ct::Tensor basis_;
    ct::Tensor image_;
    ct::Tensor packed_;
    ct::Tensor hpacked_;
    ct::Tensor coefficients_;
    ct::Tensor products_;
    ct::Tensor permutation_;
    std::vector<double> small_vectors_;
    std::vector<double> small_values_;
    std::vector<double> small_work_;
    int batch_size(const OrbLayout& layout) const;
    void combine(const OrbLayout& layout, int count, const T* basis, const std::vector<T>& coefficients, T* output);
    std::vector<std::complex<double>> basis_gram(const OrbLayout& layout, int count, const T* basis);
    bool diagonalize(const std::vector<double>& diagonal, const std::vector<double>& off_diagonal, int band, int bands, int count);
    ExpResult attempt(const HSOperator<T, Device>& op,
                      const OrbLayout& layout,
                      double alpha,
                      double rel_tol,
                      const T* input,
                      T* output,
                      T* full,
                      bool reorthogonalize,
                      bool* lost_orthogonality);
    ExpResult interval(const HSOperator<T, Device>& op,
                       const OrbLayout& layout,
                       double alpha,
                       double rel_tol,
                       const T* input,
                       T* output,
                       T* full,
                       int depth);
    ExpResult evaluate(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output, T* full);

  public:
    ExpLanczos(const ExpControls& exp_ctrl, const diag_comm_info& comm);
    ExpResult apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output) override;
    /** @brief Share one Krylov construction for alpha and twice alpha. */
    ExpResult apply_pair(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* half, T* full) override;
};
} // namespace hsolver
#endif

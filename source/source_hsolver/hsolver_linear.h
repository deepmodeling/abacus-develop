#ifndef HSOLVER_LINEAR_H
#define HSOLVER_LINEAR_H
#include "source_hsolver/linear_bicgstab.h"
#include "source_hsolver/linear_cgs.h"

#include <memory>

namespace hsolver
{
/** @brief Method selection and reusable numerical workspace, independent of physics. */
template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolverLinear
{
  private:
    std::unique_ptr<LinearBiCGSTAB<T, Device>> bicgstab_;
    std::unique_ptr<LinearCGS<T, Device>> cgs_;
    double tolerance_ = 0.0;

  public:
    HSolverLinear(const LinearSolveOptions& options, const diag_comm_info& comm);
    /** @brief Return the effective tolerance after resolving the precision-dependent default. */
    double tolerance() const
    {
        return tolerance_;
    }
    /** @brief Solve independent columns using the configured method and identity preconditioning. */
    LinearSolveResult solve(const LinearOperator<T, Device>& op, const int ld, const int nvec, const int dim, T* x, const T* b);
    /** @brief Solve with an explicit right inverse preconditioner. */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            const int ld,
                            const int nvec,
                            const int dim,
                            T* x,
                            const T* b);
};
} // namespace hsolver
#endif

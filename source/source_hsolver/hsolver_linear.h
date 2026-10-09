#ifndef HSOLVER_LINEAR_H
#define HSOLVER_LINEAR_H
#include "source_hsolver/linear_bicgstab.h"
#include "source_hsolver/linear_cgs.h"
#include "source_hsolver/linear_gmres.h"

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
    std::unique_ptr<LinearGMRES<T, Device>> gmres_;
    double tolerance_ = 0.0;
    const LinearSolveControl default_control_;

  public:
    HSolverLinear(const LinearSolveOptions& options, const diag_comm_info& comm);
    /** @brief Return the effective tolerance after resolving the precision-dependent default. */
    double tolerance() const
    {
        return tolerance_;
    }
    /** @brief Reuse an optional CN residual and control independent GMRES verification. */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            int ld,
                            int nvec,
                            int dim,
                            T* x,
                            const T* b,
                            const T* initial_residual,
                            bool force_check);
    /** @brief Reuse solver storage with an explicit per-call budget and reconstruction policy. */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            int ld,
                            int nvec,
                            int dim,
                            T* x,
                            const T* b,
                            const T* initial_residual,
                            bool force_check,
                            const LinearSolveControl& control);
};
} // namespace hsolver
#endif

#include "source_hsolver/module_linalg/hsolver_linear.h"

#include "source_base/tool_quit.h"

namespace hsolver
{
const char* linear_status_name(const LinearSolveStatus status)
{
    switch (status)
    {
    case LinearSolveStatus::converged:
        return "converged";
    case LinearSolveStatus::breakdown:
        return "breakdown";
    case LinearSolveStatus::max_iterations:
        return "maximum iterations";
    case LinearSolveStatus::residual_mismatch:
        return "true residual check failed";
    case LinearSolveStatus::preconditioner_failure:
        return "preconditioner application failed";
    }
    return "unknown linear solve status";
}

template <typename T, typename Device>
HSolverLinear<T, Device>::HSolverLinear(const LinearControls& linear_ctrl, const diag_comm_info& comm)
    : max_iter_(linear_ctrl.max_iterations), reconstruct_(linear_ctrl.reconstruct)
{
    using Real = typename GetTypeReal<T>::type;
    const double roundoff_tolerance = 100.0 * std::numeric_limits<Real>::epsilon();
    tolerance_ = linear_ctrl.tolerance == 0.0 ? std::max(1e-10, roundoff_tolerance) : linear_ctrl.tolerance;
    if (linear_ctrl.solver == LinearMethod::bicgstab)
    {
        bicgstab_.reset(new LinearBiCGSTAB<T, Device>(tolerance_, comm));
    }
    else if (linear_ctrl.solver == LinearMethod::cgs)
    {
        cgs_.reset(new LinearCGS<T, Device>(tolerance_, comm));
    }
    else if (linear_ctrl.solver == LinearMethod::gmres)
    {
        gmres_.reset(new LinearGMRES<T, Device>(tolerance_, linear_ctrl.restart, linear_ctrl.max_iterations, comm));
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLinear", "Unsupported linear solver method.");
    }
}

template <typename T, typename Device>
LinearSolveResult HSolverLinear<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                  const LinearOperator<T, Device>& preconditioner,
                                                  int ld,
                                                  int nvec,
                                                  int dim,
                                                  T* x,
                                                  const T* b,
                                                  const T* initial_residual,
                                                  bool force_check)
{
    return solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, force_check, max_iter_, reconstruct_);
}

template <typename T, typename Device>
LinearSolveResult HSolverLinear<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                  const LinearOperator<T, Device>& preconditioner,
                                                  int ld,
                                                  int nvec,
                                                  int dim,
                                                  T* x,
                                                  const T* b,
                                                  const T* initial_residual,
                                                  bool force_check,
                                                  int max_iter,
                                                  bool reconstruct)
{
    if (bicgstab_)
    {
        return bicgstab_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, max_iter);
    }
    else if (cgs_)
    {
        return cgs_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, max_iter);
    }
    else if (gmres_)
    {
        return gmres_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, force_check, max_iter, reconstruct);
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLinear::solve", "No linear solver has been initialized.");
    }
}

template class HSolverLinear<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverLinear<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class HSolverLinear<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverLinear<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver

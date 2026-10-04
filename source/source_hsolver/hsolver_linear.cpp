#include "source_hsolver/hsolver_linear.h"

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
HSolverLinear<T, Device>::HSolverLinear(const LinearSolveOptions& options, const diag_comm_info& comm)
    : default_control_{options.max_iterations, options.reconstruct}
{
    using Real = typename GetTypeReal<T>::type;
    tolerance_ = options.tolerance == 0.0 ? std::max(1e-10, 100.0 * std::numeric_limits<Real>::epsilon()) : options.tolerance;
    if (options.method == LinearMethod::bicgstab)
    {
        bicgstab_.reset(new LinearBiCGSTAB<T, Device>(tolerance_, comm));
    }
    else if (options.method == LinearMethod::cgs)
    {
        cgs_.reset(new LinearCGS<T, Device>(tolerance_, comm));
    }
    else if (options.method == LinearMethod::gmres)
    {
        gmres_.reset(new LinearGMRES<T, Device>(tolerance_, options, comm));
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLinear", "Unsupported linear solver method.");
    }
}

template <typename T, typename Device>
LinearSolveResult HSolverLinear<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                  const LinearOperator<T, Device>& preconditioner,
                                                  const int ld,
                                                  const int nvec,
                                                  const int dim,
                                                  T* x,
                                                  const T* b)
{
    return solve(op, preconditioner, ld, nvec, dim, x, b, nullptr, true);
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
    return solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, force_check, default_control_);
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
                                                  const LinearSolveControl& control)
{
    if (bicgstab_)
    {
        return bicgstab_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, control.max_iterations);
    }
    else if (cgs_)
    {
        return cgs_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, control.max_iterations);
    }
    else if (gmres_)
    {
        return gmres_->solve(op, preconditioner, ld, nvec, dim, x, b, initial_residual, force_check, control);
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

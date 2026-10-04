#include "source_hsolver/hsolver_linear.h"

#include "source_base/module_device/memory_op.h"
#include "source_base/tool_quit.h"

namespace hsolver
{
namespace
{
template <typename T, typename Device>
class IdentityOperator final : public LinearOperator<T, Device>
{
  private:
    const int dim_;

  public:
    explicit IdentityOperator(const int dim) : dim_(dim)
    {
    }
    bool is_identity() const override
    {
        return true;
    }
    void apply(const T* x, T* y, const int ld, const int nvec) const override
    {
        if (dim_ > 0)
        {
            base_device::memory::synchronize_memory_2d_op<T, Device, Device>()(y, ld, x, ld, dim_, nvec);
        }
    }
};
} // namespace

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
    }
    return "unknown linear solve status";
}

template <typename T, typename Device>
HSolverLinear<T, Device>::HSolverLinear(const LinearSolveOptions& options, const diag_comm_info& comm)
{
    using Real = typename GetTypeReal<T>::type;
    tolerance_ = options.tolerance == 0.0 ? std::max(1e-10, 100.0 * std::numeric_limits<Real>::epsilon()) : options.tolerance;
    if (options.method == LinearMethod::bicgstab)
    {
        bicgstab_.reset(new LinearBiCGSTAB<T, Device>(tolerance_, options.max_iterations, comm));
    }
    else if (options.method == LinearMethod::cgs)
    {
        cgs_.reset(new LinearCGS<T, Device>(tolerance_, options.max_iterations, comm));
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLinear", "Unsupported linear solver method.");
    }
}

template <typename T, typename Device>
LinearSolveResult HSolverLinear<T, Device>::solve(const LinearOperator<T, Device>& op,
                                                  const int ld,
                                                  const int nvec,
                                                  const int dim,
                                                  T* x,
                                                  const T* b)
{
    const IdentityOperator<T, Device> identity(dim);
    return solve(op, identity, ld, nvec, dim, x, b);
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
    if (bicgstab_)
    {
        return bicgstab_->solve(op, preconditioner, ld, nvec, dim, x, b);
    }
    else if (cgs_)
    {
        return cgs_->solve(op, preconditioner, ld, nvec, dim, x, b);
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

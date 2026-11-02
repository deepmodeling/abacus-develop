#ifndef HSOLVER_LINEAR_OPERATOR_H
#define HSOLVER_LINEAR_OPERATOR_H

#include "source_base/module_device/types.h"

#include <exception>

namespace hsolver
{

/** @brief Stop the current recurrence when a preconditioner cannot be applied. */
class LinearPreconditionerError final : public std::exception
{
  public:
    const char* what() const noexcept override
    {
        return "Preconditioner application failed.";
    }
};

/**
 * @brief Device-resident block operator for independent right-hand sides.
 *
 * Column i starts at x + i*ld. The operator knows the valid row count.
 * Inputs are read-only; input/output buffers must not overlap. The same
 * operator acts on every column, including reordered or compacted blocks.
 * A zero local row count must still participate in operator collectives.
 */
template <typename T, typename Device = base_device::DEVICE_CPU>
class LinearOperator
{
  public:
    virtual ~LinearOperator() = default;
    virtual void apply(const T* x, T* y, const int ld, const int nvec) const = 0;
    /** @brief Allow solvers to bypass identity preconditioner copies. */
    virtual bool is_identity() const { return false; }
};

} // namespace hsolver
#endif

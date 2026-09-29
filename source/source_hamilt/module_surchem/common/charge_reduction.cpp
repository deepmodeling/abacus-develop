#include "charge_reduction.h"

#include "source_base/parallel_reduce.h"

#include <cmath>
#include <stdexcept>

namespace ModuleSurchem
{

void SerialChargeReduction::reduce_sum(double& value) const
{
    if (!std::isfinite(value))
    {
        throw std::domain_error("SCCS charge integral must be finite before reduction");
    }
}

void SerialChargeReduction::reduce_sum(double* values, const int count) const
{
    if (count < 0 || (count > 0 && values == nullptr))
    {
        throw std::invalid_argument("SCCS charge-array reduction requires valid storage and size");
    }
    for (int index = 0; index < count; ++index)
    {
        if (!std::isfinite(values[index]))
        {
            throw std::domain_error("SCCS charge array must be finite before reduction");
        }
    }
}

void SerialChargeReduction::reduce_max(double& value) const
{
    if (!std::isfinite(value))
    {
        throw std::domain_error("SCCS grid maximum must be finite before reduction");
    }
}

PoolChargeReduction::PoolChargeReduction(const int process_count)
    : process_count_(process_count)
{
    if (process_count_ <= 0)
    {
        throw std::invalid_argument("SCCS pool reduction requires a positive process count");
    }
}

void PoolChargeReduction::reduce_sum(double& value) const
{
    Parallel_Reduce::reduce_pool(value);
}

void PoolChargeReduction::reduce_sum(double* values, const int count) const
{
    if (count < 0 || (count > 0 && values == nullptr))
    {
        throw std::invalid_argument("SCCS pool reduction requires valid array storage and size");
    }
    if (count == 0)
    {
        return;
    }
    Parallel_Reduce::reduce_pool(values, count);
}

void PoolChargeReduction::reduce_max(double& value) const
{
    Parallel_Reduce::reduce_max_pool(process_count_, value);
}

} // namespace ModuleSurchem

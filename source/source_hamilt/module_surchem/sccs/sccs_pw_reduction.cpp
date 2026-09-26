#include "sccs_pw_reduction.h"

#include "source_base/parallel_reduce.h"

#include <stdexcept>

namespace ModuleSccs
{

PoolPolarizationReduction::PoolPolarizationReduction(const int process_count)
    : process_count_(process_count)
{
    if (process_count_ <= 0)
    {
        throw std::invalid_argument("SCCS pool reduction requires a positive process count");
    }
}

void PoolPolarizationReduction::reduce_residual(double& square_sum,
                                                double& maximum,
                                                double& point_count) const
{
    Parallel_Reduce::reduce_pool(square_sum);
    Parallel_Reduce::reduce_max_pool(process_count_, maximum);
    Parallel_Reduce::reduce_pool(point_count);
}

void PoolPolarizationReduction::reduce_sum(double& value) const
{
    Parallel_Reduce::reduce_pool(value);
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

} // namespace ModuleSccs

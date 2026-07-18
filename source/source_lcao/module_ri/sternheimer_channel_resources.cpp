#include "source_lcao/module_ri/sternheimer_channel_resources.h"

#include <algorithm>
#include <complex>
#include <limits>
#include <stdexcept>

namespace ModuleRI
{
namespace
{

constexpr std::uint64_t kComplexVectorsPerWorker = 120;

std::uint64_t checked_multiply(const std::uint64_t lhs, const std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs)
    {
        throw std::overflow_error("Sternheimer channel memory estimate overflow.");
    }
    return lhs * rhs;
}

std::uint64_t three_quarters(const std::uint64_t bytes)
{
    return (bytes / 4) * 3 + ((bytes % 4) * 3) / 4;
}

std::uint64_t checked_remaining(const std::uint64_t target, const std::uint64_t current)
{
    if (current >= target)
    {
        throw std::runtime_error(
            "Sternheimer channel memory baseline leaves no room for one worker within the 75 percent target.");
    }
    return target - current;
}

} // namespace

std::uint64_t estimate_sternheimer_channel_worker_bytes(const std::size_t grid_size)
{
    return checked_multiply(checked_multiply(kComplexVectorsPerWorker, static_cast<std::uint64_t>(grid_size)),
                            sizeof(std::complex<double>));
}

SternheimerChannelWorkerPlan plan_sternheimer_channel_workers(const int num_channels,
                                                              const int omp_threads,
                                                              const std::size_t grid_size,
                                                              const int user_cap,
                                                              const SternheimerMemorySnapshot& memory)
{
    if (num_channels <= 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a positive channel count.");
    }
    if (omp_threads <= 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a positive OpenMP thread count.");
    }
    if (grid_size == 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a non-empty real-space grid.");
    }
    if (user_cap < 0)
    {
        throw std::invalid_argument("Sternheimer channel worker user cap must be non-negative.");
    }
    if (memory.local_mpi_ranks <= 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a positive local MPI rank count.");
    }

    SternheimerChannelWorkerPlan plan;
    plan.memory_per_worker_bytes = estimate_sternheimer_channel_worker_bytes(grid_size);
    if (memory.mode == SternheimerMemoryAccountingMode::fallback_one)
    {
        return plan;
    }
    if (memory.limit_bytes == 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a positive detected memory value.");
    }

    switch (memory.mode)
    {
        case SternheimerMemoryAccountingMode::node_aggregate:
            plan.target_bytes = three_quarters(memory.limit_bytes);
            plan.increment_bytes_per_rank
                = checked_remaining(plan.target_bytes, memory.current_bytes)
                  / static_cast<std::uint64_t>(memory.local_mpi_ranks);
            break;
        case SternheimerMemoryAccountingMode::per_rank:
        {
            const std::uint64_t rank_limit
                = memory.limit_bytes / static_cast<std::uint64_t>(memory.local_mpi_ranks);
            plan.target_bytes = three_quarters(rank_limit);
            plan.increment_bytes_per_rank = checked_remaining(plan.target_bytes, memory.current_bytes);
            break;
        }
        case SternheimerMemoryAccountingMode::available:
            plan.target_bytes = three_quarters(memory.limit_bytes);
            plan.increment_bytes_per_rank
                = plan.target_bytes / static_cast<std::uint64_t>(memory.local_mpi_ranks);
            break;
        case SternheimerMemoryAccountingMode::fallback_one:
            break;
    }

    const std::uint64_t memory_worker_count = plan.increment_bytes_per_rank / plan.memory_per_worker_bytes;
    if (memory_worker_count == 0)
    {
        throw std::runtime_error(
            "Sternheimer channel memory budget cannot accommodate one worker within the 75 percent target.");
    }
    const std::uint64_t bounded_memory_workers
        = std::min(memory_worker_count, static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    plan.automatic_workers
        = std::min({num_channels, omp_threads, static_cast<int>(bounded_memory_workers)});
    plan.effective_workers = user_cap > 0 ? std::min(plan.automatic_workers, user_cap) : plan.automatic_workers;
    return plan;
}

} // namespace ModuleRI

#ifndef STERNHEIMER_CHANNEL_RESOURCES_H
#define STERNHEIMER_CHANNEL_RESOURCES_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace ModuleRI
{

enum class SternheimerMemoryAccountingMode
{
    node_aggregate,
    per_rank,
    available,
    fallback_one
};

struct SternheimerMemorySnapshot
{
    SternheimerMemoryAccountingMode mode = SternheimerMemoryAccountingMode::fallback_one;
    std::uint64_t limit_bytes = 0;
    std::uint64_t current_bytes = 0;
    int local_mpi_ranks = 1;
    std::string source = "unavailable";
};

struct SternheimerChannelWorkerPlan
{
    int automatic_workers = 1;
    int effective_workers = 1;
    std::uint64_t target_bytes = 0;
    std::uint64_t increment_bytes_per_rank = 0;
    std::uint64_t memory_per_worker_bytes = 0;
};

std::uint64_t estimate_sternheimer_channel_worker_bytes(std::size_t grid_size);

SternheimerChannelWorkerPlan plan_sternheimer_channel_workers(int num_channels,
                                                              int omp_threads,
                                                              std::size_t grid_size,
                                                              int user_cap,
                                                              const SternheimerMemorySnapshot& memory);

} // namespace ModuleRI

#endif

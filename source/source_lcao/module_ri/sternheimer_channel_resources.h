#ifndef STERNHEIMER_CHANNEL_RESOURCES_H
#define STERNHEIMER_CHANNEL_RESOURCES_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

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
    int channel_batch_width = 1;
    int batch_tasks = 1;
};

struct SternheimerChannelBatch
{
    int begin = 0;
    int size = 0;
};

std::vector<SternheimerChannelBatch> make_sternheimer_channel_batches(int num_channels, int batch_width);

std::uint64_t estimate_sternheimer_channel_worker_bytes(std::size_t grid_size);

SternheimerChannelWorkerPlan plan_sternheimer_channel_workers(int num_channels,
                                                              int omp_threads,
                                                              std::size_t grid_size,
                                                              int user_cap,
                                                              const SternheimerMemorySnapshot& memory,
                                                              int channel_batch_width = 1);

SternheimerChannelWorkerPlan plan_sternheimer_owned_channel_workers(int global_num_channels,
                                                                    int owned_num_channels,
                                                                    int omp_threads,
                                                                    std::size_t grid_size,
                                                                    int user_cap,
                                                                    const SternheimerMemorySnapshot& memory,
                                                                    int channel_batch_width = 1);

SternheimerMemorySnapshot detect_sternheimer_memory_snapshot();

void initialize_sternheimer_memory_detection();

std::string sternheimer_memory_accounting_mode_name(SternheimerMemoryAccountingMode mode);

std::string format_sternheimer_channel_worker_diagnostic(const SternheimerMemorySnapshot& memory,
                                                         const SternheimerChannelWorkerPlan& plan,
                                                         std::size_t grid_size,
                                                         int user_cap);

namespace detail
{

template <typename T>
class SternheimerOptionalValue
{
  public:
    SternheimerOptionalValue() = default;
    SternheimerOptionalValue(const T& value) : has_value_(true), value_(value)
    {
    }

    SternheimerOptionalValue& operator=(const T& value)
    {
        has_value_ = true;
        value_ = value;
        return *this;
    }

    explicit operator bool() const
    {
        return has_value_;
    }

    bool has_value() const
    {
        return has_value_;
    }

    const T& operator*() const
    {
        return value_;
    }

    T& operator*()
    {
        return value_;
    }

    const T* operator->() const
    {
        return &value_;
    }

    T* operator->()
    {
        return &value_;
    }

    void reset()
    {
        has_value_ = false;
        value_ = T{};
    }

    template <typename U>
    bool operator==(const U& rhs) const
    {
        return has_value_ && value_ == rhs;
    }

  private:
    bool has_value_ = false;
    T value_{};
};

struct SternheimerMemoryCandidates
{
    SternheimerOptionalValue<std::uint64_t> cgroup_limit_bytes;
    SternheimerOptionalValue<std::uint64_t> cgroup_current_bytes;
    SternheimerOptionalValue<std::uint64_t> slurm_limit_bytes;
    SternheimerOptionalValue<std::uint64_t> mem_available_bytes;
    SternheimerOptionalValue<std::uint64_t> process_rss_bytes;
    SternheimerOptionalValue<std::uint64_t> physical_memory_bytes;
    std::string cgroup_source;
};

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_memory_bytes(const std::string& text);

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_slurm_mem_per_node(const std::string& text);

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_kib_field(const std::string& text, const std::string& key);

SternheimerOptionalValue<std::string> parse_sternheimer_cgroup_v2_path(const std::string& text);

SternheimerOptionalValue<std::string> parse_sternheimer_cgroup_v1_memory_path(const std::string& text);

SternheimerMemorySnapshot select_sternheimer_memory_snapshot(const SternheimerMemoryCandidates& candidates,
                                                             int local_mpi_ranks);

} // namespace detail

} // namespace ModuleRI

#endif

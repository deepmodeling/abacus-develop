#include "source_lcao/module_ri/sternheimer_channel_resources.h"

#include <algorithm>
#include <cctype>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleRI
{
namespace
{

constexpr std::uint64_t kComplexVectorsPerWorker = 120;

#ifdef __MPI
MPI_Comm sternheimer_local_communicator()
{
    static MPI_Comm local_communicator = MPI_COMM_NULL;
    if (local_communicator == MPI_COMM_NULL)
    {
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_communicator);
    }
    // Keep the communicator alive until MPI_Finalize so late memory probes do not
    // repeat topology discovery after the large response grids are allocated.
    return local_communicator;
}
#endif

std::uint64_t checked_multiply(const std::uint64_t lhs, const std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs)
    {
        throw std::overflow_error("Sternheimer channel memory estimate overflow.");
    }
    return lhs * rhs;
}

std::uint64_t checked_add(const std::uint64_t lhs, const std::uint64_t rhs)
{
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs)
    {
        throw std::overflow_error("Sternheimer shared memory reservation overflow.");
    }
    return lhs + rhs;
}

std::uint64_t three_quarters(const std::uint64_t bytes)
{
    return (bytes / 4) * 3 + ((bytes % 4) * 3) / 4;
}

std::runtime_error memory_budget_error(const char* reason,
                                      const SternheimerMemorySnapshot& memory,
                                      const std::uint64_t target)
{
    std::ostringstream message;
    message << reason << " resource_source=" << memory.source
            << " accounting_mode=" << sternheimer_memory_accounting_mode_name(memory.mode)
            << " node_memory_limit_bytes=" << memory.limit_bytes
            << " memory_current_bytes=" << memory.current_bytes
            << " memory_target_bytes=" << target << " local_mpi_ranks=" << memory.local_mpi_ranks;
    return std::runtime_error(message.str());
}

std::uint64_t checked_remaining(const std::uint64_t target, const SternheimerMemorySnapshot& memory)
{
    if (memory.current_bytes >= target)
    {
        throw memory_budget_error(
            "Sternheimer channel memory baseline leaves no room for one worker within the 75 percent target.",
            memory,
            target);
    }
    return target - memory.current_bytes;
}

std::string trim_copy(const std::string& text)
{
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first])))
    {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1])))
    {
        --last;
    }
    return text.substr(first, last - first);
}

detail::SternheimerOptionalValue<std::uint64_t> parse_unsigned_decimal(const std::string& text)
{
    if (text.empty())
    {
        return {};
    }
    std::uint64_t value = 0;
    for (const char character: text)
    {
        if (!std::isdigit(static_cast<unsigned char>(character)))
        {
            return {};
        }
        const std::uint64_t digit = static_cast<std::uint64_t>(character - '0');
        if (value > (std::numeric_limits<std::uint64_t>::max() - digit) / 10)
        {
            return {};
        }
        value = value * 10 + digit;
    }
    return value;
}

detail::SternheimerOptionalValue<std::string> read_text_file(const std::string& path)
{
    std::ifstream input(path);
    if (!input)
    {
        return {};
    }
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

std::string cgroup_file_path(const std::string& mount, const std::string& relative, const std::string& filename)
{
    std::string path = mount;
    if (!relative.empty() && relative != "/")
    {
        if (relative.front() != '/')
        {
            path += '/';
        }
        path += relative;
    }
    path += '/';
    path += filename;
    return path;
}

void read_live_memory_candidates(detail::SternheimerMemoryCandidates& candidates)
{
    const detail::SternheimerOptionalValue<std::string> proc_cgroup = read_text_file("/proc/self/cgroup");
    if (proc_cgroup)
    {
        const detail::SternheimerOptionalValue<std::string> v2_path
            = detail::parse_sternheimer_cgroup_v2_path(*proc_cgroup);
        if (v2_path)
        {
            const detail::SternheimerOptionalValue<std::string> limit
                = read_text_file(cgroup_file_path("/sys/fs/cgroup", *v2_path, "memory.max"));
            const detail::SternheimerOptionalValue<std::string> current
                = read_text_file(cgroup_file_path("/sys/fs/cgroup", *v2_path, "memory.current"));
            if (limit)
            {
                candidates.cgroup_limit_bytes = detail::parse_sternheimer_memory_bytes(*limit);
            }
            if (current)
            {
                candidates.cgroup_current_bytes = detail::parse_sternheimer_memory_bytes(*current);
            }
            candidates.cgroup_source = "cgroup_v2";
        }
        else
        {
            const detail::SternheimerOptionalValue<std::string> v1_path
                = detail::parse_sternheimer_cgroup_v1_memory_path(*proc_cgroup);
            if (v1_path)
            {
                const std::string mount = "/sys/fs/cgroup/memory";
                const detail::SternheimerOptionalValue<std::string> limit
                    = read_text_file(cgroup_file_path(mount, *v1_path, "memory.limit_in_bytes"));
                const detail::SternheimerOptionalValue<std::string> current
                    = read_text_file(cgroup_file_path(mount, *v1_path, "memory.usage_in_bytes"));
                if (limit)
                {
                    candidates.cgroup_limit_bytes = detail::parse_sternheimer_memory_bytes(*limit);
                }
                if (current)
                {
                    candidates.cgroup_current_bytes = detail::parse_sternheimer_memory_bytes(*current);
                }
                candidates.cgroup_source = "cgroup_v1";
            }
        }
    }

    if (const char* const slurm_memory = std::getenv("SLURM_MEM_PER_NODE"))
    {
        candidates.slurm_limit_bytes = detail::parse_sternheimer_slurm_mem_per_node(slurm_memory);
    }
    const detail::SternheimerOptionalValue<std::string> meminfo = read_text_file("/proc/meminfo");
    if (meminfo)
    {
        candidates.mem_available_bytes = detail::parse_sternheimer_kib_field(*meminfo, "MemAvailable");
        candidates.physical_memory_bytes = detail::parse_sternheimer_kib_field(*meminfo, "MemTotal");
    }
    const detail::SternheimerOptionalValue<std::string> status = read_text_file("/proc/self/status");
    if (status)
    {
        candidates.process_rss_bytes = detail::parse_sternheimer_kib_field(*status, "VmRSS");
    }
}

#ifdef __MPI
detail::SternheimerOptionalValue<std::uint64_t> allreduce_optional_min(
    const detail::SternheimerOptionalValue<std::uint64_t>& value,
    MPI_Comm communicator)
{
    const unsigned long long local
        = value ? static_cast<unsigned long long>(*value) : std::numeric_limits<unsigned long long>::max();
    unsigned long long reduced = 0;
    MPI_Allreduce(&local, &reduced, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, communicator);
    if (reduced == std::numeric_limits<unsigned long long>::max())
    {
        return {};
    }
    return static_cast<std::uint64_t>(reduced);
}

detail::SternheimerOptionalValue<std::uint64_t> allreduce_optional_max(
    const detail::SternheimerOptionalValue<std::uint64_t>& value,
    MPI_Comm communicator)
{
    const unsigned long long local = value ? static_cast<unsigned long long>(*value) : 0ULL;
    unsigned long long reduced = 0;
    MPI_Allreduce(&local, &reduced, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, communicator);
    if (reduced == 0)
    {
        return {};
    }
    return static_cast<std::uint64_t>(reduced);
}
#endif

} // namespace

SternheimerHeapTrimStatus trim_sternheimer_process_heap() noexcept
{
#if defined(__GLIBC__)
    return malloc_trim(0) != 0 ? SternheimerHeapTrimStatus::pages_released
                               : SternheimerHeapTrimStatus::no_pages_released;
#else
    return SternheimerHeapTrimStatus::unsupported;
#endif
}

SternheimerMemorySnapshot reserve_sternheimer_shared_memory(const SternheimerMemorySnapshot& memory,
                                                           const std::uint64_t bytes_per_rank)
{
    if (memory.mode == SternheimerMemoryAccountingMode::fallback_one)
    {
        return memory;
    }
    if (memory.local_mpi_ranks <= 0)
    {
        throw std::invalid_argument("Sternheimer shared memory reservation requires a positive local MPI rank count.");
    }
    SternheimerMemorySnapshot adjusted = memory;
    switch (memory.mode)
    {
    case SternheimerMemoryAccountingMode::node_aggregate:
        adjusted.current_bytes = checked_add(
            memory.current_bytes, checked_multiply(bytes_per_rank, static_cast<std::uint64_t>(memory.local_mpi_ranks)));
        break;
    case SternheimerMemoryAccountingMode::per_rank:
        adjusted.current_bytes = checked_add(memory.current_bytes, bytes_per_rank);
        break;
    case SternheimerMemoryAccountingMode::available: {
        const std::uint64_t reserved_bytes
            = checked_multiply(bytes_per_rank, static_cast<std::uint64_t>(memory.local_mpi_ranks));
        if (reserved_bytes > memory.limit_bytes)
        {
            throw memory_budget_error("Sternheimer shared memory reservation exceeds available memory.",
                                      memory,
                                      three_quarters(memory.limit_bytes));
        }
        adjusted.limit_bytes -= reserved_bytes;
        break;
    }
    case SternheimerMemoryAccountingMode::fallback_one:
        break;
    }
    return adjusted;
}

std::uint64_t estimate_sternheimer_channel_worker_bytes(const std::size_t grid_size)
{
    return checked_multiply(checked_multiply(kComplexVectorsPerWorker, static_cast<std::uint64_t>(grid_size)),
                            sizeof(std::complex<double>));
}

std::vector<SternheimerChannelBatch> make_sternheimer_channel_batches(const int num_channels, const int batch_width)
{
    if (num_channels < 0)
    {
        throw std::invalid_argument("Sternheimer channel batching requires a non-negative channel count.");
    }
    if (batch_width <= 0)
    {
        throw std::invalid_argument("Sternheimer channel batching requires a positive batch width.");
    }
    std::vector<SternheimerChannelBatch> batches;
    for (int begin = 0; begin < num_channels; begin += batch_width)
    {
        batches.push_back({begin, std::min(batch_width, num_channels - begin)});
    }
    return batches;
}

std::vector<SternheimerChannelBatch> make_sternheimer_owned_channel_batches(const int num_channels,
                                                                           const int ownership_batch_width,
                                                                           const int worker_batch_width,
                                                                           const int replica_count,
                                                                           const int replica_index,
                                                                           const int occupied_band_index)
{
    if (num_channels < 0 || ownership_batch_width <= 0 || worker_batch_width <= 0
        || replica_count <= 0 || replica_index < 0 || replica_index >= replica_count
        || occupied_band_index < 0)
    {
        throw std::invalid_argument("Invalid Sternheimer owned-channel batch dimensions.");
    }
    std::vector<SternheimerChannelBatch> batches;
    if (num_channels == 0)
    {
        return batches;
    }
    const std::int64_t chunk_count = 1 + (num_channels - 1) / ownership_batch_width;
    // Match sternheimer_channel_batch_replica_owner using the global requested-width count.
    const std::int64_t band_offset = static_cast<std::int64_t>(occupied_band_index) * chunk_count % replica_count;
    const std::int64_t first_chunk = (replica_index - band_offset + replica_count) % replica_count;
    for (std::int64_t chunk = first_chunk; chunk < chunk_count; chunk += replica_count)
    {
        int begin = static_cast<int>(chunk * ownership_batch_width);
        const int end = begin + std::min(ownership_batch_width, num_channels - begin);
        while (begin < end)
        {
            const int size = std::min(worker_batch_width, end - begin);
            batches.push_back({begin, size});
            begin += size;
        }
    }
    return batches;
}

SternheimerChannelWorkerPlan plan_sternheimer_channel_workers(const int num_channels,
                                                              const int omp_threads,
                                                              const std::size_t grid_size,
                                                              const int user_cap,
                                                              const SternheimerMemorySnapshot& memory,
                                                              const int channel_batch_width)
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
    if (channel_batch_width <= 0)
    {
        throw std::invalid_argument("Sternheimer channel worker planning requires a positive batch width.");
    }

    SternheimerChannelWorkerPlan plan;
    const std::uint64_t base_worker_bytes = estimate_sternheimer_channel_worker_bytes(grid_size);
    plan.channel_batch_width = channel_batch_width;
    plan.batch_tasks = (num_channels + channel_batch_width - 1) / channel_batch_width;
    plan.memory_per_worker_bytes
        = checked_multiply(base_worker_bytes, static_cast<std::uint64_t>(channel_batch_width));
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
        plan.increment_bytes_per_rank = checked_remaining(plan.target_bytes, memory)
                                        / static_cast<std::uint64_t>(memory.local_mpi_ranks);
        break;
    case SternheimerMemoryAccountingMode::per_rank: {
        const std::uint64_t rank_limit = memory.limit_bytes / static_cast<std::uint64_t>(memory.local_mpi_ranks);
        plan.target_bytes = three_quarters(rank_limit);
        plan.increment_bytes_per_rank = checked_remaining(plan.target_bytes, memory);
        break;
    }
    case SternheimerMemoryAccountingMode::available:
        plan.target_bytes = three_quarters(memory.limit_bytes);
        plan.increment_bytes_per_rank = plan.target_bytes / static_cast<std::uint64_t>(memory.local_mpi_ranks);
        break;
    case SternheimerMemoryAccountingMode::fallback_one:
        break;
    }

    const std::uint64_t memory_batch_width = plan.increment_bytes_per_rank / base_worker_bytes;
    if (memory_batch_width == 0)
    {
        throw memory_budget_error(
            "Sternheimer channel memory budget cannot accommodate one worker within the 75 percent target.",
            memory,
            plan.target_bytes);
    }
    plan.channel_batch_width = std::min(
        channel_batch_width,
        static_cast<int>(std::min(memory_batch_width,
                                  static_cast<std::uint64_t>(std::numeric_limits<int>::max()))));
    plan.batch_tasks = (num_channels + plan.channel_batch_width - 1) / plan.channel_batch_width;
    plan.memory_per_worker_bytes
        = checked_multiply(base_worker_bytes, static_cast<std::uint64_t>(plan.channel_batch_width));
    const std::uint64_t memory_worker_count = plan.increment_bytes_per_rank / plan.memory_per_worker_bytes;
    const std::uint64_t bounded_memory_workers
        = std::min(memory_worker_count, static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    plan.automatic_workers
        = std::min({plan.batch_tasks, omp_threads, static_cast<int>(bounded_memory_workers)});
    const int capped_workers = user_cap > 0 ? std::min(plan.automatic_workers, user_cap) : plan.automatic_workers;
    // Use inner grid OpenMP only when independent channel tasks fill less than
    // one quarter of the available thread team.  Above that point the FD
    // stencil is memory-bandwidth bound and channel parallelism is preferable.
    const int minimum_channel_workers = (omp_threads + 3) / 4;
    plan.effective_workers = capped_workers >= minimum_channel_workers ? capped_workers : 1;
    return plan;
}

SternheimerChannelWorkerPlan plan_sternheimer_owned_channel_workers(const int global_num_channels,
                                                                    const int owned_num_channels,
                                                                    const int omp_threads,
                                                                    const std::size_t grid_size,
                                                                    const int user_cap,
                                                                    const SternheimerMemorySnapshot& memory,
                                                                    const int channel_batch_width)
{
    if (global_num_channels <= 0)
    {
        throw std::invalid_argument("Sternheimer owned-channel planning requires a positive global channel count.");
    }
    if (owned_num_channels <= 0 || owned_num_channels > global_num_channels)
    {
        throw std::invalid_argument(
            "Sternheimer owned-channel planning requires a positive local count no larger than the global count.");
    }
    return plan_sternheimer_channel_workers(owned_num_channels,
                                            omp_threads,
                                            grid_size,
                                            user_cap,
                                            memory,
                                            channel_batch_width);
}

std::string sternheimer_memory_accounting_mode_name(const SternheimerMemoryAccountingMode mode)
{
    switch (mode)
    {
    case SternheimerMemoryAccountingMode::node_aggregate:
        return "node_aggregate";
    case SternheimerMemoryAccountingMode::per_rank:
        return "per_rank";
    case SternheimerMemoryAccountingMode::available:
        return "available";
    case SternheimerMemoryAccountingMode::fallback_one:
        return "fallback_one";
    }
    throw std::invalid_argument("Unknown Sternheimer memory accounting mode.");
}

std::string format_sternheimer_channel_worker_diagnostic(const SternheimerMemorySnapshot& memory,
                                                         const SternheimerChannelWorkerPlan& plan,
                                                         const std::size_t grid_size,
                                                         const int user_cap)
{
    std::ostringstream diagnostic;
    diagnostic << "resource_source=" << memory.source
               << " accounting_mode=" << sternheimer_memory_accounting_mode_name(memory.mode)
               << " node_memory_limit_bytes=" << memory.limit_bytes << " memory_current_bytes=" << memory.current_bytes
               << " memory_target_bytes=" << plan.target_bytes << " local_mpi_ranks=" << memory.local_mpi_ranks
               << " grid_size=" << grid_size << " memory_per_worker_bytes=" << plan.memory_per_worker_bytes
               << " channel_batch_width=" << plan.channel_batch_width << " batch_tasks=" << plan.batch_tasks
               << " automatic_workers=" << plan.automatic_workers << " user_cap=" << user_cap
               << " effective_workers=" << plan.effective_workers;
    return diagnostic.str();
}

SternheimerMemorySnapshot detect_sternheimer_memory_snapshot()
{
    detail::SternheimerMemoryCandidates candidates;
    int local_mpi_ranks = 1;
#ifdef __MPI
    const MPI_Comm local_communicator = sternheimer_local_communicator();
    MPI_Comm_size(local_communicator, &local_mpi_ranks);
    MPI_Barrier(local_communicator);
    read_live_memory_candidates(candidates);
    candidates.cgroup_limit_bytes = allreduce_optional_min(candidates.cgroup_limit_bytes, local_communicator);
    candidates.cgroup_current_bytes = allreduce_optional_max(candidates.cgroup_current_bytes, local_communicator);
    candidates.slurm_limit_bytes = allreduce_optional_min(candidates.slurm_limit_bytes, local_communicator);
    candidates.mem_available_bytes = allreduce_optional_min(candidates.mem_available_bytes, local_communicator);
    candidates.process_rss_bytes = allreduce_optional_max(candidates.process_rss_bytes, local_communicator);
    candidates.physical_memory_bytes = allreduce_optional_min(candidates.physical_memory_bytes, local_communicator);
#else
    read_live_memory_candidates(candidates);
#endif
    return detail::select_sternheimer_memory_snapshot(candidates, local_mpi_ranks);
}

void initialize_sternheimer_memory_detection()
{
#ifdef __MPI
    (void)sternheimer_local_communicator();
#endif
}

namespace detail
{

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_memory_bytes(const std::string& input)
{
    std::string text = trim_copy(input);
    if (text.empty() || text == "max")
    {
        return {};
    }

    std::uint64_t multiplier = 1;
    const char suffix = text.back();
    if (!std::isdigit(static_cast<unsigned char>(suffix)))
    {
        switch (static_cast<char>(std::toupper(static_cast<unsigned char>(suffix))))
        {
        case 'K':
            multiplier = 1024ULL;
            break;
        case 'M':
            multiplier = 1024ULL * 1024ULL;
            break;
        case 'G':
            multiplier = 1024ULL * 1024ULL * 1024ULL;
            break;
        default:
            return {};
        }
        text.pop_back();
    }
    const SternheimerOptionalValue<std::uint64_t> value = parse_unsigned_decimal(text);
    if (!value)
    {
        return {};
    }
    try
    {
        return checked_multiply(*value, multiplier);
    }
    catch (const std::overflow_error&)
    {
        return {};
    }
}

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_slurm_mem_per_node(const std::string& input)
{
    const SternheimerOptionalValue<std::uint64_t> mebibytes = parse_unsigned_decimal(trim_copy(input));
    if (!mebibytes)
    {
        return {};
    }
    try
    {
        return checked_multiply(*mebibytes, 1024ULL * 1024ULL);
    }
    catch (const std::overflow_error&)
    {
        return {};
    }
}

SternheimerOptionalValue<std::uint64_t> parse_sternheimer_kib_field(const std::string& text, const std::string& key)
{
    std::size_t first = 0;
    while (first <= text.size())
    {
        const std::size_t end = text.find('\n', first);
        const std::string line = text.substr(first, end == std::string::npos ? text.size() - first : end - first);
        const std::size_t colon = line.find(':');
        if (colon != std::string::npos && trim_copy(line.substr(0, colon)) == key)
        {
            std::string value = trim_copy(line.substr(colon + 1));
            const std::size_t space = value.find_first_of(" \t");
            const std::string number = value.substr(0, space);
            value = space == std::string::npos ? std::string{} : trim_copy(value.substr(space));
            if (value != "kB")
            {
                return {};
            }
            const SternheimerOptionalValue<std::uint64_t> kibibytes = parse_unsigned_decimal(number);
            if (!kibibytes)
            {
                return {};
            }
            try
            {
                return checked_multiply(*kibibytes, 1024ULL);
            }
            catch (const std::overflow_error&)
            {
                return {};
            }
        }
        if (end == std::string::npos)
        {
            break;
        }
        first = end + 1;
    }
    return {};
}

SternheimerOptionalValue<std::string> parse_sternheimer_cgroup_v2_path(const std::string& text)
{
    std::size_t first = 0;
    while (first <= text.size())
    {
        const std::size_t end = text.find('\n', first);
        const std::string line = text.substr(first, end == std::string::npos ? text.size() - first : end - first);
        if (line.compare(0, 3, "0::") == 0 && line.size() > 3)
        {
            return std::string(line.substr(3));
        }
        if (end == std::string::npos)
        {
            break;
        }
        first = end + 1;
    }
    return {};
}

SternheimerOptionalValue<std::string> parse_sternheimer_cgroup_v1_memory_path(const std::string& text)
{
    std::size_t first = 0;
    while (first <= text.size())
    {
        const std::size_t end = text.find('\n', first);
        const std::string line = text.substr(first, end == std::string::npos ? text.size() - first : end - first);
        const std::size_t first_colon = line.find(':');
        const std::size_t second_colon
            = first_colon == std::string::npos ? std::string::npos : line.find(':', first_colon + 1);
        if (second_colon != std::string::npos)
        {
            const std::string controllers = line.substr(first_colon + 1, second_colon - first_colon - 1);
            std::size_t controller_first = 0;
            while (controller_first <= controllers.size())
            {
                const std::size_t comma = controllers.find(',', controller_first);
                const std::string controller = controllers.substr(
                    controller_first,
                    comma == std::string::npos ? controllers.size() - controller_first : comma - controller_first);
                if (controller == "memory" && second_colon + 1 < line.size())
                {
                    return std::string(line.substr(second_colon + 1));
                }
                if (comma == std::string::npos)
                {
                    break;
                }
                controller_first = comma + 1;
            }
        }
        if (end == std::string::npos)
        {
            break;
        }
        first = end + 1;
    }
    return {};
}

SternheimerMemorySnapshot select_sternheimer_memory_snapshot(const SternheimerMemoryCandidates& candidates,
                                                             const int local_mpi_ranks)
{
    if (local_mpi_ranks <= 0)
    {
        throw std::invalid_argument("Sternheimer memory detection requires a positive local MPI rank count.");
    }

    SternheimerOptionalValue<std::uint64_t> cgroup_limit = candidates.cgroup_limit_bytes;
    if (cgroup_limit)
    {
        const bool huge_without_physical = !candidates.physical_memory_bytes && *cgroup_limit >= (1ULL << 60);
        const bool huge_relative_to_physical = candidates.physical_memory_bytes && *candidates.physical_memory_bytes > 0
                                               && *cgroup_limit / *candidates.physical_memory_bytes >= 16;
        if (huge_without_physical || huge_relative_to_physical)
        {
            cgroup_limit.reset();
        }
    }

    SternheimerOptionalValue<std::uint64_t> enforced_limit;
    if (cgroup_limit)
    {
        enforced_limit = cgroup_limit;
    }
    if (candidates.slurm_limit_bytes)
    {
        enforced_limit
            = enforced_limit ? std::min(*enforced_limit, *candidates.slurm_limit_bytes) : *candidates.slurm_limit_bytes;
    }

    SternheimerMemorySnapshot snapshot;
    snapshot.local_mpi_ranks = local_mpi_ranks;
    if (enforced_limit && candidates.cgroup_current_bytes)
    {
        snapshot.mode = SternheimerMemoryAccountingMode::node_aggregate;
        snapshot.limit_bytes = *enforced_limit;
        snapshot.current_bytes = *candidates.cgroup_current_bytes;
        if (cgroup_limit && candidates.slurm_limit_bytes)
        {
            snapshot.source = candidates.cgroup_source + "+slurm";
        }
        else if (!candidates.cgroup_source.empty())
        {
            snapshot.source = candidates.cgroup_source;
        }
        else
        {
            snapshot.source = "cgroup_current+slurm";
        }
        return snapshot;
    }
    if (enforced_limit && candidates.process_rss_bytes)
    {
        snapshot.mode = SternheimerMemoryAccountingMode::per_rank;
        snapshot.limit_bytes = *enforced_limit;
        snapshot.current_bytes = *candidates.process_rss_bytes;
        if (cgroup_limit && candidates.slurm_limit_bytes)
        {
            snapshot.source = candidates.cgroup_source + "+slurm+proc_status";
        }
        else if (cgroup_limit)
        {
            snapshot.source = candidates.cgroup_source + "+proc_status";
        }
        else
        {
            snapshot.source = "slurm+proc_status";
        }
        return snapshot;
    }
    if (candidates.mem_available_bytes)
    {
        snapshot.mode = SternheimerMemoryAccountingMode::available;
        snapshot.limit_bytes = *candidates.mem_available_bytes;
        snapshot.current_bytes = 0;
        snapshot.source = "proc_meminfo";
        return snapshot;
    }
    return snapshot;
}

} // namespace detail

} // namespace ModuleRI

#include "source_base/memory_budget.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

namespace ModuleBase
{

namespace
{

std::size_t saturating_mul(const std::size_t a, const std::size_t b)
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    if (a > std::numeric_limits<std::size_t>::max() / b)
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return a * b;
}

std::size_t saturating_add(const std::size_t a, const std::size_t b)
{
    if (a > std::numeric_limits<std::size_t>::max() - b)
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return a + b;
}

// Returns true when the file was readable.  "max" means no cgroup limit and
// is returned as value 0.
bool read_size_file(const std::string& path, std::size_t& value)
{
    std::ifstream ifs(path);
    if (!ifs)
    {
        return false;
    }
    std::string token;
    if (!(ifs >> token))
    {
        return false;
    }
    if (token == "max")
    {
        value = 0;
        return true;
    }
    try
    {
        const unsigned long long parsed = std::stoull(token);
        value = parsed > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())
                    ? std::numeric_limits<std::size_t>::max()
                    : static_cast<std::size_t>(parsed);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

std::size_t read_inactive_file(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs)
    {
        return 0;
    }
    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream iss(line);
        std::string key;
        std::size_t value = 0;
        if (iss >> key >> value && (key == "inactive_file" || key == "total_inactive_file"))
        {
            return value;
        }
    }
    return 0;
}

std::size_t read_mem_available()
{
    std::ifstream ifs("/proc/meminfo");
    if (!ifs)
    {
        return 0;
    }
    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.compare(0, 13, "MemAvailable:") == 0)
        {
            std::istringstream iss(line.substr(13));
            std::size_t value_kb = 0;
            if (iss >> value_kb)
            {
                return saturating_mul(value_kb, 1024);
            }
        }
    }
    return 0;
}

// Resolve the cgroup directory.  The explicit root is used first, then the
// cgroup v2 path from /proc/self/cgroup.  Nested v1 paths are intentionally
// left to the caller-provided root.
std::string resolve_cgroup_dir(const std::string& root)
{
    {
        std::ifstream v2(root + "/memory.max");
        if (v2)
        {
            return root;
        }
    }
    {
        std::ifstream v1(root + "/memory/memory.limit_in_bytes");
        if (v1)
        {
            return root;
        }
    }

    std::ifstream cgroup("/proc/self/cgroup");
    if (!cgroup)
    {
        return root;
    }
    std::string line;
    while (std::getline(cgroup, line))
    {
        if (line.compare(0, 3, "0::") == 0)
        {
            const std::string path = line.substr(3);
            if (!path.empty() && path[0] == '/')
            {
                const std::string candidate = root + path;
                std::ifstream v2(candidate + "/memory.max");
                if (v2)
                {
                    return candidate;
                }
            }
        }
    }
    return root;
}

} // namespace

MemoryStatus get_memory_status()
{
    return get_memory_status("/sys/fs/cgroup");
}

MemoryStatus get_memory_status(const std::string& cgroup_root)
{
    MemoryStatus status;

    const std::string dir = resolve_cgroup_dir(cgroup_root);
    std::size_t limit = 0;
    std::size_t current = 0;
    std::size_t reclaimable = 0;

    if (read_size_file(dir + "/memory.max", limit))
    {
        read_size_file(dir + "/memory.current", current);
        reclaimable = read_inactive_file(dir + "/memory.stat");
        status.from_cgroup = true;
    }
    else if (read_size_file(dir + "/memory/memory.limit_in_bytes", limit))
    {
        read_size_file(dir + "/memory/memory.usage_in_bytes", current);
        reclaimable = read_inactive_file(dir + "/memory/memory.stat");
        status.from_cgroup = true;
    }

    if (!status.from_cgroup)
    {
        status.available_bytes = read_mem_available();
        return status;
    }

    // cgroup v1 reports a huge value instead of "max" for an unlimited cgroup.
    if (limit > (std::size_t(1) << 62))
    {
        limit = 0;
    }

    status.limit_bytes = limit;
    status.current_bytes = current;
    status.reclaimable_bytes = reclaimable;

    if (limit == 0)
    {
        status.available_bytes = read_mem_available();
    }
    else
    {
        const std::size_t free_bytes = limit > current ? limit - current : 0;
        status.available_bytes = std::min(limit, saturating_add(free_bytes, reclaimable));
    }
    return status;
}

LcaoDenseEstimate estimate_lcao_dense(const std::size_t nrloc,
                                      const std::size_t ncloc,
                                      const std::size_t nbands,
                                      const std::size_t nks,
                                      const int nspin,
                                      const bool gamma_only,
                                      const bool gamma_direct,
                                      const bool use_elpa,
                                      const std::size_t element_size)
{
    LcaoDenseEstimate estimate;
    const std::size_t local_elements = saturating_mul(nrloc, ncloc);
    const std::size_t matrix_bytes = saturating_mul(local_elements, element_size);

    estimate.hk_bytes = matrix_bytes;
    estimate.sk_bytes = matrix_bytes;
    estimate.dmk_bytes = saturating_mul(matrix_bytes, nks);

    std::size_t psi_elements = 0;
    if (gamma_only && gamma_direct)
    {
        // Gamma-only direct solvers keep ncol = ncloc, so Psi is a full
        // local N x N block rather than nbands x nrloc.
        psi_elements = saturating_mul(local_elements, static_cast<std::size_t>(nspin > 0 ? nspin : 1));
    }
    else
    {
        psi_elements = saturating_mul(saturating_mul(nks, nbands), nrloc);
    }
    estimate.psi_bytes = saturating_mul(psi_elements, element_size);

    if (use_elpa)
    {
        // ELPA allocates one additional local work matrix of the same shape.
        estimate.elpa_work_bytes = matrix_bytes;
    }

    estimate.per_rank_total_bytes = saturating_add(
        saturating_add(saturating_add(estimate.hk_bytes, estimate.sk_bytes), estimate.dmk_bytes),
        saturating_add(estimate.psi_bytes, estimate.elpa_work_bytes));
    return estimate;
}

MemoryDecision decide_memory(const std::size_t required_lower_bytes,
                             const std::size_t required_projected_bytes,
                             const std::size_t available_bytes)
{
    if (required_lower_bytes > available_bytes)
    {
        return MemoryDecision::Abort;
    }
    if (required_projected_bytes > available_bytes)
    {
        return MemoryDecision::Warn;
    }
    return MemoryDecision::Ok;
}

} // namespace ModuleBase

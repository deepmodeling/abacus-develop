#include "source_base/memory_budget.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
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

// Returns true when the file was readable and parses to a usable limit.
// "max" means no cgroup limit; it is reported as value 0 with
// unlimited set to true so callers can distinguish it from a valid
// finite limit of 0 bytes (a real, if unusual, cgroup configuration).
bool read_size_file(const std::string& path, std::size_t& value, bool& unlimited)
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
        unlimited = true;
        return true;
    }
    try
    {
        const unsigned long long parsed = std::stoull(token);
        value = parsed > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())
                    ? std::numeric_limits<std::size_t>::max()
                    : static_cast<std::size_t>(parsed);
        unlimited = false;
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
    bool limit_unlimited = false;

    bool dummy_unlimited = false;
    if (read_size_file(dir + "/memory.max", limit, limit_unlimited))
    {
        read_size_file(dir + "/memory.current", current, dummy_unlimited);
        reclaimable = read_inactive_file(dir + "/memory.stat");
        status.from_cgroup = true;
    }
    else if (read_size_file(dir + "/memory/memory.limit_in_bytes", limit, limit_unlimited))
    {
        read_size_file(dir + "/memory/memory.usage_in_bytes", current, dummy_unlimited);
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
        limit_unlimited = true;
    }

    status.limit_bytes = limit;
    status.current_bytes = current;
    status.reclaimable_bytes = reclaimable;

    // An explicit "max" (or the v1 huge sentinel) means no cgroup ceiling, so
    // fall back to host memory.  Any other parsed value -- including a finite
    // 0 -- is a real limit the caller must respect.
    if (limit_unlimited)
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

    const std::size_t hs_bytes = saturating_add(estimate.hk_bytes, estimate.sk_bytes);
    const std::size_t hsdmk_bytes = saturating_add(hs_bytes, estimate.dmk_bytes);
    const std::size_t psi_elpa_bytes = saturating_add(estimate.psi_bytes, estimate.elpa_work_bytes);
    estimate.per_rank_total_bytes = saturating_add(hsdmk_bytes, psi_elpa_bytes);
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

bool is_gamma_direct_ks_solver(const std::string& ks_solver)
{
    return ks_solver == "genelpa" || ks_solver == "elpa" || ks_solver == "lapack"
           || ks_solver == "pexsi" || ks_solver == "cusolver" || ks_solver == "cusolvermp";
}

bool is_elpa_ks_solver(const std::string& ks_solver)
{
    return ks_solver == "genelpa" || ks_solver == "elpa";
}

std::size_t per_rank_memory_budget(const std::size_t available_node_bytes, const std::size_t ranks_per_node)
{
    if (ranks_per_node == 0)
    {
        return available_node_bytes;
    }
    return available_node_bytes / ranks_per_node;
}

std::string format_memory_bytes(const std::size_t bytes)
{
    static const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1024.0 && unit < 5)
    {
        value /= 1024.0;
        ++unit;
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << value << " " << units[unit];
    return oss.str();
}

} // namespace ModuleBase

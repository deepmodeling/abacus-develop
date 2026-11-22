#ifndef MEMORY_BUDGET_H
#define MEMORY_BUDGET_H

#include <cstddef>
#include <string>

namespace ModuleBase
{

/**
 * @brief Memory limit and current usage of the current cgroup / host.
 *
 * The values come from cgroup v2 (memory.max / memory.current /
 * memory.stat) or cgroup v1 (memory.limit_in_bytes / memory.usage_in_bytes
 * / memory.stat), with /proc/meminfo MemAvailable as a fallback.  A zero
 * limit means "no cgroup limit"; available_bytes is then the host value.
 */
struct MemoryStatus
{
    std::size_t limit_bytes = 0;
    std::size_t current_bytes = 0;
    std::size_t reclaimable_bytes = 0;
    std::size_t available_bytes = 0;
    bool from_cgroup = false;
};

/// Read the memory status from the default cgroup root.
MemoryStatus get_memory_status();

/// Read the memory status.  cgroup_root is injectable for unit tests.
MemoryStatus get_memory_status(const std::string& cgroup_root);

/**
 * @brief Per-rank estimate of the dense LCAO matrices used by the SCF path.
 *
 * This is a lower-bound model, not an allocator: it covers H(k), S(k), DMK,
 * Psi and the ELPA work matrix.  Real-space H(R)/S(R)/DMR, gint buffers and
 * charge grids are added by the caller when they are known.
 */
struct LcaoDenseEstimate
{
    std::size_t hk_bytes = 0;
    std::size_t sk_bytes = 0;
    std::size_t dmk_bytes = 0;
    std::size_t psi_bytes = 0;
    std::size_t elpa_work_bytes = 0;
    std::size_t per_rank_total_bytes = 0;
};

/**
 * @brief Estimate the local dense LCAO matrices in bytes.
 *
 * nrloc/ncloc are the local block dimensions.  For gamma-only direct
 * solvers Psi stores a full local N x N block, so gamma_direct must be true;
 * otherwise Psi is nks * nbands * nrloc.  element_size is sizeof(TK).
 */
LcaoDenseEstimate estimate_lcao_dense(std::size_t nrloc,
                                      std::size_t ncloc,
                                      std::size_t nbands,
                                      std::size_t nks,
                                      int nspin,
                                      bool gamma_only,
                                      bool gamma_direct,
                                      bool use_elpa,
                                      std::size_t element_size);

enum class MemoryDecision
{
    Ok,
    Warn,
    Abort
};

/// Abort when the lower bound does not fit; warn when the projection does not.
MemoryDecision decide_memory(std::size_t required_lower_bytes,
                             std::size_t required_projected_bytes,
                             std::size_t available_bytes);

} // namespace ModuleBase

#endif // MEMORY_BUDGET_H

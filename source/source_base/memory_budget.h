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
 * / memory.stat), with /proc/meminfo MemAvailable as a fallback.
 *
 * available_bytes is meaningful only when known is true.  A readable cgroup
 * whose usage has reached its limit with no reclaimable cache reports
 * available_bytes == 0 with known == true; that is a real, enforceable zero
 * and must not be confused with the "could not read anything" state, which
 * sets known == false (and leaves available_bytes at 0).  unlimited is true
 * only when the cgroup reports no ceiling ("max" or the v1 sentinel); a
 * finite limit of 0 is reported as unlimited == false so it stays
 * distinguishable from the unlimited case.
 */
struct MemoryStatus
{
    std::size_t limit_bytes = 0;
    std::size_t current_bytes = 0;
    std::size_t reclaimable_bytes = 0;
    std::size_t available_bytes = 0;
    bool from_cgroup = false;
    bool known = false;
    bool unlimited = false;
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

/// Return true when ks_solver is a gamma-only direct solver that keeps Psi as
/// a full local N x N block (genelpa/elpa/lapack/pexsi/cusolver/cusolvermp).
/// This only checks solver membership; the caller must also require gamma_only,
/// matching the predicate in source_psi/setup_psi.cpp.
bool is_gamma_direct_ks_solver(const std::string& ks_solver);

/// Return true when ks_solver uses ELPA, which allocates an additional local
/// work matrix of the same shape (genelpa/elpa).
bool is_elpa_ks_solver(const std::string& ks_solver);

/// Per-rank memory budget derived from a node-level available value.
/// ranks_per_node == 0 (serial build before init) returns the node value
/// unchanged; otherwise returns available / ranks_per_node.  Division cannot
/// overflow, so this is safe for any input.
std::size_t per_rank_memory_budget(std::size_t available_node_bytes, std::size_t ranks_per_node);

/// Format a byte count as a human-readable string (B/KB/MB/GB/TB/PB) with
/// two decimal places.
std::string format_memory_bytes(std::size_t bytes);

} // namespace ModuleBase

#endif // MEMORY_BUDGET_H

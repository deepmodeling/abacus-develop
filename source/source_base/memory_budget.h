#ifndef MEMORY_BUDGET_H
#define MEMORY_BUDGET_H

#include <cstddef>
#include <string>

namespace ModuleBase
{

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

/// Return true when ks_solver is a gamma-only direct solver that keeps Psi as
/// a full local N x N block (genelpa/elpa/lapack/pexsi/cusolver/cusolvermp).
/// This only checks solver membership; the caller must also require gamma_only,
/// matching the predicate in source_psi/setup_psi.cpp.
bool is_gamma_direct_ks_solver(const std::string& ks_solver);

/// Return true when ks_solver uses ELPA, which allocates an additional local
/// work matrix of the same shape (genelpa/elpa).
bool is_elpa_ks_solver(const std::string& ks_solver);

/// Format a byte count as a human-readable string (B/KB/MB/GB/TB/PB) with
/// two decimal places.
std::string format_memory_bytes(std::size_t bytes);

} // namespace ModuleBase

#endif // MEMORY_BUDGET_H

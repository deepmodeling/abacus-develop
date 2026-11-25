#ifndef LCAO_MEMORY_CHECK_H
#define LCAO_MEMORY_CHECK_H

#include <cstddef>
#include <ostream>

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/input_parameter.h"

namespace LCAO_domain
{

/**
 * @brief Estimate the dense LCAO matrix memory footprint, compare it with the
 *        per-rank available memory, and print a diagnostic block before any
 *        large allocation is performed.
 *
 * Invoked from LCAO_domain::set_psi_occ_dm_chg (once, during
 * ESolver_KS_LCAO::before_all_runners), after init_basis_lcao has configured
 * the 2D block distribution (so nrloc/ncloc are valid) and before
 * Setup_Psi::allocate_psi allocates the first big dense matrix.
 *
 * This function reads no process-global state; every cross-layer value is
 * injected through the parameter list so the body stays free of hidden
 * coupling and is unit-testable.
 *
 * Behavior:
 *  - rank 0 prints the estimate to @p log (the running log), not stdout;
 *  - every rank evaluates the decision; on Abort every rank calls
 *    WARNING_QUIT so no rank is left hanging on a subsequent MPI call;
 *  - when available memory cannot be determined (0 bytes) the check only
 *    prints the estimate and never aborts, avoiding false positives on
 *    non-Linux or cgroup-less environments.
 *
 * @tparam TK k-space type: double for gamma-only, std::complex<double>
 *            otherwise.  Fixes the element size of the dense matrices.
 * @param pv             parallel orbitals, providing local block dimensions
 *                       nrloc/ncloc and the global basis count.
 * @param kv             k-point metadata, providing nks.
 * @param inp            input parameters (ks_solver, nspin, nbands).
 * @param gamma_only     effective gamma-only flag (the runtime value, not the
 *                       raw input), matching what Setup_Psi::allocate_psi sees.
 * @param my_rank        this rank's index in the global communicator.
 * @param ranks_per_node number of MPI ranks sharing one node's memory.
 * @param log            running-log stream to write the estimate to (rank 0).
 */
template <typename TK>
void check_lcao_dense_memory(const Parallel_Orbitals& pv,
                             const K_Vectors& kv,
                             const Input_para& inp,
                             bool gamma_only,
                             int my_rank,
                             int ranks_per_node,
                             std::ostream& log);

} // namespace LCAO_domain

#endif

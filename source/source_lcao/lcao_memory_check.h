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
 * @brief Estimate the total memory of one full set of dense LCAO matrices and
 *        print the estimate to the running log on rank 0 before any large
 *        allocation is performed.
 *
 * Invoked from LCAO_domain::set_psi_occ_dm_chg (once, during
 * ESolver_KS_LCAO::before_all_runners), before Setup_Psi::allocate_psi
 * allocates the first big dense matrix.
 *
 * The estimate uses the global NLOCAL x NLOCAL matrix dimensions, so it
 * reports the memory of a single non-distributed copy of H, S, DMK and Psi
 * (plus the ELPA work matrix), independent of the 2D block-cyclic
 * distribution.  It is purely informational: the function never reads
 * cgroup/proc memory state and never aborts, so it is safe on non-Linux or
 * cgroup-less environments.
 *
 * This function reads no process-global state; every cross-layer value is
 * injected through the parameter list so the body stays free of hidden
 * coupling and is unit-testable.
 *
 * @tparam TK k-space type: double for gamma-only, std::complex<double>
 *            otherwise.  Fixes the element size of the dense matrices.
 * @param pv             parallel orbitals, providing the global basis and band
 *                       counts.
 * @param kv             k-point metadata, providing nks.
 * @param inp            input parameters (ks_solver, nspin).
 * @param gamma_only     effective gamma-only flag (the runtime value, not the
 *                       raw input), matching what Setup_Psi::allocate_psi sees.
 * @param my_rank        this rank's index in the global communicator; only
 *                       rank 0 prints.
 * @param log            running-log stream to write the estimate to (rank 0).
 */
template <typename TK>
void check_lcao_dense_memory(const Parallel_Orbitals& pv,
                             const K_Vectors& kv,
                             const Input_para& inp,
                             bool gamma_only,
                             int my_rank,
                             std::ostream& log);

} // namespace LCAO_domain

#endif

#ifndef STERNHEIMER_GALERKIN_AUDIT_H
#define STERNHEIMER_GALERKIN_AUDIT_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace ModuleRI
{

struct SternheimerGalerkinAuditMatrices
{
    std::size_t dimension = 0;
    // Column-major matrices in the ORIGINAL input-state coordinates, including
    // zero/dependent states. C = J* Phi_f; every operator term is C* H_c,term C.
    std::vector<SternheimerFDHamiltonian::Complex> overlap;
    std::vector<SternheimerFDHamiltonian::Complex> kinetic;
    std::vector<SternheimerFDHamiltonian::Complex> local_potential;
    std::vector<SternheimerFDHamiltonian::Complex> nonlocal;
    std::vector<SternheimerFDHamiltonian::Complex> hamiltonian;
    std::vector<double> fine_norm_squared;
    std::vector<double> projected_norm_squared;
    // ||J C_i - Phi_f,i||_f / ||Phi_f,i||_f; defined as zero for a zero input state.
    std::vector<double> reconstruction_relative_error;
    // Maximum absolute H-matrix discrepancy against independently contracted
    // (J C)* H_f (J C), using direct fine-H applications and fine-grid BLAS.
    // This compares FILTERED fine states, not the unfiltered Phi_f* H_f Phi_f.
    double filtered_matrix_identity_max_abs_error = 0.0;
    std::size_t workspace_bytes = 0;
};

// Allocation-free peak numerical-payload bound, including returned matrices,
// norms, nonidentity caches for C and J C, and blocked/streamed workspace.
// Excludes the borrowed H_f and input states, FFTW/BLAS opaque internal storage,
// allocator/container and runtime overhead. Not a process-RSS limit.
std::size_t sternheimer_galerkin_audit_workspace_bytes(
    const SternheimerFDHamiltonian& fine_hamiltonian,
    const SternheimerFDHamiltonian::Grid& coarse_grid,
    std::size_t dimension);

// Synchronous spectrum-only audit. All views must be nonnull finite full fine-grid
// vectors in H_f's Bloch sector and remain alive/immutable throughout the call.
// Values alone cannot verify KS provenance or source/target-sector assignments.
// No normalization, state selection, energy shifts, diagonalization, or file I/O.
// With equal grids, borrows the fine views without making another full-state copy;
// that call yields the UNFILTERED fine reference. With different grids, caches C
// and J C once each; operator images are streamed in fixed-size state blocks.
// Insufficient budgets are rejected before numerical allocation, never truncated.
// Each call owns its mutable contexts; a shared H_f must remain immutable.
// Neither these matrices nor a positive selected subspace prove full-Pc positivity.
SternheimerGalerkinAuditMatrices audit_sternheimer_galerkin_matrices(
    std::shared_ptr<const SternheimerFDHamiltonian> fine_hamiltonian,
    const std::vector<const SternheimerFDHamiltonian::Vector*>& fine_states,
    SternheimerFDHamiltonian::Grid coarse_grid,
    std::size_t max_workspace_bytes);

} // namespace ModuleRI

#endif

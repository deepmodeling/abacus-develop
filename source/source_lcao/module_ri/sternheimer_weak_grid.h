#ifndef STERNHEIMER_WEAK_GRID_H
#define STERNHEIMER_WEAK_GRID_H

#include "source_lcao/module_ri/sternheimer_delta.h"
#include "source_lcao/module_ri/sternheimer_grid_transfer.h"

#include <memory>

namespace ModuleRI
{

// Prefix-preserving normalization U <- U R^-1, where S=dV U*U=R*R and R
// is upper triangular. Applies the SAME transform to values and all three
// supplied analytic gradients; no state dropping, shifts or full fine-U copy.
// Rejects lambda_min(S) <= min_metric*max(1,lambda_max(S)). The budget covers
// two nstate-square matrices, eigensolver work and one all-state grid tile;
// excludes caller fields, allocator and opaque BLAS/LAPACK runtime storage.
// Invalid inputs, metric or insufficient budget leave states untouched.
// Arithmetic overflow during the in-place transform may leave partial changes.
// Returns the recomputed max absolute entry of dV U*U-I; callers set their
// own acceptance gate on this diagnostic. Each original prefix span is retained.
double orthonormalize_sternheimer_weak_states_in_place(
    std::vector<SternheimerDeltaGridFunction>& states,
    double volume_element,
    double min_metric = 1e-12,
    std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);

struct SternheimerWeakGridBlocks
{
    std::size_t state_count = 0;
    std::size_t coarse_count = 0;
    // Column-major H_U=h(U,U), L=<U,E>, K=h(U,E). U must be fine-metric
    // orthonormal, including identical transformations of analytic gradients.
    SternheimerFDHamiltonian::Vector state_hamiltonian;
    SternheimerFDHamiltonian::Vector state_coarse_overlap;
    SternheimerFDHamiltonian::Vector state_coarse_hamiltonian;
    double maximum_metric_error = 0.0;
};

// Fine-integrated weak form, not the FullGrid-FD8 Galerkin reference.
// E=J/sqrt(dVc) has Euclidean-normalized coefficients: E* Wf E=I.
// Kinetic uses analytic i(G+k), while local/nonlocal terms retain fine H data.
// Each concurrent worker needs its own instance and immutable shared fine H.
class SternheimerWeakGridOperator
{
  public:
    using Hamiltonian = SternheimerFDHamiltonian;
    using Grid = Hamiltonian::Grid;
    using Vector = Hamiltonian::Vector;
    using Function = SternheimerDeltaGridFunction;

    SternheimerWeakGridOperator(std::shared_ptr<const Hamiltonian> fine_hamiltonian,
                               const Grid& coarse_grid,
                               std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);
    SternheimerWeakGridOperator(const SternheimerWeakGridOperator&) = delete;
    SternheimerWeakGridOperator& operator=(const SternheimerWeakGridOperator&) = delete;

    const Grid& grid() const;
    const Hamiltonian& fine_hamiltonian() const;
    double coarse_volume_element() const;
    double fine_volume_element() const;

    // No dense coarse matrix or explicit E is allocated. Input/output may alias.
    void apply(const Vector& coefficients, Vector& result);
    void lift(const Vector& coefficients, Vector& fine_values);
    void project(const Vector& fine_values, Vector& coefficients);

    SternheimerWeakGridBlocks assemble_blocks(const std::vector<Function>& states,
                                             double metric_tolerance = 1e-8,
                                             std::size_t max_matrix_bytes = 512ULL * 1024 * 1024);

    // Owned numerical payload plus peak scalar nonlocal coefficients. Excludes
    // H data, caller states/results, FFTW opaque plans, BLAS/runtime allocations.
    // Block assembly additionally uses the blocked weak-integral routine's
    // temporary tiles; max_matrix_bytes bounds its dense matrix payload only.
    static std::size_t workspace_bytes_required(const Hamiltonian& fine, const Grid& coarse);

  private:
    std::shared_ptr<const Hamiltonian> fine_;
    std::size_t workspace_bytes_;
    SternheimerGridTransfer transfer_;
    double coarse_dv_;
    double fine_dv_;
    Vector fine_input_, fine_output_, coarse_temporary_, coarse_result_;
};

} // namespace ModuleRI
#endif

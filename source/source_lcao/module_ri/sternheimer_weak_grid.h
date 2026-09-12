#ifndef STERNHEIMER_WEAK_GRID_H
#define STERNHEIMER_WEAK_GRID_H

#define ABACUS_STERNHEIMER_WEAK_EXACT_APPLY_CACHE 1

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
    ~SternheimerWeakGridOperator();
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

    // Explicit opt-in, normally after assemble_blocks. Neither cache changes
    // fine blocks, lift/project or analytic kinetic energy. false/false disables.
    // NL retains every projector and D entry. Local retains exactly the fine
    // discrete-quadrature difference frequencies on min(Nf,2*Nc) per axis.
    // Budget bounds EXTRA numerical payload at construction peak, including
    // any old cache until replacement succeeds; excludes the existing operator,
    // shared H, allocator/object overhead and opaque FFTW/BLAS allocations.
    // A failed enable leaves the old cache active. Concurrent calls need separate
    // operator instances; the cache owns mutable FFT and application scratch.
    void enable_exact_apply_cache(bool cache_nonlocal, bool cache_local,
                                 std::size_t max_extra_workspace_bytes = 512ULL * 1024 * 1024);
    bool has_exact_nonlocal_cache() const;
    bool has_exact_local_cache() const;
    std::size_t exact_cache_storage_bytes() const;
    // Fresh-cache peak, not RSS. For replacement add exact_cache_storage_bytes().
    static std::size_t exact_cache_workspace_bytes_required(const Hamiltonian& fine, const Grid& coarse,
                                                           bool cache_nonlocal, bool cache_local);

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
    struct ExactApplyCache;
    std::unique_ptr<ExactApplyCache> exact_cache_;
};

} // namespace ModuleRI
#endif

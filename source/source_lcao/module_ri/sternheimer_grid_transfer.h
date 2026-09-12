#ifndef STERNHEIMER_GRID_TRANSFER_H
#define STERNHEIMER_GRID_TRANSFER_H

#include "source_lcao/module_ri/sternheimer_fd_hamiltonian.h"

#include <array>
#include <cstddef>
#include <memory>

namespace ModuleRI
{

// Transfers full complex Bloch fields on periodic, z-fast grids of the same cell and k.
// Modes use i <= N/2 ? i : i-N: the unique even-grid Nyquist mode is positive.
// No conjugate pairing, Nyquist merging, or wrapping of k to another sector is performed.
// J is isometric for dV = |det(cell)|/N; restrict_adjoint is its weighted adjoint J*.
// This interpolation does not recover modes outside the chosen coarse Fourier space.
// Each worker must own its instance: execution is non-const and is NOT reentrant,
// even with distinct output vectors. A shared Hamiltonian must use worker-owned
// transfer contexts, not one shared mutable transfer. FFT execution takes no lock.
// Const grid accessors and the allocation-free workspace estimate do not mutate state.
// Plan management uses the shared OpenMP FFTW critical region. Without OpenMP,
// callers must serialize construction/destruction across native threads.
class SternheimerGridTransfer
{
  public:
    using Grid = SternheimerFDHamiltonian::Grid;
    using Vector = SternheimerFDHamiltonian::Vector;

    // Fine dimensions must be >= coarse dimensions in every direction. Effective
    // lattice entries and k must agree to relative/absolute 1e-12, respectively.
    // Invalid geometry, nonperiodic grids and point counts exceeding INT_MAX throw.
    SternheimerGridTransfer(const Grid& coarse_grid, const Grid& fine_grid);
    ~SternheimerGridTransfer();
    SternheimerGridTransfer(SternheimerGridTransfer&&) noexcept;
    SternheimerGridTransfer& operator=(SternheimerGridTransfer&&) noexcept;
    SternheimerGridTransfer(const SternheimerGridTransfer&) = delete;
    SternheimerGridTransfer& operator=(const SternheimerGridTransfer&) = delete;

    const Grid& coarse_grid() const;
    const Grid& fine_grid() const;

    // Output is resized; input and output may refer to the same vector.
    // Equal-size grids take an exact copy path, without FFTs or phase roundoff.
    void interpolate(const Vector& coarse, Vector& fine);
    void restrict_adjoint(const Vector& fine, Vector& coarse);

    // D_a = partial_a J uses i(G+k)_a in Cartesian coordinates, with reciprocal
    // lattice vectors including 2*pi. Gradients differentiate the FULL Bloch
    // field, not just its periodic part. Same-grid values are copied exactly,
    // but same-grid derivatives still use Fourier transforms. The positive
    // Nyquist choice makes no q-star/time-reversal covariance claim.
    // Input may alias fine_values or one gradient. fine_values must not alias
    // any fine_gradients entry. All outputs are resized. New spectral operations
    // reject wrong-sized/nonfinite inputs and insufficient budgets before writing
    // outputs. The budget bounds internal numeric payload, excluding outputs.
    // Arithmetic overflow throws overflow_error; outputs may be partial then.
    void interpolate_with_gradients(const Vector& coarse,
                                    Vector& fine_values,
                                    std::array<Vector, 3>& fine_gradients,
                                    std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);

    // Sum_a D_a* f_a with the dV-weighted adjoint: multiplier -i(G+k)_a
    // and fine-grid Fourier normalization 1/Nf. Output may alias any input.
    void restrict_gradient_adjoint(const std::array<Vector, 3>& fine_fields,
                                   Vector& coarse,
                                   std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);

    // Analytic coarse-space -Laplacian, with symbol |G+k|^2 and no kinetic
    // prefactor. Equals sum_a D_a* D_a, NOT the FD kinetic operator. Uses only
    // coarse FFT execution, no fine gradient arrays. Input/output may alias.
    void apply_negative_laplacian(const Vector& coarse,
                                  Vector& result,
                                  std::size_t max_workspace_bytes = 512ULL * 1024 * 1024);

    // Allocation-free validation and numeric payload estimate: two FFT buffers,
    // two Bloch-phase arrays, and one coarse-to-fine index map. Identity needs zero.
    // Excludes output vectors, opaque FFTW plans, allocator and small object overhead;
    // this is not a process-memory bound or a constructor memory-budget guarantee.
    static std::size_t workspace_bytes_required(const Grid& coarse_grid, const Grid& fine_grid);

    // Total internal numeric payload after enabling ANY spectral operation above.
    // Nonidentity reuses the value-transfer workspace with no extra arrays.
    // Identity lazily allocates the same two FFT buffers, phase arrays and index
    // map; those stay owned by the instance. Use this estimate instead of the
    // value-only estimate thereafter. Exclusions are the same as above; four
    // fine output fields are caller-owned and are NOT part of this estimate.
    static std::size_t gradient_workspace_bytes_required(const Grid& coarse_grid, const Grid& fine_grid);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ModuleRI

#endif

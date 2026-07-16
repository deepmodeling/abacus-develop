# Sternheimer Auxiliary-Channel OpenMP Design

## Scope

This change parallelizes the independent auxiliary-basis response equations
owned by one frequency MPI rank. The existing MPI decomposition over imaginary
frequencies remains unchanged: with six frequencies and six ranks, each rank
continues to own one frequency, while its OpenMP team solves different ABFS
columns concurrently.

The change does not alter the Sternheimer or Delta-Sternheimer equations, the
frequency grid, spin accumulation, matrix symmetrization, MPI file ownership,
or the LibRPA reader-v1 format.

## Parallel Unit

For fixed spin, frequency, and occupied band, channel \(\nu\) solves

\[
 Q_\sigma(H_\sigma-\epsilon_{i\sigma}+i\omega)Q_\sigma
 \delta\psi_{i\sigma}^{\nu}
 =-Q_\sigma V_\nu^{\mathrm H}\psi_{i\sigma}.
\]

All inputs are read-only for different \(\nu\). The resulting wavefunction
contributes only to column \(\nu\) of the unsymmetrized branch matrix,

\[
 B_{\mu\nu}^{\sigma}(i\omega)
 =\sum_{i\in\mathrm{occ},\sigma} f_{i\sigma}
 \langle\psi_{i\sigma}|V_\mu^{\mathrm H}|
 \delta\psi_{i\sigma}^{\nu}(i\omega)\rangle.
\]

Distinct OpenMP iterations therefore read shared Hamiltonian, projector,
occupied states, Delta subspace, perturbations, and probe potentials, but write
different matrix columns and different indexed result records.

## Thread-Safety Boundary

The current finite-difference Hamiltonian and nonlocal projector store only
immutable grid, potential, projector, and D-matrix data during `apply`. Their
temporary vectors are local to each call. The standard and Delta response
solvers allocate GMRES vectors, Hessenberg matrices, denominators, coefficients,
and reconstructed wavefunctions per call. These objects may be shared as const
inputs but must never be reused as writable scratch between channels.

Each OpenMP iteration returns one result containing the solver diagnostics and
the explicit equation residual. Exceptions are captured per channel and
re-thrown after the parallel region in increasing channel order. Global
diagnostics, equation counters, progress output, and maximum residuals are
updated only by a serial ordered pass. This preserves deterministic progress
records and avoids concurrent stream or status-file writes.

The response matrix is a `std::vector<std::complex<double>>` in row-major
ordering. Channel \(\nu\) writes elements `row * naux + nu`; no two concurrent
iterations write the same element. Spin and occupied-band accumulation remain
serial outside the channel loop, so additions from different bands and spins
retain their previous order.

## Runtime Policy

The existing `OMP_NUM_THREADS` controls channel parallelism. Production jobs
use one MPI rank per node and all 30 node cores. `MKL_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1` prevent nested threading in the small per-GMRES LAPACK
least-squares calls. OpenMP dynamic scheduling balances channels with different
iteration counts.

## Verification

1. A unit test must prove that all channel results are returned in channel
   order even when worker completion order differs.
2. A unit test must prove that more than one OpenMP worker can execute channel
   tasks and that the first channel-indexed exception is re-thrown after the
   parallel region.
3. Focused Sternheimer tests and the `abacus_3p` executable must build on a
   `df_dcu` normal node.
4. The same small physical input must be run with `OMP_NUM_THREADS=1` and 30.
   Reader-v1 frequency metadata must be identical, response matrices must agree
   to the existing numerical tolerance, convergence flags and equation counts
   must match, and the 30-thread run must show multi-core CPU utilization.
5. The current N2 12 A run remains the serial baseline and is not cancelled.


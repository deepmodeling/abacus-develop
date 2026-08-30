# Sternheimer Multi-RHS Acceleration Design

## Objective

Reduce the remaining Delta-Sternheimer first-order response wall time by
processing a small group of independent auxiliary-channel right-hand sides in
one operator traversal.  The physical Hamiltonian, Delta-ST decomposition,
GMRES tolerance, frequency grid, response definition, and LibRPA interface
remain unchanged.

## Selected Approach

Use a default micro-batch width of two for channels that share one response
k point, occupied state, frequency, Hamiltonian, fixed subspace, and spectral
preconditioner.  Each column retains its own Krylov basis, Hessenberg matrix,
least-squares solve, residual, iteration count, and convergence decision.  The
implementation is therefore a lockstep collection of independent flexible
GMRES solves, not a coupled block-GMRES method.

Columns that converge early leave the active set immediately.  Remaining
columns continue without waiting for the completed column.  At restart
boundaries, only active columns are submitted to the batched callbacks.

## Batched Primitives

`SternheimerFDHamiltonian::apply_batch` accepts a vector of grid vectors and
applies exactly the same FD, local-potential, and nonlocal-projector operator to
each column.  The grid traversal reuses neighbor indices, Bloch phases, stencil
weights, and local-potential values across columns while retaining the scalar
summation order within each column.

`SternheimerSubspaceProjector::project_batch` projects the same immutable fixed
subspace from each column.  Basis vectors and cached norms are reused; each
column's dot product and sequential basis-vector order remain the same as the
scalar path.

`SternheimerFDSpectralPreconditioner::apply_batch` initially invokes the
accepted scalar FFT preconditioner once per column.  This preserves the
validated FFT path and avoids introducing a second numerical change in the
first multi-RHS implementation.  A separate FFTW `plan_many` optimization may
be considered only after the Hamiltonian/projector batching is measured.

The nonlocal projector computes all column coefficients while one projector
block is resident, then applies the same D matrix and reconstruction for every
column.  The accumulation order for each coefficient and output grid point is
unchanged.

## Batched Delta-ST Solve

A batch entry point receives one RHS and one virtual-state perturbation vector
per channel.  It constructs one fixed projector, denominator table, and
spectral preconditioner for the shared equation.  RHS projection, low-rank
Delta correction, Hamiltonian application, final projection, and physical
residual evaluation use batch primitives.  Reconstruction of SOS, Pulay, and
out-of-subspace response components remains column-local and uses the existing
scalar formulas.

The ABACUS producer groups the existing `owned_channels` in stable channel
order.  Complete groups of two use the batch solve; a final short group uses
its actual width.  `ABACUS_STERNHEIMER_CHANNEL_BATCH_WIDTH=1` is retained only
as a rollback and controlled A/B route.  Missing configuration selects width
two, and output provenance records the effective width.

## Parallelism and Memory

One batch worker owns approximately `batch_width` times the scalar Krylov
storage.  The channel-worker plan keeps the accepted outer OpenMP worker budget
instead of dividing it by batch width; batch columns are not independent outer
workers.  A separate per-rank memory estimate limits the effective width and
worker count.  If memory detection permits fewer than one full batch, or fewer
than two channels are locally owned, the scalar route is used.

This design does not add MPI communication, change q/frequency ownership, or
alter response-file ordering.

## Numerical Gates

Development uses test-first red-green-refactor cycles.  Local acceptance
requires:

- scalar and batch Hamiltonian applications agree within `1e-13` for FD2 and
  FD8, Gamma and twisted periodic grids, skew cells, and nonlocal projectors;
- scalar and batch fixed-subspace projections agree within `1e-13`;
- independent batched GMRES reproduces scalar convergence flags, iteration
  counts, residuals, and solutions within `1e-12` on systems with unequal
  convergence histories;
- scalar and batched Delta-ST responses agree within `1e-12`, including SOS,
  Pulay, out, reconstructed wavefunctions, coefficients, and physical
  residuals;
- all directly affected Sternheimer tests pass.

Remote physical calculations run only on df_dcu `/work1` with immutable
artifacts.  The accepted Si one-frequency case is run with batch widths one,
two, and four using identical inputs and layout.  Acceptance requires
`all_converged=yes`, complete response/provenance artifacts, relative response
Frobenius difference no greater than `1e-8`, no residual regression, and
measured response-wall improvement.  Only then is the complete accepted Si
canonical-q workload used for the final energy and resource gate, with total
energy difference no greater than `1e-8 Ha`.

## Integration Decision

The same-artifact Si A/B/C selected width two as the no-configuration
production default.  It reduced total wall and node-hours by 4.13% relative to
width one, with direct response-matrix relative Frobenius difference
`7.59e-13` and 24.8% higher per-rank MaxRSS.  Width four was numerically
equivalent but improved total wall by only 0.99% while increasing MaxRSS by
about 74%, so it is not the default.  Accepted implementation, executable
hashes, A/B differences, wall time, node-hours, and remaining bottlenecks are
recorded in the development test document and production skill.

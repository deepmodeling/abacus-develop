# Sternheimer Operator Acceleration Design

## Objective

Reduce Delta-Sternheimer first-order response wall time without changing the
physical operator, solver tolerance, response matrices, or RPA energy.  The
already validated finite-difference spectral preconditioner becomes the default
solver path.  Two additional exact optimizations remove repeated work from the
finite-difference Hamiltonian and fixed-subspace projection.

## Scope

This change covers three tightly related parts of the existing solver:

1. default-on selection of the finite-difference spectral preconditioner;
2. caching and exact zero-work elimination in `SternheimerFDHamiltonian`;
3. cached norms for repeated projection onto the fixed occupied plus Delta
   virtual subspace.

The change does not introduce block GMRES, Krylov recycling, altered solver
tolerances, approximate operators, reduced precision, or new physical input
parameters.  Those approaches require separate numerical designs after the
exact optimizations are measured.

## Default Preconditioner Semantics

`SternheimerRPA::SolverOptions::use_fd_spectral_preconditioner` defaults to
`true`.  ABACUS production entry points use the same default when
`ABACUS_STERNHEIMER_FD_ST_SPECTRAL_PRECONDITIONER` is absent.

The environment variable remains an explicit rollback and A/B control:

- absent or a recognized true value: enable the spectral preconditioner;
- a recognized false value: disable it;
- any other nonempty value: fail with a clear invalid-value error.

The regularization default changes from the legacy `0.2` to `0.0`, matching the
accepted HF and Si A/B benchmarks.  Output provenance continues to write
`sternheimer_preconditioner fd_spectral` or `none`, so a production result is
auditable even though no explicit enabling variable is required.

## Finite-Difference Hamiltonian Optimization

### Cached immutable data

Construction of `SternheimerFDHamiltonian` computes and stores all quantities
that depend only on the grid, lattice, and finite-difference order:

- finite-difference radius and first- and second-derivative weights;
- dual lattice vectors and the symmetric Laplacian coefficient matrix;
- the diagonal center coefficient;
- the list of mixed-derivative pairs whose coefficient is numerically nonzero;
- whether a Gamma-only phase shortcut is valid.

`apply_grid_terms` consumes these immutable values and performs no repeated
geometry construction.  The cache is private to the Hamiltonian object and is
read-only during OpenMP execution.

### Exact zero-work elimination

Mixed derivative work is performed only for pairs retained in the cached list.
A pair is removed only when its coefficient is exactly zero after construction;
no tolerance-based truncation is allowed in the production operator.  For the
separable high-order periodic route, first-derivative buffers are constructed
only for directions required by at least one retained pair.

For a periodic grid point that does not cross a cell boundary, the Bloch phase
is exactly `1+0i` and no trigonometric function is evaluated.  At Gamma, all
translations also return `1+0i`.  Boundary translations at non-Gamma k points
continue to use the existing `sternheimer_bloch_phase` function.

These changes preserve operation order for every nonzero stencil term.  Tests
cover periodic and nonperiodic FD2, FD4, FD6, and FD8 operators, orthogonal and
skew cells, Gamma and twisted boundary conditions, and Hermiticity.

## Fixed-Subspace Projection Optimization

Repeated Delta-Sternheimer operator and preconditioner calls project vectors out
of one immutable fixed subspace.  A `SternheimerSubspaceProjector` object owns a
reference to the fixed vectors and stores their dot-product norms once at
construction.  Its `project` method retains the current sequential projection
order and computes only the basis-vector-to-target dot product on each use.

Zero-norm vectors retain the current behavior and are skipped.  Size validation
and the missing-dot-callback error remain explicit.  The existing static
`project_out_subspace` entry point remains available for callers that project
only once; the Delta linear-response solve uses the cached projector for its
repeated operator and preconditioner applications.

The fixed-subspace vectors must remain unchanged for the projector lifetime.
The projector is created inside one linear-response solve and captured by value
through a shared immutable object, so worker threads do not share mutable state.

## Validation

Development follows test-first red-green-refactor cycles.  Each optimization is
committed separately after its focused and neighboring tests pass.

Local gates:

- missing preconditioner environment variable selects `fd_spectral`;
- explicit false selects `none`, and invalid text is rejected;
- cached and reference FD applications agree within `1e-13` for all covered
  grid/order/boundary combinations;
- cached and uncached projection agree within `1e-13`, including non-unit and
  zero-norm basis vectors;
- all Sternheimer FD Hamiltonian, preconditioner, RPA, Delta, periodic solver,
  and production-entry tests pass.

Remote df_dcu gates use `/work1` and a single immutable executable artifact.
No physical calculation runs locally.  HF and Si A/B cases reuse the established
inputs, MPI/OpenMP layout, frequency, grid, pseudopotential, NAO/ABFS, and solver
tolerances.  Acceptance requires:

- application success and `all_converged=yes` for every equation;
- identical equation counts and iteration-audit coverage;
- response relative Frobenius difference no greater than `1e-8`;
- response maximum absolute difference recorded;
- final energy difference no greater than `1e-8 Ha` where LibRPA energy is
  produced;
- no regression in maximum reported physical residual.

## Performance Measurement

Incremental benchmarks measure the spectral-preconditioner baseline, the FD
operator optimization, and the projector optimization separately.  Queue time
is excluded.  Reports include ABACUS process wall time, first-order response wall
time when available, MPI ranks, OpenMP threads, nodes, peak memory, node-hours,
equation count, and iteration count.

After both exact optimizations pass, the final end-to-end benchmark uses the
previously validated Si FD8 solid protocol with analytic head/wing enabled and
the complete canonical-q workload.  The optimized run is compared with the
same-input accepted baseline artifact; completed baseline physics is not rerun
solely to recreate timing.  ABACUS producer and LibRPA postprocessing times are
reported separately and together.

## Integration and Documentation

Implementation commits use Codex as author and AroundPeking as committer.  A
validated feature branch is merged into `master_ghj` only after local and remote
precision gates pass.  The development test document records the algorithms,
commit and executable hashes, A/B response differences, wall times, resources,
speedups, and remaining bottlenecks.  The production skill is updated to state
that the spectral preconditioner is default-on and that explicit false is only
for rollback or controlled comparison.

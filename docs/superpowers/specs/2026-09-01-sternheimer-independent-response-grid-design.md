# Independent Sternheimer Response Grid Design

## Goal

Allow an LCAO PBE calculation to retain the real-space integration grid required
for its converged ground-state energy while evaluating the Delta-Sternheimer
response on a separately converged, lower-cutoff grid.  The first production
use is gated by a small periodic Si calculation; the graphene plus water system
is not released until that gate passes.

## Problem Statement

The current periodic Delta-Sternheimer path reuses the PBE `PW_Basis` grid for
all response objects:

- the local potential,
- sampled occupied and virtual LCAO states,
- nonlocal projectors,
- auxiliary-basis perturbation channels, and
- finite-difference Hamiltonian applications.

This coupling is unnecessarily expensive when the PBE density and local
integration need a fine grid but the response matrix is already converged on a
coarser grid.  It also makes a large vacuum direction especially costly because
the number of response-grid values is replicated in several equation and
subspace objects.

`lcao_ecut` does not solve this problem.  It controls two-center LCAO table
accuracy, while `ecutwfc` controls the PBE density/potential grid.  The response
grid therefore needs an explicit input contract rather than overloading either
existing parameter.

## Scope

Add `sternheimer_response_ecutwfc` to the LCAO Delta-Sternheimer route.

- `0.0` keeps the current PBE grid and is the compatibility default.
- A positive value requests an independent response grid using the same ABACUS
  grid rule as `ecutwfc`, namely a density-grid cutoff of four times the input
  cutoff.
- The requested cutoff must not exceed `ecutwfc`.  This first implementation
  performs spectral restriction only; it does not invent missing high-frequency
  PBE information by upsampling.
- The independent grid is supported only for `sternheimer_delta=true`.  Other
  Sternheimer paths continue to use the PBE grid.
- Reader-v1 data format, frequency grids, auxiliary-basis definitions, Coulomb
  treatment, k/q sampling, and solver tolerances are unchanged.

This change does not make a low response grid a physical default.  Every system
class still requires a response-grid convergence test.

## Grid Construction

Create two serial, replicated `PW_Basis` objects on every response rank.  The
first has the exact PBE grid dimensions and exists only to transform the already
gathered full-grid potential without mixing it with the distributed PBE FFT
layout.  The second uses the converged PBE lattice, a cutoff of
`4 * sternheimer_response_ecutwfc`, a complex positive/negative-G
representation, and the normal ABACUS FFT grid selection.  The resulting
dimensions and actual grid spacing are the response provenance; the nominal
cutoff alone is insufficient.

When `sternheimer_response_ecutwfc=0` or the generated dimensions equal the PBE
dimensions, the code must take the existing path.  This preserves the current
floating-point operation order and avoids an unnecessary FFT round trip.

## Fine-To-Coarse Potential Transfer

The PBE local potential is the only response input that exists first on the fine
PBE grid.  Transfer it by a reciprocal-space projection:

1. gather the complete fine-grid real potential as the current k-point MPI path
   already does;
2. transform it with a serial complex-G basis having the exact PBE dimensions;
3. copy coefficients whose integer reciprocal index exists on the response
   grid;
4. set all other response coefficients to zero; and
5. inverse transform on the response grid.

This is the exact orthogonal projection of the represented fine-grid field onto
the common Fourier subspace.  It preserves constants and all retained plane
waves and removes modes outside the response grid without aliasing.

The resulting finite-difference Hamiltonian is nevertheless not algebraically
identical to applying an explicit Galerkin operator `R H_f P`: its kinetic and
nonlocal terms are rebuilt on the coarse grid.  It is therefore a systematically
convergent response discretization, not an exact change of representation.

## Response Objects

Build all remaining response objects directly on the response grid:

- sample occupied and virtual LCAO eigenvectors from their unchanged LCAO
  coefficients;
- sample nonlocal pseudopotential projectors;
- sample auxiliary-basis perturbation channels; and
- construct the finite-difference Hamiltonian and Delta subspace.

The PBE eigenvalues, LCAO coefficients, occupations, overlap-metric conventions,
and k/q mapping remain unchanged.  No coefficients are interpolated from the PBE
grid.

## Interface and Provenance

Status output must record:

- `pbe_ecutwfc_Ry`,
- `sternheimer_response_ecutwfc_Ry`,
- `sternheimer_response_grid_source` as `pbe` or `independent`,
- PBE grid dimensions and point count, and
- response grid dimensions and point count.

Progress output continues to report the current stages.  `full_grid_ready` means
only that the replicated response grid and fine-to-coarse potential are ready;
it is not evidence that channels, the Delta subspace, or any equation succeeded.

## Validation Gates

### Unit and Interface Gates

1. The parser defaults to `0.0`, accepts a positive cutoff not above `ecutwfc`,
   and rejects negative values or upsampling.
2. Explicit same-grid selection reproduces the compatibility path without
   invoking spectral transfer.
3. Spectral restriction preserves a constant field to `1e-12` relative error.
4. A retained low-G plane wave is preserved to `1e-11` relative error.
5. A fine-grid-only high-G plane wave is removed to `1e-11` of its input norm.
6. Mixed low- and high-G fields preserve the low-G coefficient without aliasing.
7. Existing input, finite-difference, Delta, periodic k/q, and reader-v1 tests
   pass unchanged.

### Small Periodic Si Gate

Use the existing two-atom Si primitive-cell Delta-Sternheimer input because it
is fast and exercises periodic `k -> k+q` sampling.  Keep PP, NAO, auxiliary
basis, PBE cutoff, k/q mesh, FD8, fixed frequency, solver tolerances, and binary
identical.  Change only the response cutoff.

Run the compatibility grid and a cost-ordered response-grid scan.  For the
first gate, a bounded one-q, one-frequency, reduced-channel case is sufficient;
the final selected cutoff must then be checked with the complete channel set.
Require:

- all equations converged with finite residuals;
- matching reader-v1 dimensions, q/frequency metadata, and weights;
- same-grid relative Frobenius difference at or below `1e-10`;
- reduced-grid response-matrix relative Frobenius difference at or below
  `1e-6` for the candidate promoted to the large-system pilot;
- no material worsening of maximum equation residual; and
- measured wall time, node-hours, and peak memory.

The `1e-6` matrix gate is an engineering gate for trying the large system, not a
universal physical convergence threshold.  Final RPA energies require their own
cutoff convergence.

### Graphene Transfer Gate

Only after Si passes, run one 2x2 graphene plus water q/frequency/channel pilot
with the selected response grid.  Keep the PBE grid, geometry, PP, NAO, auxiliary
basis, full state space, symmetry, full 2D Ewald, and no-shrink settings fixed.
Require `channels_ready`, `delta_subspace_ready`, at least two finite converged
equation batches, no OOM, and a viable runtime estimate before releasing further
q points.

The 2x2 pilot is feasibility evidence only.  A physical adsorption energy still
requires complete adsorbed/slab q coverage, the common validated H2O reference,
finite reader-v1 data, LibRPA completion, and a supercell-size convergence test.

## Failure Policy

- If same-grid equivalence fails, stop before any remote physics run.
- If spectral transfer tests fail, do not replace it with real-space point
  sampling or trilinear interpolation.
- If Si exceeds the matrix gate, retain the fine response grid and investigate
  the operator decomposition before trying graphene.
- If the 2x2 pilot OOMs before `channels_ready`, preserve stage and memory
  evidence; do not silently reduce channels, virtual states, q points, or
  Coulomb physics.

## Delivery

Commit implementation and tests to `master_ghj` only after unit and same-grid
regressions pass.  Record the exact source commit, executable hash, grid
dimensions, Si numerical comparison, and any graphene pilot job in the existing
AITP TeX note.  Keep code completion, producer completion, numerical acceptance,
and physical acceptance as separate statements.

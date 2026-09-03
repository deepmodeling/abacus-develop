# Explicit Anisotropic Sternheimer Response Grid Design

## Goal

Allow Delta-Sternheimer calculations to select independent real-space response
grid dimensions in each lattice direction without changing the converged PBE
grid.  The first application fixes the graphene plus water and clean-graphene
PBE calculations at 80 Ry and 13 Angstrom vacuum, then converges the response
grid from low to high cost against the RPA adsorption energy.

## Input Interface

Add three integer INPUT keys:

```text
sternheimer_response_nx  0
sternheimer_response_ny  0
sternheimer_response_nz  0
```

All three zero values preserve current behavior.  All three positive values
request an explicit Delta-Sternheimer response grid.  A mixed zero/positive
triplet is rejected because it leaves the effective grid ambiguous.

The explicit dimensions are valid only for the Delta-Sternheimer reader-v1
producer and are mutually exclusive with a positive
`sternheimer_response_ecutwfc`.  Every requested dimension must be no larger
than the corresponding PBE grid dimension and must be large enough for the
selected finite-difference stencil.  The PBE `ecutwfc` and the global
`nx/ny/nz` inputs are not modified.

The existing scalar `sternheimer_response_ecutwfc` remains supported.  It is
useful for isotropic physical spacing and backward compatibility.  The new
dimensions are the preferred interface for slab calculations that require
different in-plane and out-of-plane resolution.

## Grid Construction And Transfer

Extend `make_sternheimer_response_grid` to accept the explicit response-grid
triplet.  When the triplet equals the PBE dimensions, retain the current PBE
path and operation order.  Otherwise create the serial response `PW_Basis`
with the requested dimensions while retaining the PBE lattice and the current
reciprocal-space representation.

The existing scalar-cutoff route retains its spherical plane-wave restriction.
The explicit-dimension route instead uses a complete rectangular complex FFT:
each Cartesian reciprocal component admitted by the requested dimension is
copied independently, and coefficients outside that rectangular box are
discarded.  When an even target dimension merges the positive and negative
Nyquist modes, both fine-grid coefficients are combined.  This prevents a
sparse z dimension from incorrectly deleting retained in-plane modes and
keeps the transfer alias-free.  Occupied states, virtual states, nonlocal
projectors, and auxiliary perturbations are sampled directly from their
unchanged LCAO definitions on the response grid.

An anisotropic response grid is a new response discretization, not an exact
identity transformation of the finite-difference Hamiltonian.  Its accuracy
must therefore be established by convergence of the observable RPA energy.

## Provenance And Failure Rules

The producer status records the requested response dimensions, actual response
dimensions, PBE dimensions, point counts, and grid source (`pbe`, `cutoff`, or
`explicit`).  Invalid combinations fail before allocating response channels.
Reader-v1 metadata, q/frequency indexing, Coulomb treatment, auxiliary basis,
and solver tolerances remain unchanged.

## Verification

Development follows test-first red-green-refactor cycles.

1. Input tests cover the zero default, a valid positive triplet, mixed
   zero/positive rejection, negative values, and conflict with a positive
   scalar response cutoff.
2. Grid tests prove that an explicit triplet creates exactly the requested
   dimensions, rejects dimensions larger than the PBE grid, and reduces to the
   existing path when the dimensions match the PBE grid.
3. Spectral-transfer tests use unequal dimensions in all three directions and
   prove that retained modes are preserved while a mode removed only in the z
   direction is not aliased into the coarse grid.  They also cover a retained
   in-plane mode outside the sphere implied by a sparse z dimension and the
   even-grid Nyquist boundary.
4. The focused response-grid, input, adapter, Delta, periodic-solver, k/q, and
   ABACUS smoke tests must pass before a production binary is accepted.
5. A same-grid physical control must reproduce the scalar response path within
   its established numerical tolerance before anisotropic production tests.

## Graphene Plus Water Convergence Campaign

Use the 2x2 graphene plus water cell and matching clean slab with 13 Angstrom
vacuum, k441 sampling, the fixed 80 Ry PBE grid, TZDP orbitals, complete
occupied and virtual spaces, complete auxiliary space, FD8, strict 2D Coulomb,
and the same executable.  Reuse the independently converged isolated-water
result; do not rerun it for response-grid convergence.

Run the response grids in increasing point-count order.  First screen a small
set with the same representative q/frequency definition and complete response
space.  Refine x/y resolution and z resolution separately so that an apparent
plateau caused by compensating errors cannot pass.  Only candidates with
finite, complete, all-converged reader-v1 output proceed to complete q and
frequency coverage for both adsorbed and clean-slab systems.

The physical acceptance observable is

```text
E_ads^RPA = E_graphene+H2O^RPA - E_graphene^RPA - E_H2O^RPA.
```

A response grid is accepted when its complete RPA adsorption energy differs by
no more than 10 meV from a strictly denser grid and independent one-direction
refinements do not exceed the same 10 meV threshold.  Matrix norms, residuals,
wall time, and memory are diagnostics; a fixed matrix relative-error threshold
is not substituted for the energy gate.

## Scope Boundaries

This change does not lower the PBE grid, alter the 13 Angstrom geometry, rerun
the water reference, truncate virtual states, reduce auxiliary channels, or
change the strict-2D Coulomb method.  Passing a representative q/frequency
screen is a feasibility result, not a converged adsorption energy.  A physical
result is reported only after complete producer, numerical, LibRPA, and
10 meV adsorption-energy gates pass.

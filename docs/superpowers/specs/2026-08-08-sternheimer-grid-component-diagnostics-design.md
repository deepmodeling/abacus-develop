# Sternheimer Grid-Component Diagnostics Design

## Goal

Identify which grid-dependent term causes the observed Ecut dependence in
the N2 and N Delta-Sternheimer RPA correlation energies, and quantify how much
each term changes the final LibRPA energy.  The diagnostic must run together
with the normal response calculation and must not alter the existing total
reader-v1 matrices.

## Scope

The first production diagnosis uses PCA `1e-4` with the established
12-frequency N2/N definition at 40 and 50 Ry.  The smaller 236/118 auxiliary
dimensions reduce cost, while the two cutoffs use genuinely different FFT
grids.  The diagnostic jobs do not wait for the complete PCA `1e-4` scan.  The
implementation remains restricted to the current molecular, Gamma-point LCAO
Delta-ST route and works with `global_equation` MPI.

This change does not implement separate integration and response grids.  It
produces the evidence needed to decide whether that optimization is valid.

## User Interface

Add one Boolean ABACUS input:

```text
sternheimer_grid_diagnostics 0
```

The default is false.  When false, allocations, output names, and numerical
results remain unchanged.  When true, the option requires
`out_sternheimer_librpa=1` and `sternheimer_delta=1`.

## Response Decomposition

For every spin, occupied state, auxiliary perturbation, and imaginary
frequency, retain the existing Delta-ST decomposition

\[
\Delta\psi=\Delta\psi_{\rm SOS}+\Delta\psi_{\rm Pulay}
            +\Delta\psi_Q.
\]

Accumulate three additional response matrices with the same occupation,
complex-conjugate symmetrization, auxiliary ordering, units, and reader-v1
metadata as the existing total matrix:

\[
M=M_{\rm SOS}+M_{\rm Pulay}+M_Q.
\]

The existing total file keeps its current name.  Additional files use:

```text
v1_sternheimer_component_sos_iq_<iq>_ifreq_<ifreq>_rank0.dat
v1_sternheimer_component_pulay_iq_<iq>_ifreq_<ifreq>_rank0.dat
v1_sternheimer_component_qspace_iq_<iq>_ifreq_<ifreq>_rank0.dat
```

All global-equation ranks accumulate their owned equations and use the same
matrix reduction route as the total response.  Rank zero writes the files.
Before output, require

\[
\frac{\lVert M-M_{\rm SOS}-M_{\rm Pulay}-M_Q\rVert_F}
     {\max(\lVert M\rVert_F,10^{-30})}<10^{-10}.
\]

Failure of this invariant aborts the diagnostic run instead of writing
inconsistent matrices.

## Stable-Basis Grid Matrices

The current Delta Hamiltonian is assembled after grid orthonormalization.  Its
matrix elements cannot be compared element by element across grids because the
orthonormalized vectors may rotate.  Diagnostics therefore also assemble the
following matrices in the original, ordered KS virtual-band basis before grid
Gram-Schmidt:

\[
S^A,\qquad T^A,\qquad V_{\rm loc}^A,\qquad V_{\rm nl}^A,
\qquad H^A=T^A+V_{\rm loc}^A+V_{\rm nl}^A.
\]

Also output the occupied-virtual overlap matrix.  This distinguishes grid
quadrature loss of KS orthogonality from changes in the Hamiltonian terms.
Use one versioned text file per physical spin:

```text
STERNHEIMER_DELTA_GRID_MATRICES_spin_<spin>.dat
```

The file records grid dimensions, volume element, spin, occupied and virtual
dimensions, matrix labels, row/column indices, and complex values.  Require
`H = T + Vloc + Vnl` to relative Frobenius tolerance `1e-12` before writing.

## Perturbation Matrix

Write the frequency-independent matrix

\[
B_{a\mu i}=\langle\eta_a|P_{\mu,h}|\psi_i\rangle_h
\]

in the same original KS virtual-band order.  The versioned file

```text
STERNHEIMER_DELTA_PERTURBATION_spin_<spin>.dat
```

contains dimensions, indices, and complex values.  This localizes grid
sensitivity from sampling the auxiliary Hartree potential and integrating its
product with the occupied and virtual states.  The matrix is computed once per
spin and reused by all frequencies rather than recomputed for diagnostic
output.

## Energy Attribution

A standalone postprocessor reads the 40 and 50 Ry component files and builds,
for each component `X`,

\[
M_{40\rightarrow50}^{X}=M_{40}-M_{40}^{X}+M_{50}^{X}.
\]

Each hybrid matrix is passed through the normal LibRPA `sternheimer_rpa` task.
The report contains

\[
\Delta E_c^X=E_c[M_{40\rightarrow50}^{X}]-E_c[M_{40}],
\]

the complete 40-to-50 Ry energy change, and the non-additive remainder.  The
remainder is required because RPA energy is nonlinear in the response matrix;
individual component shifts are diagnostic influences, not an additive energy
partition.

## Decision Rule for Ecut Strategy

- If `M_Q` carries the dominant matrix and energy change, the response grid is
  accuracy-critical and cannot be made substantially coarser without a better
  finite-difference operator.
- If `M_Q` is stable while `S`, `T`, `Vloc`, `Vnl`, `B`, `M_SOS`, or
  `M_Pulay` changes, retain a fine integration grid and decouple a coarser
  response grid.
- If only the additional tight PCA directions vary, choose the integration
  grid from auxiliary radial resolution and keep the response-grid cutoff as a
  separate convergence parameter.

Both absolute component energies `Ec(N2)` and `Ec(N)` must pass independently;
the dissociation contribution is evaluated only after those two checks.

## Validation

1. Unit-test `H = T + Vloc + Vnl` on a small grid with and without a nonlocal
   projector.
2. Unit-test `M = M_SOS + M_Pulay + M_Q` using deterministic complex vectors.
3. Verify diagnostic-off output is byte-identical to the frozen baseline.
4. Verify one-rank and global-equation MPI component files agree numerically.
5. Run a small H2 smoke before the N2/N 40/50 Ry production jobs.
6. Build and test only on the configured remote server; local compilation is
   excluded by the project execution policy.

## Artifacts and Documentation

Archive source commit, executable hashes, exact inputs, component files,
hybrid matrices, LibRPA logs, and a CSV summary.  Add the component-energy and
operator-matrix conclusions to the existing Sternheimer TeX note before any
two-grid implementation is proposed as production behavior.

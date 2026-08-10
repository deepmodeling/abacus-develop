# Efficient Eighth-Order Sternheimer Finite-Difference Design

## Goal

Determine the lowest-cost finite-difference operator that makes the N and N2
Delta-Sternheimer RPA response insensitive to the 40-to-50 Ry grid change at
the 0.1 kcal/mol target.  Sixth order remains the current production setting.
Eighth order is added as a bounded diagnostic and is promoted only if its
measured accuracy gain justifies its node-hour cost.

## Confirmed Starting Point

With the same TZDP-8au orbitals, PCA `1e-4` auxiliary space, 12 GreenX
frequencies, and fixed full-Ewald Coulomb matrix, changing only the Cartesian
Laplacian order reduced the 40-to-50 Ry total-energy sensitivity as follows:

| System | Second order | Fourth order | Sixth order |
|---|---:|---:|---:|
| N | 0.409318 | 0.069273 | 0.016689 |
| N2 | 1.302188 | 0.237951 | 0.071841 |

The corresponding N2 Q-space hybrid sensitivity is 1.292097, 0.311673, and
0.107662 kcal/mol.  Sixth order passes the total-energy gate but the internal
N2 Q-space diagnostic is slightly above it.  SOS is order-independent to about
0.002 kcal/mol, so this design changes only the grid kinetic operator.

## Scope

Implement `sternheimer_fd_order=8` for the existing molecular, Gamma-point,
uniform-grid Sternheimer route.  Keep the global default at order 2 for input
compatibility and select order 6 explicitly for current production jobs.
Do not implement spectral kinetic energy in this change.  A spectral operator
is the next reference only if eighth order does not pass the Q-space gate.

## Eighth-Order Operator

For one Cartesian direction,

\[
 f''_i = \frac{1}{h^2}\left[
 -\frac{f_{i-4}+f_{i+4}}{560}
 +\frac{8(f_{i-3}+f_{i+3})}{315}
 -\frac{f_{i-2}+f_{i+2}}{5}
 +\frac{8(f_{i-1}+f_{i+1})}{5}
 -\frac{205 f_i}{72}
 \right]+O(h^8).
\]

The three-dimensional Laplacian is the sum of the x, y, and z operators.  The
periodic wrapping and nonperiodic out-of-domain behavior remain identical to
the existing second-, fourth-, and sixth-order implementations.  The operator
remains real symmetric before adding the imaginary frequency shift.

## Efficient Kernel Structure

The current inner grid loop tests the selected order and constructs shifted
indices at every point.  Replace that structure with these bounded changes:

1. Select an order-specialized stencil once at `apply()` entry.  Use
   compile-time radii 1, 2, 3, and 4 with fixed coefficient arrays, so the
   innermost loop has no finite-difference-order branches.
2. Precompute positive and negative shifted coordinate indices for each axis
   and offset in the Hamiltonian constructor.  Storage is
   `2 * radius * (nx + ny + nz)` integers, not another full-grid field.
3. Keep one fused x/y/z accumulation pass and the existing OpenMP
   `collapse(2)` decomposition.  Do not allocate a Laplacian work vector.
4. Preserve the second-order accumulation order so the existing frozen H2
   baseline remains unchanged to its current numerical tolerance.
5. Keep the nonlocal projector call after the local kinetic-plus-potential
   kernel.  No MPI distribution, reader-v1 output, or Delta decomposition is
   changed.

In the current `global_equation` layout each MPI rank owns complete equation
grids, so the radius-four stencil adds local memory traffic but no spatial halo
exchange.  The stencil alone reads 25 points instead of 19 for sixth order;
the expected total wall-time increase is much smaller because projection,
orthogonalization, nonlocal, and response assembly costs are unchanged.

## Input and Compatibility

Accept `sternheimer_fd_order` values 2, 4, 6, and 8.  Reject all other values.
The default remains 2.  Status output records the selected order exactly as it
does now.  Existing order-2/4/6 inputs and output formats remain compatible.

## Validation Gates

### Unit and Interface Gates

1. The input parser accepts 8 and rejects 10.
2. A periodic plane wave has a smaller kinetic-energy error at order 8 than at
   order 6; the order-8 error must be below 10% of the order-6 error for the
   existing 32-point, mode-3 test.
3. Small periodic and nonperiodic dense matrices remain Hermitian.
4. The existing second-, fourth-, sixth-order, solver, Delta, diagnostics, and
   input tests all pass on `df_iopcas_ghj`.

### Real H2 Gate

Run orders 2, 6, and 8 with the existing one-frequency H2 smoke definition.
Require:

- order 2 reproduces the frozen response to relative error below `1e-11`;
- each order reconstructs `SOS + Pulay + Q = total` below `1e-10`;
- order 8 differs finitely from order 6;
- order-8 wall time is recorded but is not rejected from this small noisy case.

### Cost-Ordered Production Gate

1. Run N 40/50 Ry at order 8 first.  If its Q-space sensitivity does not
   decrease relative to order 6, stop without submitting N2 and investigate
   the remaining potential or projection terms.
2. If N passes, run N2 40/50 Ry at order 8 with the same 18-rank by 40-thread
   layout and fixed full-Ewald postprocessing.
3. Require separate N and N2 total changes below 0.1 kcal/mol.  The target for
   the N2 Q-space hybrid is also below 0.1 kcal/mol.
4. Report node-hours for each order and cutoff.  Compare accuracy per node-hour,
   not stencil point count alone.

## Production Decision

- Keep order 6 for production if order 8 changes the N2 total result by less
  than 0.02 kcal/mol and its Q-space improvement is not operationally useful.
- Recommend order 8 for the molecular Delta-ST route if it passes the Q-space
  gate with no more than 20% wall-time overhead relative to order 6 at matched
  cutoff.
- If order 8 remains above the Q-space gate, do not add still higher orders
  blindly.  Implement a small-system spectral-kinetic reference to separate
  residual finite-difference error from local-potential, nonlocal-projector,
  and Q-projection quadrature errors.

## Documentation and Provenance

Record the source commit, remote build job, binary hash, H2 gate, N/N2 job IDs,
wall times, node-hours, total energies, and hybrid component energies in the
existing Sternheimer TeX note.  Preserve all failed and successful remote
artifacts.  ABACUS compilation and physics runs remain server-only.

# Isolated molecular Delta-Sternheimer Coulomb coupling

This opt-in path is for an isolated molecule or atom, not a periodic solid
sampled at Gamma. The global default remains `sternheimer_molecular_coulomb none`.
Ordinary RI-SOS and the periodic Sternheimer Poisson path are not changed.

## Inputs and provenance

Use a serial producer with the same PP, NAO, finalized ABFS, geometry and
numerical settings as the response. Its complete free-space Center2 Coulomb
matrix is selected by `exx_singularity_correction limits`; do not substitute
a periodic Ewald, direct-reciprocal, truncated, or multi-rank local matrix.
Export the AO tensor with `tools/molecular_sternheimer` and archive the input
hashes. The tensor header carries dimensions and Hartree units, not a proof of
radial-basis identity. The workflow must verify that identity separately.

Response INPUT:

```text
symmetry                         -1
exx_singularity_correction        limits
out_sternheimer_librpa            1
sternheimer_delta                 1
sternheimer_q_index               0
sternheimer_molecular_coulomb     isolated_ri
sternheimer_ao_potential_file     ao_potentials.dat
sternheimer_delta_virtual_source ks_bands
sternheimer_fd_order              8
rpa_abfs_preorth                  onsite_coulomb
rpa_abfs_preorth_threshold        1e-2
```

Use an actual single Gamma point, full LCAO bands, a shared molecular/atomic
frequency file, and zero-regularization FD spectral preconditioning. The
experimental AO/tail environment switches are rejected; they cannot silently
activate only half of the new contract. Rotated, moment-split, shrunken, mixed
external-plus-product ABFS, SIAB, forces and periodic q output are rejected
until their shared finalized-basis contract has been validated. External-only
ABFS require `exx_pca_threshold >= 1`; a product-only PCA basis has no
`ABFS_ORBITAL` entries. Runtime guards do not replace scientific convergence.

## What changes

With the paper's notation, the physical kernel is unchanged:

\[
V_\nu(\mathbf r)=\int P_\nu(\mathbf r')/|\mathbf r-\mathbf r'|\,d\mathbf r',
\qquad V_{\mu\nu}=\langle P_\mu|V_\nu\rangle.
\]

The finite AO-space perturbation is evaluated in the producer RI
representation, rather than by a second, independently sampled grid integral:

\[
\langle\phi_p|V_\nu|\phi_q\rangle
\simeq\sum_\mu C_{pq}^{\mu} V_{\mu\nu},
\qquad V_{ai}^{\nu}=\langle\eta_a|V_\nu|\psi_i\rangle.
\]

AO coordinates are recovered for the actual sampled and normalized states
using their sampled AO Gram matrix. Both the right source and the conjugated
left response contraction use these vertices. Hartree vertices are multiplied
by two on entry to the Rydberg linear equation; the response contraction uses
Hartree. A unit test compares a complex multi-state, multi-RHS solve with an
independently assembled full block matrix at both frequency signs.

The sampled AO Gram matrix must have an eigenvalue ratio of at least 1e-8.
Otherwise the calculation stops: a tiny wavefunction reconstruction residual
does not bound the AO coefficient error along a nearly dependent direction.
The ratio is reported as `analytic_ao_gram_reciprocal_condition`. This is an
admission guard, not an auxiliary Coulomb cutoff or silent AO truncation.

The same paper equations are solved again with the new vertices:

\[
c_{ai}^{\nu}=\frac{V_{ai}^{\nu}+\langle D_a|x_i^\nu\rangle}{d_{ai}},
\quad
\left[Q_\Delta L_i Q_\Delta+\sum_a\frac{|D_a\rangle\langle D_a|}{d_{ai}}\right]
|x_i^\nu\rangle
=Q_\Delta|b_i^\nu\rangle-\sum_a\frac{V_{ai}^{\nu}}{d_{ai}}|D_a\rangle.
\]

Here `Q_Delta b` still uses the sampled grid potential. Replacing an SOS term
after solving the old grid equation is not this algorithm. The reported
residual certifies this mixed AO/grid discrete problem, not an independently
assembled pure-grid equation or the continuum limit. RI, grid, box and
partition errors still require convergence tests.

## Exterior potential

For the radial component of the auxiliary charge, the isolated Poisson
solution outside its support is

\[
V_l(r>R_c)=\frac{4\pi}{2l+1}\frac{1}{r^{l+1}}
\int_0^{R_c}s^{l+2}P_l(s)\,ds.
\]

The previous grid sampler set a tabulated potential to zero outside the
table, even when its auxiliary charge had a nonzero multipole moment. The
new sampler retains interior values and continues the exterior using the
actual finalized charge's Simpson-integrated moment. An earlier diagnostic
instead used the last potential value times `R^(l+1)`; that can attach a
small false tail to a nominally zero-moment channel and is not the production
coefficient. No small moment is rounded to zero. `STERNHEIMER_COULOMB_TAIL.dat`
records both coefficients for auditing. The potential table must enclose the
charge support. There is no periodic image sum for this perturbation; the
existing FD Hamiltonian boundary condition is unchanged.

Compact charge is not zero multipole moment. A radial rotation that isolates
one moment-carrying channel does not remove that channel; subsequent mixing
can redistribute its moment. Never impose zero tails on all channels merely
because an input charge is spatially compact.

## Acceptance and historical data

Use the complete producer V with the new response B. Apply the same on-site
congruence scaling to both, with a primary hard threshold of 1e-4; archive
1e-5, 1e-6 and zero as diagnostics. Check exact frequency/occupation/basis
metadata, equation counts, all-converged status, residuals and finite output.
Compare absolute atomic and molecular Ec as well as their dissociation
combination. A finite-NAO SOS value need not equal Delta-ST.

Keep historical Ewald, truncated-tail and endpoint-tail data immutable and
label them by definition. A successful corrected molecular test neither
invalidates every periodic-solid result nor validates a solid's atomic
reference. FHI-aims' diatomic iterative Poisson route supports the same
isolated-kernel requirement but does not use this AO/grid partition.

## Promotion regression (2026-09-07)

Remote p1 build 3267969 compiled the full executable and passed seven focused
CTest targets: input reader, runtime options, RPA, Delta, grid diagnostics,
Coulomb tail and ABACUS ST smoke. Binary SHA256:
`e06c37bfd706bb4a835c826679eda6b5674455c6add75a646641943332ce15b7`.
The sampled-Gram rejection was first observed to fail in red test 3267914.
Seventeen Python exporter tests passed; a real Li2 producer export matched the
earlier tensor to relative Frobenius 3.12619e-15. This is focused validation,
not a full-suite claim: the pre-existing separate ABFS perturbation test
references a missing stream-transform API and was left unchanged.

Before the additional input/Gram admission guards, the same solve arithmetic
was tested with full32 H2/H (3267688/3267690) and single-frequency Li2/GTH-N2
(3267692/3267694). Complete producer V, AO order, frequencies and both-sided
congruence scaling were fixed. Maximum B relative Frobenius differences from
the endpoint-tail controls were 1.65e-12, 6.76e-12, 3.11e-9 and 2.05e-9;
all four archived cutoff energy differences were below 1e-7 Ha. These tests
establish bounded implementation regression, not molecular basis convergence.

Final-artifact H2 pilot 3268100 also completed with 164 converged equations,
maximum residual 8.74030e-7 and sampled AO Gram eigenvalue ratio 3.87252e-4.
Its single-frequency energy changed by less than 5e-17 Ha relative to the
pre-guard full-frequency run's matching point at all four cutoffs.

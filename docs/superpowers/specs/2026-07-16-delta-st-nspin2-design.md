# Delta-Sternheimer Collinear `nspin=2` Design

## Scope

This change extends the molecular Gamma-point LCAO Delta-Sternheimer response
path to collinear `nspin=2` calculations in which both spin channels contain
occupied states. The first production target is the quartet N atom used as the
atomic reference for the N2 binding energy.

The change does not add multi-k-point response, noncollinear spin, `nspin=4`,
spin-orbit coupling, or spin-flip response.

## Reference Behavior

The supplied OpenPFEM Delta-Sternheimer implementation reads one `i_spin` per
execution and uses spin-specific occupied orbitals, eigenvalues, and effective
potential. Separate spin responses are combined only at the polarizability
level. The general FHI-aims independent-particle response follows the same
physical decomposition by summing independent same-spin contributions.

ABACUS will preserve this decomposition inside one run:

\[
 Q_\sigma(H_\sigma-\epsilon_{i\sigma}+i\omega)Q_\sigma
 \delta\psi_{i\sigma}^{\nu}
 =-Q_\sigma V_\nu^{\mathrm H}\psi_{i\sigma},
\]

with

\[
 Q_\sigma=1-\sum_{j\in\mathrm{occ},\sigma}
 |\psi_{j\sigma}\rangle\langle\psi_{j\sigma}|.
\]

There are no cross-spin matrix elements in the collinear density response. The
matrix written for LibRPA is

\[
 M(i\omega)=M^\uparrow(i\omega)+M^\downarrow(i\omega).
\]

## Data Flow

The following data are shared by both spin channels and are built once:

- the uniform real-space grid and volume element;
- auxiliary-basis Hartree perturbations and v1 channel ordering;
- sampled LCAO values and analytic gradients;
- the requested or GreenX minimax imaginary-frequency grid.

Each occupied spin channel is then processed independently:

1. Read its LCAO coefficients, eigenvalues, and occupations.
2. Construct the spin-specific finite-difference Hamiltonian from
   `Potential::get_eff_v(spin_index)` and the shared nonlocal projectors.
3. Construct and orthonormalize only that spin channel's occupied projector.
4. Build a separate reference-value-gradient Delta subspace using the same AO
   candidates but the spin-specific projector and Hamiltonian.
5. Solve all occupied-band and auxiliary-channel response equations owned by
   the current frequency MPI rank.
6. Add the result to the rank-local frequency accumulator.

The spin loop is outside the frequency loop. Once one spin channel has been
accumulated, its Hamiltonian, occupied states, and Delta subspace can be
released before constructing the next channel. This avoids holding two large
Delta subspaces simultaneously. The accumulated matrices are small
`naux x naux` objects, one for each frequency owned by the MPI rank.

After all spin channels are complete, each owned matrix is symmetrized once and
written through the existing LibRPA v1 writer. The v1 interface and LibRPA do
not require a spin dimension.

## Occupation and Conjugate Factors

ABACUS `elec_state.wg` already contains the physical occupation of each state:
normally two for a closed-shell `nspin=1` spatial orbital and one for an
occupied state in each `nspin=2` channel. The existing response accumulation
continues to multiply by `wg`. No additional `2/nspin` factor is applied.

The existing positive-frequency branch and conjugate symmetrization continue
to supply the `2 Re` contribution. Spin summation occurs before this final
symmetrization, and the symmetrization is applied exactly once.

## Frequency Window

When GreenX minimax frequencies are generated internally, the transition
window is the union of all occupied spin-channel windows:

\[
 E_{\min}=\min_\sigma E_{\min}^{\sigma},\qquad
 E_{\max}=\max_\sigma E_{\max}^{\sigma}.
\]

An explicitly supplied frequency-grid file remains authoritative and is not
modified.

## Diagnostics and Errors

The status output records the number and indices of processed spin channels,
occupied bands per spin, and total occupied bands. Existing one-channel
`nspin=1` and spin-polarized H output remains valid.

The code rejects duplicate or out-of-range spin indices, empty channel lists,
inconsistent basis dimensions, non-Gamma input, and spin modes other than
`nspin=1` or collinear `nspin=2`.

For a five-valence-electron N pseudopotential calculation with `nspin=2` and
`nupdown=3`, the expected occupied-band counts are four spin-up and one
spin-down. The expected number of response equations is

\[
 N_{\mathrm{eq}}=N_\omega N_{\mathrm{ABFS}}(4+1).
\]

## Verification

Implementation follows test-driven development:

1. Add a unit test that requires two valid occupied spin channels to be
   accepted and represented independently.
2. Add tests for unioning spin transition windows and additive spin response
   matrices.
3. Run the focused Sternheimer unit tests and full ABACUS build on the remote
   server.
4. Re-run the one-channel spin-polarized H smoke case as a regression.
5. Run a quartet N atom Gamma-point `nspin=2`, `nupdown=3` smoke case and verify
   spin-resolved diagnostics, equation count, convergence, and six v1 files.
6. Compare serial and frequency-MPI v1 matrices element by element.
7. Compare the N atom Delta-Sternheimer response and RPA correlation energy
   against a same-input SOS calculation before starting the N2 campaign.

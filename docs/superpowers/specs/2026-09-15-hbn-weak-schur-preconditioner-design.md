# Weak Schur Spectral Preconditioner Design

## Goal

Reduce the wall time of the fine-integral/coarse-response periodic Delta-Sternheimer weak-form solver without changing its physical operator, response-grid definition, or reader-v1 output contract. Production solves may use a `1e-6` original-equation residual only after a fixed-input comparison quantifies the response and RPA-energy change relative to `1e-8`.

## Existing Equation

The worker solves the normalized Schur equation

\[
  T Q(z) T p = T b, \qquad x = T p, \qquad T=S_W^{-1/2},
\]

where `Q(z)` contains the fine-grid weak Hamiltonian, the complement metric, and the exact virtual-state Schur correction. The current GMRES has no preconditioner. The earlier finite-difference spectral preconditioner cannot be attached directly because it represents an FD stencil on the regular-grid vector, while the repaired weak path uses an analytic spectral coarse kinetic operator and normalized complement coordinates.

## Chosen Design

Add an analytic coarse-grid kinetic preconditioner whose Fourier denominator is

\[
  D(G)=\alpha |G+k|^2-\epsilon+i\omega+r,
\]

using the same lattice, Bloch label, kinetic prefactor, and exact spectral wave vector as the weak-grid kinetic term. `r` is a non-negative numerical regularization and remains zero for admitted production inputs.

The analytic coarse-grid preconditioner is a physical approximate inverse of the shifted operator and acts directly on the regular-grid vector represented by the normalized Schur coordinate `p`. Therefore the normalized GMRES right preconditioner is

\[
  M_p^{-1}=D^{-1}.
\]

It must not be wrapped in additional `S_W^(1/2)` factors. The first hBN+H2O A/B test exposed why: `min eig(S_W)=3.0814849782956344e-08`, so two extra square-root factors attenuated the corresponding direction by approximately `3e-8` and made the right-preconditioned system effectively singular. The `S_W^(-1/2)` construction itself was independently checked against a non-diagonal dense reference; a remote RED test then isolated the callback-coordinate mismatch. The preconditioner changes only the Krylov basis. Final convergence remains determined from the recomputed original augmented residual.

## Runtime Contract

- Weak q production defaults to `spectral` preconditioning.
- `ABACUS_STERNHEIMER_WEAK_PRECONDITIONER=none` disables it for regression comparisons.
- `ABACUS_STERNHEIMER_WEAK_RESIDUAL_TOL` selects the requested original-equation residual and defaults to `1e-6` for weak q production.
- Invalid modes, non-finite tolerances, dimension mismatches, non-finite vectors, and singular Fourier denominators fail explicitly.
- Audit output records the preconditioner, regularization, requested original residual, achieved Schur residual, achieved augmented residual, and iteration count.

## Verification Gates

1. Unit tests prove the complement inverse-square-root map, direct normalized-coordinate callback, output validation, analytic Fourier symbol, Bloch phase, and signed-frequency behavior.
2. Remote ABACUS build on `df_iopcas_ghj` passes the affected tests and existing Sternheimer regression set. No native ABACUS or LibRPA compilation occurs locally.
3. At `1e-8`, preconditioned and unpreconditioned fixed-input hBN+H2O responses must have matching metadata, finite values, and relative Frobenius difference at most `1e-8`; all original augmented residuals must pass `1e-8`.
4. At `1e-6`, the same fixed-input response is compared with the `1e-8` result. All original residuals must pass `1e-6`, and the measured single-q/frequency weighted RPA contribution change must be at most `1 meV`.
5. The preconditioned run must reduce both mean iteration count and wall time. Otherwise the feature remains available but is not used for the remaining production q points.

## Production Release

The currently running q001 recovery is not modified or repeated. Once the gates pass, the implementation is committed with the requested attribution, integrated into `master_ghj`, built remotely with recorded provenance and SHA256, and used for the remaining representative q points. The stricter q001 result is compatible with later `1e-6` solves because every accepted response satisfies at least the production residual gate. Final adsorption-energy acceptance still requires complete q/frequency coverage, reader-v1 validation, LibRPA completion, and a verified subtraction.

# RPA On-Site Coulomb Preorthogonalization Design

## Purpose

Large product auxiliary bases can contain nearly dependent radial functions in
the isolated-atom Coulomb metric.  For molecular full-Ewald RPA calculations,
those directions make the global Coulomb eigenspectrum sensitive to the final
square-root cutoff.  The N2 validation campaign showed that an FHI-aims-style
on-site Coulomb preorthogonalization removes the discontinuous p5-to-p6 failure
when it is combined with a fixed global hard cutoff.

This change moves the validated basis transformation into ABACUS so the RPA
producer and Sternheimer response paths use exactly the same transformed basis.
The legacy path remains available and remains the program default.

## Scope

The new behavior applies only to:

- ABACUS RPA dataset production, including reader-v1 auxiliary-basis and
  Coulomb output; and
- Sternheimer response output for LibRPA.

Ordinary EXX, GW, SCF, force, stress, and unrelated auxiliary-basis paths are
unchanged.  If an RPA producer also reports RI-EXX as part of that RPA dataset,
the reported value uses the same transformed RPA auxiliary basis.

The implementation is general over atom types and angular-momentum channels,
but its accepted production profile is initially limited to isolated molecular
Gamma-point Delta-Sternheimer calculations.  Periodic-solid use remains an
explicit experimental opt-in until a separate solid regression is accepted.

## Runtime Interface

Add two ABACUS input parameters:

```text
rpa_abfs_preorth              none
rpa_abfs_preorth_threshold    1e-2
```

`rpa_abfs_preorth` accepts:

- `none`: preserve the current auxiliary basis exactly;
- `onsite_coulomb`: apply the transformation described below.

The global ABACUS default is `none` for backward compatibility.  Molecular
Delta-Sternheimer production guidance will explicitly select
`onsite_coulomb`, use `rpa_abfs_preorth_threshold 1e-2`, and set LibRPA
`sqrt_coulomb_threshold = 1e-4`.

ABACUS does not infer whether a cell is a molecule from its vacuum or k-point
layout.  Automatic molecule detection would be ambiguous and is intentionally
out of scope.

## Transformation

After product-PCA generation or explicit ABFS loading, process each atom type
and angular-momentum radial channel independently.  For normalized real
spherical harmonics, construct

\[
 V^{(l)}_{nn'} = \frac{4\pi}{2l+1}
 \int dr\,r^2 f_n(r)\int dr'\,{r'}^2 f_{n'}(r')
 \frac{r_<^l}{r_>^{l+1}}.
\]

Traverse radial functions in their existing order.  Normalize each candidate
in this metric, project it against all retained functions, and compute the
remaining squared Coulomb norm.  Reject the candidate when that norm is less
than or equal to `rpa_abfs_preorth_threshold`.  Otherwise perform repeated
modified Gram-Schmidt until the accepted vectors are Coulomb orthonormal.

Retention and rejection operate on complete `(2*l+1)` multiplets because the
decision is made once per radial function.  A non-positive input or residual
Coulomb norm is a fatal input error; the code must not silently repair it.

The transformed in-memory radial functions retain unit Coulomb norm.  They are
not L2-renormalized after the transformation.  Consequently, all downstream
RPA coefficients, Coulomb matrices, and Sternheimer response matrices are
already expressed in the same preorthogonalized basis.  The external diagonal
congruence scaling used by the N2 prototype is therefore unnecessary for data
produced by this runtime mode.

## Code Boundaries

Create a shared module under `source/source_lcao/module_ri/` that:

- computes one radial-channel isolated Coulomb metric;
- applies repeated modified Gram-Schmidt and threshold rejection;
- rebuilds the `Numerical_Orbital_Lm` channel without L2 renormalization; and
- returns per-type/per-channel provenance and numerical diagnostics.

Both the RPA producer setup and `build_sternheimer_abfs` call this shared module
after their existing ABFS construction or loading step.  No duplicate numerical
implementation is allowed between the two paths.

The producer running log and Sternheimer status/provenance output record:

- selected mode and threshold;
- input and output radial counts for every atom type and `l`;
- rejected input indices and minimum residual norm; and
- maximum on-site identity error after transformation.

`basis_aux_out` continues to describe the basis actually used.  Existing hash
and reader-v1 dimension checks remain authoritative.

## Hard Cutoff Contract

The ABACUS switch controls only on-site basis preorthogonalization.  The final
global Coulomb eigenspace cutoff remains a LibRPA runtime parameter.  For the
validated molecular profile it is fixed to:

```text
sqrt_coulomb_threshold = 1e-4
```

The code must not silently rewrite a LibRPA input or choose a cutoff from the
reference energy.  Thresholds `1e-5`, `1e-6`, and zero remain diagnostic
archives, not alternate production defaults.

## Failure Handling

The new mode stops with a clear error if:

- the threshold is not finite or is outside `(0, 1)`;
- a channel has a non-positive Coulomb norm;
- reorthogonalization does not reach its tolerance;
- producer and response transformed dimensions disagree; or
- the transformed metric fails its identity check.

The `none` path performs none of these new operations and must retain legacy
behavior.

## Verification

Development follows test-first implementation.

1. Input tests prove the default is `none`, both modes parse, and invalid modes
   or thresholds fail with a useful message.
2. Unit tests use synthetic radial channels to prove that a nearly duplicate
   radial is rejected, an independent radial is retained, complete multiplet
   counts follow radial decisions, and the retained Coulomb metric is identity.
3. A legacy regression proves `none` leaves radial values and dimensions
   unchanged.
4. Producer and Sternheimer tests prove both paths call the same transform and
   emit matching basis dimensions and provenance.
5. The existing focused module-RI test suite and build must pass.
6. A molecular N2 p5 lowest-frequency regression uses the validated
   `fd_spectral`, zero-regularization executable.  Against the archived external
   preorthogonalization reference, reader-v1 metadata must match and the
   transformed Coulomb and response matrices must each have relative Frobenius
   error no larger than `1e-6`.  Solver residual remains no larger than `1e-6`.

Only after these gates pass is the implementation pushed to `master_ghj`.

## Molecular Skill Update

Update `abacus-delta-st-molecular-convergence` and its parameter contract so
future new molecular campaigns:

- use a `master_ghj` binary containing this feature;
- explicitly set `rpa_abfs_preorth onsite_coulomb` and threshold `1e-2` in both
  producer and response inputs;
- use `fd_spectral` with zero regularization;
- use LibRPA hard threshold `1e-4` as the primary definition;
- require byte-identical producer/response `basis_aux_out` and consistent full
  Coulomb provenance;
- retain p5-to-p6 checks for absolute molecular and atomic correlation energies
  as well as the binding contribution; and
- use `none` only for historical reproduction or an explicitly requested A/B
  diagnostic.

The skill must state that the method is not intrinsically molecule-only, but
that automatic use is limited to molecules because that is the validated
production domain.

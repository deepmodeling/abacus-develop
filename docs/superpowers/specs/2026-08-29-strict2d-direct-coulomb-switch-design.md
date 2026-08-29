# Strict-2D direct Coulomb reader-v1 design

## Scope

Add an explicit, default-off ABACUS input route that writes the strict-2D full auxiliary Coulomb matrices for LibRPA reader-v1 by direct mixed-Fourier quadrature.  The existing 2D Ewald split remains the default and continues to generate the producer exchange and RI coefficients.  Cut Coulomb is not a production alternative.

The direct route is intended for cases where the legacy full-2D-Ewald auxiliary metric is strongly indefinite or cutoff-sensitive even after Hermiticity, radial-mesh, tail, and auxiliary-basis checks.  The accepted WS2 case is the motivating regression: the direct construction changed the full metric from a minimum eigenvalue of about -1529 Ha to 74 Hermitian positive-definite matrices, while keeping the accepted producer exchange unchanged.

## User interface

Four INPUT keys are added:

```text
out_librpa_2d_coulomb_method       ewald
out_librpa_2d_direct_ecut          0
out_librpa_2d_direct_kz_order      0
out_librpa_2d_direct_gamma_order   8
```

`out_librpa_2d_coulomb_method` accepts `ewald` and `direct_mixed_fourier`.  `ewald` preserves all current behavior.  `direct_mixed_fourier` is valid only with `rpa 1`, `out_librpa_reader_version 1`, and `exx_ewald_dimension 2`.  It also requires a positive reciprocal cutoff in Ry, a positive kz quadrature order, and a positive even Gamma-plane order.  Invalid combinations fail before the expensive producer stage; there is no automatic fallback.

The WS2 values `110 Ry`, `64`, and `8` are a tested starting point, not universal defaults.  Keeping cutoff and kz order at zero unless the direct method is selected prevents an apparently high-accuracy but unconverged choice from being applied silently.

## Numerical construction

For each full-grid q point, form the auxiliary Fourier amplitudes B_mu(G_parallel,kz) and accumulate

```text
V_mu,nu(q) = sum_p w_p conj(B_mu(p)) B_nu(p),  w_p > 0.
```

Regular non-Gamma planes use reciprocal-lattice enumeration up to `out_librpa_2d_direct_ecut` and positive kz quadrature.  The Gamma plane uses an explicit positive polar quadrature of even order.  The matrix is therefore a Gram matrix up to floating-point roundoff.  Hermiticity is enforced only by averaging V with its adjoint after accumulation; negative eigenvalues are not clipped.

The direct matrices replace only `v1_coulomb_full_iq_*`.  The accepted ABFS, shrink transformation, Cs/sinvS tensors, exchange contraction, q ordering, and the analytic strict-2D head normalization remain on the established producer route.  This separation makes a legacy-versus-direct comparison a single-variable test of the exported full Coulomb metric.

## Provenance and validation

Rank zero writes `librpa_2d_coulomb_method.dat` with method, cutoff, kz order, Gamma order, q count, auxiliary dimension, and executable/source identity available at runtime.  The existing `librpa_2d_coulomb_head.dat` remains unchanged so LibRPA reader-v1 compatibility is preserved.

Acceptance requires all of the following:

- exact reader-v1 headers, byte sizes, q ordering, and auxiliary dimension;
- finite matrices with Hermiticity residual reported for every q;
- no materially negative eigenvalue; numerical roundoff is reported, never hidden by clipping;
- convergence of the minimum eigenvalue and near-Fermi EXX/GW bands with respect to reciprocal cutoff, kz order, and Gamma-plane order;
- unchanged producer exchange and Cs/sinvS artifacts relative to the matched legacy producer;
- separate reporting of high-empty-state Padé/QPE anomalies, which the direct Coulomb route does not claim to cure.

## Compatibility

Default input produces byte-for-byte the existing Ewald route.  Legacy reader output is not extended.  Environment variables used by the WS2 prototype are deliberately not part of the production interface.

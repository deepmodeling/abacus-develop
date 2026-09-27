# SCCS discrete energy derivatives

The polarization equation is discretized before differentiation. On a fixed
uniform real-space grid, let `q` be the solute charge, `t = q + p` the screened
charge, `C` the symmetric Coulomb operator, and `G` its implemented field-gradient
operator. `G` includes the analytic PCC gradient when PCC is enabled. Let `D` be
the periodic spectral gradient applied to the cavity electron density `n`.

The converged polarization equation is

```
f(n) = log(epsilon(n))
g = f'(n) D n
a = g / (4 pi)
L = I - a dot G
L t = q / epsilon.
```

The electrostatic correction differentiated here is

```
E = 1/2 <q, C(t-q)> + c <1, t-q>,
```

where the inner product includes the real-space volume element. The coefficient
`c` is zero except for the PCC2D ionic-shape correction. The independent vacuum
PCC energy and non-electrostatic cavity terms retain their separate derivatives.

Solve the discrete adjoint equation

```
L^T lambda = Cq + 2c,
L^T lambda = lambda - G^T(a lambda).
```

The exact derivatives of this discrete energy at a converged polarization state
are

```
dE/dq = U = (Ct + lambda/epsilon)/2 - Cq - c
dE/dn = -U - lambda q epsilon'(n)/(2 epsilon^2)
         + [D^T(f'(n) lambda Gt) + f''(n) lambda (Dn dot Gt)]/(8 pi).
```

`charge_potential` stores `U` for the explicit ionic force;
`electron_potential` stores `dE/dn` for the Hamiltonian and LCAO Pulay force.
`reaction_potential` remains `C(t-q)` for physical diagnostics. At finite grid
resolution these are distinct: replacing `U` by the reaction potential or using
the continuum field-squared cavity derivative assumes identities that do not
hold for the discrete polarization equation.

For PCC2D, the additional explicit ionic derivative of `c` must also be included
in the force. Differentiating only the polarization charge in `c <1,p>` omits
the derivative of the centered smooth and point ionic moments.

The periodic part of `G^T` is `C_periodic (-div_periodic)`. The analytic PCC part
is transposed separately using integrated vector-field moments. Applying a
periodic FFT derivative to the PCC polynomial would not give this transpose.

The transpose system is nonsymmetric and is solved with BiCGSTAB. Its true
residual must satisfy both configured SCCS residual tolerances; the configured
iteration limit also applies to this solve. The existing mixing options control
the original polarization iteration. A compatible previous adjoint solution is
reused as the initial guess. `sccs_debug 2` reports the adjoint iteration count
and residuals. The ordinary SCCS time includes both solves, while `SCCS_ITER`
continues to count polarization iterations.

This change adds no INPUT parameters and changes no defaults. Existing residual
and iteration thresholds also govern the adjoint solve. It changes the finite-grid
electronic potential and therefore can change self-consistent densities and
energies. It also adds computational cost. It does not establish convergence
with respect to the grid or eliminate residual grid-dependent translation error.
For LCAO, the grid convergence parameter in the SN2 validation is `ecutwfc`.

Regression coverage includes the transpose inner-product identity and energy
directional derivatives for periodic, PCC0D, and PCC2D operators, with a sharp
water dielectric cavity. Full SCF finite differences are additionally required
to verify the Hamiltonian/Pulay/ionic-force integration.

## Gradient-only polarization iterations

The polarization residual uses `grad(C t)` and does not use the scalar
potential `C t`. `CoulombOperator::apply_gradient` therefore permits a
specialized path that omits the scalar-potential inverse FFT and the PCC
scalar-potential evaluation during inner iterations. Periodic, PCC0D and
PCC2D implementations preserve the full operator's gradient arithmetic.
Other operators can use the complete-field fallback.

Both convergence and iteration-limit exits reconstruct the complete field at
the returned polarization charge. Energy and discrete-adjoint evaluation thus
receive the same scalar potential and gradient as before. Gradient sizes and
finite values are checked every iteration; scalar-potential validity is checked
when the complete terminal field is constructed. Nonfinite-gradient exits do
not promise a usable field. Solver tolerances and adjoint equations are unchanged.


## Reusing Pulay residual products

Pulay's Gram matrix contains global inner products of immutable residual
history vectors. A solve-local cache retains products between old vectors;
only pairs involving a newly appended residual are recomputed. Removing the
oldest history removes the corresponding matrix row and column. Adaptive
and residual-growth restarts clear the cache together with both histories.
The dense regularized solve, coefficient safeguards, reduction operation and
per-product summation order are unchanged. The cache is not retained across
SCF calls, and linear/Anderson paths are unchanged.

## Chain-gradient discretization

The cavity gradient is evaluated as `(epsilon'(n)/epsilon) Dn`, rather than
`D log(epsilon)`. These operations are not interchangeable on a finite
spectral grid. Both the forward polarization equation and its discrete
adjoint derivative must change together. The divergence in the density
variation acts on `f'(n) lambda Gt`, and the local `f''(n)` term accounts for
the density dependence of the chain coefficient. It vanishes outside the
cavity transition interval. `D^T` is the negative periodic divergence;
`G` and `G^T` still include the analytic PCC terms.

Cavity density and parameters are passed explicitly to the adjoint routine;
no new global configuration access or INPUT switch is added. The ionic
source remains the original local-pseudopotential charge. This change does
not introduce Gaussian ions, an Environ CG solver, or continuum no-adjoint
force formulas. The independent PCC correction and non-electrostatic terms
retain their existing behavior.

INPUT names, meanings and defaults are unchanged, so generated INPUT
metadata needs no update. Finite-grid self-consistent densities, energies
and forces can change; older numerical references must be reassessed rather
than silently reused. Common inner RMS/max thresholds of 1e-11/1e-9 were
adequate in the tested water cases; this is not a universal grid/SCF accuracy
guarantee. Focused tests cover analytic density-mode chain gradients,
source/cavity/coupled energy variations for periodic/PCC0D/PCC2D, and
reuse of an already-converged adjoint solution.

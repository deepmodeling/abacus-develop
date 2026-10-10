# LCAO radial interpolation and force consistency

This diagnostic compares the Fe2-x analytical force with a central difference
of fully converged SCF energies. It uses a distorted periodic Fe2 cell, SOC
(`nspin 4`), **DeltaSpin off** (`sc_mag_switch 0`), no Hubbard U, a 2x2x2
k mesh, and either 100 or 200 Ry. The complete INPUT, STRU and KPT are here.
Initial moments are (0,0,2.2) Bohr magnetons; they are not constrained targets.

Run from the repository root, with an MPI + LCAO + Libxc executable:

```sh
python3 tests/integrate/test_lcao_radial_force.py \
  --abacus /absolute/path/to/abacus --work-dir /absolute/path/to/results \
  --device gpu --launcher 'mpirun -np 1' --cutoff 100
```

Repeat with `--cutoff 200`. The CPU option `--device cpu` selects
`scalapack_gvx`; GPU selects `cusolver`. The runner enforces one OpenMP and
one OpenBLAS thread and creates three independent cold-start calculations
(base, Fe2-x minus 0.0005 Angstrom, Fe2-x plus 0.0005 Angstrom). It preserves
every input/log and writes energies, raw printed forces, net forces,
convergence data and the comparison to `result.json`. It exits nonzero on
execution, parsing, convergence or force-comparison failure.

The convergence gates are the SCF-converged marker, finite observables,
DRHO < 1e-10, and absolute final `DeltaE_womix` < 1e-9 eV. The comparison is

$$
F_{\mathrm{FD}}=-\frac{E(R_{2x}+h)-E(R_{2x}-h)}{2h},\qquad h=0.0005\ \mathrm{Angstrom}.
$$

`E` is `FINAL_ETOT_IS`, including the smearing contribution, not the
extrapolated zero-temperature energy. The coordinate displacement accounts
for the STRU lattice constant and ABACUS's `BOHR_TO_A = 0.5291770`.
With symmetry disabled, ABACUS removes the mean force before printing;
the matching unshifted analytical force is reconstructed as

$$
F_{2x}=F_{2x}^{\mathrm{printed}}+F_x^{\mathrm{net}}/2.
$$

The acceptance threshold is `abs(F_FD - F_2x) < 1e-4 eV/Angstrom`.
This diagnostic is deliberately manual, not registered in the default
integration suite: it needs three tightly converged SOC SCFs per cutoff.
The small synthetic regression is part of the existing CTest target
`MODULE_AO_ORB_atomic_lm_test` and needs no pseudopotential or orbital asset.

## Why values and derivatives must come from one spline

`Numerical_Orbital_Lm::extra_uniform` previously populated `psi_uniform`
with `Uni_RadialF`, but populated `dpsi_uniform` with the derivative of a
different cubic spline. The spline routine already computed the matching
values, which were discarded in a temporary buffer. This fix stores those
values directly and retains zero value padding at and beyond the cutoff.
It preserves the existing spline boundary conditions and derivative table.

To derive the consistency condition, let the original radial knots be
$x_j,x_{j+1}$, interval width $H=x_{j+1}-x_j$, samples $f_j,f_{j+1}$ and
spline second derivatives $M_j,M_{j+1}$. Define
$a=(x_{j+1}-r)/H$ and $b=(r-x_j)/H$. The existing spline evaluator returns

$$
S(r)=a f_j+b f_{j+1}+\frac{H^2}{6}
\left[(a^3-a)M_j+(b^3-b)M_{j+1}\right].
$$

Since $a'=-1/H$ and $b'=1/H$, differentiating term by term gives

$$
S'(r)=\frac{f_{j+1}-f_j}{H}
-\frac H6(3a^2-1)M_j+\frac H6(3b^2-1)M_{j+1},
$$

exactly the derivative returned by the same evaluator. Thus storing both
outputs establishes `psi_uniform[i] = S(r_i)` and
`dpsi_uniform[i] = S'(r_i)` inside the support.

The orbital-value path in `GintAtom::set_phi` uses cubic Hermite
interpolation on a fine interval of width $h$, with $t=(r-r_i)/h$:

$$
\mathcal H(t)=h_{00}(t)y_i+h h_{10}(t)d_i
+h_{01}(t)y_{i+1}+h h_{11}(t)d_{i+1},
$$

$$
h_{00}=2t^3-3t^2+1,\quad h_{10}=t^3-2t^2+t,\quad
h_{01}=-2t^3+3t^2,\quad h_{11}=t^3-t^2.
$$

If this fine interval lies in one original spline interval, both $S$ and
$\mathcal H_{\mathrm{new}}$ are cubic polynomials with the same value and
first derivative at each endpoint. Their difference therefore has two
double roots. A polynomial of degree at most three cannot have those four
roots unless it is identically zero. Hence
$\mathcal H_{\mathrm{new}}=S$ and $\mathcal H'_{\mathrm{new}}=S'$ there.

For the old values $P(r_i)$, put $\epsilon_i=P(r_i)-S(r_i)$ while retaining
$d_i=S'(r_i)$. Linearity of Hermite interpolation gives

$$
\mathcal H_{\mathrm{old}}-S=h_{00}\epsilon_i+h_{01}\epsilon_{i+1}.
$$

Using $dt/dr=1/h$, $h'_{00}=-6t(1-t)$ and
$h'_{01}=6t(1-t)$ then gives the extra derivative term

$$
\frac{d\mathcal H_{\mathrm{old}}}{dr}-S'
=\frac{6t(1-t)}h(\epsilon_{i+1}-\epsilon_i).
$$

Its maximum magnitude is $3|\epsilon_{i+1}-\epsilon_i|/(2h)$, attained
at $t=1/2$. This term vanishes for the corrected tables. There is no claim
that it diverges as $h$ decreases: the endpoint error difference also
changes with $h$.

The GPU counterparts in `kernel/phi_operator_kernel.cuh` use the same
Hermite value and four-point value/slope formulas.

The force path `GintAtom::set_phi_dphi` separately uses four-point polynomial
interpolation of the values and slopes. When all four samples lie in one
original spline interval, cubic interpolation reproduces both the cubic
$S$ and its quadratic derivative $S'$ exactly. The old pair of tables
cannot generally satisfy that same value/derivative contract.
For an orbital $\phi(\mathbf r-\mathbf R)=S(\rho)Y_{lm}(\hat\rho)$,
$\partial\phi/\partial\mathbf R=-\nabla\phi$ contains $S'$.
For example, a fixed-potential grid term $e=wV\phi^2$ differentiates to
$\partial e/\partial R=2wV\phi\,\partial\phi/\partial R$;
using the derivative of another radial representation violates this chain
rule. This example explains the mechanism, not a derivation of every DFT
force term. The full-SCF difference is the separate physical check.

The proof is local to a single original spline interval, away from the
support boundary. It does not prove exact equality across spline knots,
zero numerical force error, cutoff convergence, or correctness of every
force/stress term. The existing derivative padding is unchanged; only
orbital **values** are explicitly zeroed outside the support.

## Regression and assets

`NumericalOrbitalLmUniformGrid.DerivativeMatchesInterpolatedValues` samples
$r^2(2-r)^2e^{-r}$ on a 0.02-Bohr mesh and resamples on a 0.001-Bohr mesh.
A five-point centered derivative, exact for a cubic within one spline
interval, independently checks the returned value table against the slope
table to 5e-10. All explicit angular-momentum boundary-condition branches
and the default branch are covered (`l = 0..5`), along with non-vacuous
zero-padding checks. Restoring the old value interpolation must fail.

The numerical diagnostic uses these exact repository assets:

| Asset | SHA256 |
| --- | --- |
| `tests/PP_ORB/Fe.upf` | `b9b7671a2c667a6f87a7228d30fcc205e339b69f5568fe3f3dba567fdce59a19` |
| `tests/PP_ORB/Fe_gga_8au_200.0Ry_4s2p2d1f.orb` | `8ff6a789b1c7636a1d8ce857edb2477aa094aa36c25cdad11fce646b852c484e` |

The pseudopotential is fully relativistic ONCVPSP, PBE, 16 valence electrons,
with nonlinear core correction. The orbital support is 8 Bohr;
`200.0Ry` in its filename describes basis generation, not the INPUT cutoff.
The same basis is used at both cutoffs. This case covers one atomic-force
component with DeltaSpin off; it is not a magnetic-force, stress, PW,
performance, or basis/cutoff convergence test.

## Reference comparison

A matched GPU comparison against upstream develop `97697c26` (v3.11.0-beta10),
with no DeltaSpin force-factor patch in either executable, gave the following.
Both binaries ran on the same V100 allocation, one MPI rank and one OpenMP
thread. All 12 cold-start SCFs met the gates above. Forces are in eV/Angstrom.

| Cutoff (Ry) | Version | Unshifted Fe2-x | Energy derivative | FD minus analytical |
| --- | --- | ---: | ---: | ---: |
| 100 | upstream | -0.1554587696 | -0.1696286481 | -1.41698785e-02 |
| 100 | consistent spline | -0.1554956557 | -0.1555013378 | -5.68205832e-06 |
| 200 | upstream | -0.1436114097 | -0.1436140328 | -2.62312237e-06 |
| 200 | consistent spline | -0.1436115722 | -0.1436087696 | +2.80262347e-06 |

The original 100-Ry result fails the 1e-4 force gate; the corrected result
passes. Both 200-Ry versions pass this component test. This does not establish
cutoff convergence: the force itself still changes between 100 and 200 Ry.
Nor does it restore zero net force. The corrected net x force is approximately
-0.03044932 eV/Angstrom at 100 Ry and -0.00831885 eV/Angstrom at 200 Ry.
Those nonzero sums remain, and their origin is not resolved by this check.
The reported success is consistency of the **unshifted Fe2-x derivative**,
not a claim of translational invariance or acceptance of all force components.

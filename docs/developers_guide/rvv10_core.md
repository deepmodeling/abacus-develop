# rVV10 implementation

`dft_functional GGA_X_RPW86+GGA_C_PBE` selects the semilocal base and
`xc_nonlocal rvv10` adds rVV10 nonlocal correlation with b=6.3 and C=0.0093.
The two choices are intentionally independent: future PBE+rVV10 and
SCAN+rVV10 combinations can change the semilocal base without changing the
nonlocal evaluator or potential component.

The VV10 parameter C=0.0093 was fitted using rPW86-PBE densities to 54
same-species C6 coefficients; the reported mean absolute percentage error was
14%, or 9% after excluding the six metal atoms Li, Na, Mg, Al, Zn and Ga.
The earlier C=0.0089 used LC-omegaPBE densities (Vydrov and Van Voorhis 2010,
Section IV and note 24). These are published model errors, not validation
results for this implementation. The original rVV10 article reports a maximum
S22 binding-energy difference from VV10 of 0.85 kcal/mol before refitting b and
0.16 kcal/mol after refitting to b=6.3; the latter does not isolate the kernel
approximation at fixed parameters.

The experimental SCF path requires Libxc, norm-conserving pseudopotentials,
CPU double precision, a PW basis, Kohn-Sham SCF and `nspin 1` or `nspin 2`.
The existing PW FFT distribution may use an MPI pool; scalar nonlocal outputs
are reduced over that pool and the potential remains local to each rank. Force,
stress, relaxation, MD, noncollinear spin and additional dispersion corrections
are unsupported. Use `gamma_only 0`: the evaluator requires a full-complex FFT
layout, although ordinary Gamma-point k sampling is supported.
For a PW input, the existing parser resets `gamma_only 1` to zero before the
functional checks; the evaluator still rejects a reduced-gamma `PW_Basis` passed
directly by a caller.
The existing PW parser normalizes `gamma_only 1` to zero and writes a Gamma KPT
file. Unsupported calculation modes are rejected rather than evaluated with a
partial correction.

## SCF integration and reused interfaces

`Potential::get_pot_type("xc")` always creates the canonical `PotXC` component.
When `xc_nonlocal=rvv10`, the PW Hamiltonian registers a second component,
`Potential::get_pot_type("rvv10")`, after `xc`. `PotRvv10` owns one
`Rvv10::Evaluator`, retaining the kernel table and interpolation splines
between SCF iterations. `XC_Functional::set_xc_type` selects the semilocal
Libxc identifiers, and `PotXC` evaluates that base. Libxc does not compute the
nonlocal correction.

For each density update, `PotRvv10` constructs total density
`rho[0] + rho_core` for `nspin=1`, or `rho[0] + rho[1] + rho_core` for
`nspin=2`, obtains the nonlocal energy and potential, and adds them to
the outputs already produced by `PotXC`. The component adds its contribution to
`etxc`, `vtxc` and `v_eff`, then the combined result satisfies

\[
vtxc=\int n_{\mathrm{val}}(\mathbf r)v_{xc}(\mathbf r)\,d\mathbf r.
\]

The nonlinear core correction (NLCC) therefore enters the energy and potential
through total density, while the potential contraction uses valence density only.
The existing SCF energy machinery consumes `etxc` and `vtxc`; no separate
post-SCF energy correction is added. Registration order is therefore part of
the contract: `xc` must run before `rvv10`, because `PotXC` initializes the
semilocal energy and potential and `PotRvv10` only adds to them.

The numerical code is in `source/source_hamilt/module_xc/xc_rvv10.{h,cpp}` and
`xc_rvv10_pw.{h,cpp}`, compiled into `xc_`. The SCF adapter is in
`source/source_estate/module_pot/pot_rvv10.{h,cpp}`. The evaluator takes explicit
dependencies:

```cpp
Evaluator(b, C).evaluate(pw, total_density, valence_density)
```

Here `pw` is an initialized density `PW_Basis`, normally `Charge::rhopw`.
The caller supplies the core-inclusive density; the result contains energy,
`vtxc` and a scalar potential on the real-space grid. Evaluation reuses
`PW_Basis::real2recip/recip2real`, `XC_Functional::grad_rho/grad_dot`, and
`ModuleBase::CubicSpline`. The spline boundary conditions are explicitly natural
(zero second derivative), rather than ABACUS's default not-a-knot conditions.
File parsing and grid distribution are caller responsibilities, not evaluator
dependencies. External density comparisons must check voxel vectors, grid
ordering and every coordinate; matching electron counts alone cannot detect
a permutation of grid values.

## Units and discrete energy

Density is in electrons/bohr\(^3\), \(\sigma=|\nabla n|^2\), and \(|G|\) is in
bohr\(^{-1}\). Nonlocal energies and potentials use Rydberg units. In this
convention the local quantities are

\[
\omega^2=\frac{16\pi n}{3}+\frac{4C\sigma^2}{n^4},\qquad
\kappa=3\pi b\left(\frac{n}{9\pi}\right)^{1/6},\qquad
q=\frac{\omega}{\kappa},\qquad
D=\frac{n}{\kappa^{3/2}},\qquad
\beta=\frac{1}{16}\left(\frac{3}{b^2}\right)^{3/4}.
\]

At b=6.3, beta is 0.00900969593138887 Ry per electron, equivalently
0.00450484796569444 Ha per electron. It multiplies electron number, not density
alone. The analytic uniform-gas cancellation refers to the unsaturated,
infinite-radius continuum kernel; it is not an exact identity of the finite
kernel table below.

The interpolation argument is

\[
q_0=\max\left\{q_{\min},\;
q_{\mathrm{cut}}\left[1-\exp\left(-\sum_{m=1}^{12}
\frac{(q/q_{\mathrm{cut}})^m}{m}\right)\right]\right\},
\quad q_{\min}=10^{-4},\quad q_{\mathrm{cut}}=0.5\;\mathrm{bohr}^{-2}.
\]

With cardinal splines \(P_\alpha\), the 20 channel fields are
\(\theta_\alpha=D P_\alpha(q_0)\). The real-space kernel paired with this
Rydberg definition of \(D\) is

\[
K(q,q',R)=\frac{-24}
{(1+qR^2)(1+q'R^2)[2+(q+q')R^2]}.
\]

Using ABACUS's forward-normalized Fourier transform, the evaluator forms
\(u_\alpha(G)=\sum_\beta\widetilde K_{\alpha\beta}(G)\theta_\beta(G)\) and

\[
E_c^{\mathrm{nl}}=\frac{\Omega}{2}\sum_{G,\alpha}
\operatorname{Re}[\theta_\alpha(G)^*u_\alpha(G)]
+\beta\Delta V\sum_r n(r),\qquad \Delta V=\frac{\Omega}{N_{\mathrm{grid}}}.
\]

No extra grid-size factor belongs in the reciprocal contraction. After inverse
transforms, the potential is

\[
v_c^{\mathrm{nl}}=\beta+
\sum_\alpha u_\alpha\frac{\partial\theta_\alpha}{\partial n}
-\nabla\!\cdot\left[
2\sum_\alpha u_\alpha\frac{\partial\theta_\alpha}{\partial\sigma}\nabla n
\right].
\]

`LocalValues::dq0_dsigma` stores \(\partial q_0/\partial\sigma\); the factor
two belongs in the vector field. The gradient and divergence use the same
density-basis G cutoff. The beta energy uses total density, and the beta
potential contributes to the valence-only `vtxc` contraction.

On the actual grid, the necessary identity is
\(\langle G_h\delta n,h\rangle=-\langle\delta n,D_hh\rangle\), with
the volume-weighted grid inner product. The implementation applies `grad_dot`
to the completed vector field; it does not move its spatially varying
coefficient outside the divergence using a continuum product rule.

## Cutoffs, kernel table and memory

For \(n\le10^{-12}\), the nonlocal weight and its derivatives are zero; the
linear beta term remains active. The energy is not differentiable across this
hard density cutoff. At the lower \(q_0\) clip, ABACUS follows the pinned QE
implementation: \(q_0\) is clamped for interpolation but its analytic
derivatives are retained. For
\(q/q_{\mathrm{cut}}\ge2\), saturation rounds to \(q_{\mathrm{cut}}\) in double
precision and its derivatives are zero. Nonzero small-gradient derivatives are
retained, including when \(\sigma\le10^{-12}\); only \(\sigma=0\) gives a zero
gradient derivative. The hard density cutoff and QE's independent small-gradient
mask remain explicit numerical conventions rather than smooth limits.

`KernelTable` tabulates 210 independent channel pairs. It evaluates
\(4\pi\int_0^{100}R^2K\,\mathrm{sinc}(GR)\,dR\) with 1024 trapezoid intervals,
including the half-weight at the final endpoint and the nonsingular \(G=0\)
limit. There are 1025 G knots with spacing \(2\pi/100\); queries require
\(0\le|G|<1024(2\pi/100)\). Extrapolation is rejected.

The fixed 100-bohr radial cutoff follows the reference discretization. At the
smallest diagonal q, \(R_{\max}\sqrt q=1\). Writing \(x=R_{\max}\sqrt q\),
the zero-G diagonal transform is

\[
\widetilde K(q,q,0;R_{\max})=-6\pi q^{-3/2}
\left[\arctan x-\frac{x(1-x^2)}{(1+x^2)^2}\right].
\]

At x=1 this is exactly half its infinite-radius value \(-3\pi^2q^{-3/2}\).
That is a kernel-element tail fraction, not a 50% total-energy error: channel
weights and the actual density determine its contribution. Tests use the
finite-radius integral. Agreement with this table does not establish
infinite-radius convergence or exact uniform-gas cancellation.

Radial point count, cutoff and q-channel count are separate numerical choices.
Increasing the radial point count at fixed cutoff refines the quadrature but
does not reduce the G-knot spacing. Increasing the cutoff also refines that
spacing, so on-knot quadrature and off-knot spline errors must be distinguished.
Refined q meshes change the channel interpolation; more nodes alone do not
prove convergence. Such experiments must preserve the production preset and
report energy and potential sensitivity, not silently change beta or references.
The table's derivative is the derivative of its spline interpolant; stress
requires additional terms and has not been implemented.

SCF energy/potential use `KernelTable::values`, which requests only function
values from the existing `CubicSpline::multi_eval`. The full `evaluate` API
retains the spline derivative for numerical checks. Both paths share the same
range check, interpolation and symmetric matrix expansion; no second kernel
implementation or cache is introduced.

Construction temporarily uses about 8 MiB for shared transform weights; the
stored splines occupy about 3.3 MiB. The table is independent of b, C and the
cell. Evaluation needs \(O(20N_G+N_r)\) scratch storage, retaining only one full
20-channel reciprocal array. As channels stream through the FFTs, each grid
point evaluates only the current channel with `CubicSpline::multi_eval(1, ...)`.
This avoids recomputing 19 unused channels on every pass without adding a
60-double-per-point real-space cache. Large-cell scaling remains to be measured.

The evaluator rejects reduced-gamma layouts, incompatible local array sizes,
nonfinite inputs, negative sigma, invalid parameters and G outside the table.
Unrepresentable local scales raise an overflow error. Distributed calls use the
existing pool collectives for scalar reductions. `PW_Basis` owns mutable FFT
buffers, so concurrent calls sharing one basis are not thread-safe.

## Tests and external comparison

Use the repository's CMake/CTest workflow. For a CPU serial build with its usual
dependencies configured:

```bash
cmake -S . -B build-rvv10 -DBUILD_TESTING=ON \
  -DENABLE_MPI=OFF -DENABLE_LCAO=OFF -DENABLE_OPENMP=OFF -DENABLE_LIBXC=ON
cmake --build build-rvv10 --parallel 2 --target abacus_pw_ser \
  MODULE_HAMILT_XCTest_RVV10_CORE MODULE_HAMILT_XCTest_RVV10_PW \
  MODULE_ESTATE_RVV10_POTENTIAL
OMP_NUM_THREADS=1 ctest --test-dir build-rvv10 --output-on-failure --no-tests=error -R RVV10
```

The core tests exercise local derivatives, clipping, spline interpolation and
the kernel table. The PW tests use actual ABACUS FFT and derivative interfaces
to check normalization, energy variations, core-inclusive density and
valence-only contractions. An independent scalar/vector-field test checks the
discrete gradient/divergence adjoint in a nonorthogonal cell, including modes
removed by the density cutoff and a nonzero contraction that tests the sign.
A periodic slab test combines material, dilute tails
and inactive vacuum on two nonorthogonal grids; its finite differences stay
away from the hard density cutoff. The potential tests check the full XC energy
derivative, repeated updates and composition of `PotXC` followed by `PotRvv10`
in one process using the actual adapters. They also check the `nspin=2`
contract: the nonlocal functional uses the summed spin density and adds the
same scalar potential to both spin rows. The SCF regression exercises the
production registration path, while the native tests keep the component
contract independent of global input state. Their expected values come from
direct calls to the real adapters; the independent energy-derivative tests
remain necessary to verify the mathematics.

The SCF driver checks He and NLCC Si energies, a PBE control, unsupported INPUT
modes and rejection after reading an actual ultrasoft UPF. It also verifies
that unsupported Libxc nonlocal-functional names stop before SCF rather than
silently receiving only a semilocal contribution. Inputs, logs and a
JSON summary are preserved. See `tests/rvv10/README.md` for tolerances and launcher
configuration. The workflow selects serial, no-Libxc and MPI-build tests. The
MPI build compares a distributed rVV10 SCF with its serial reference and checks
the controlled unsupported-mode diagnostics. Test counts and hosted CI results
belong in the PR's fresh verification record.

External comparisons performed on 2026-09-13 used the unchanged QE numerical
source identified below. On identical synthetic densities and on three
graphene SCF densities, the largest absolute nonlocal energy difference was
\(3.5\times10^{-15}\) Ry. Complete rVV10 SCF total-energy differences were
\(4.30\times10^{-5}\) eV for periodic He, \(2.14\times10^{-5}\) eV for a two-atom
Si primitive cell with NLCC, and \(2.34\text{--}3.01\times10^{-5}\) eV for two-atom
graphene cells with layer spacings of 20, 30 and 40 bohr. Both programs used the
same norm-conserving pseudopotentials and matched cells, FFT grids and k points,
with 40/160 Ry wavefunction/density cutoffs. QE's native semilocal XC and
ABACUS's Libxc base contributed to these full-program differences. Energy
comparisons first used Ry to remove differences in the programs' eV constants.

A 2026-09-15 check used the S22 water-dimer geometry distributed with
[ASE 3.26.0](https://gitlab.com/ase/ase/-/raw/3.26.0/ase/data/s22.py),
H/O ONCV PBE 1.0 pseudopotentials without NLCC, a 24-bohr cubic cell,
a 100-cubed density grid and 40/160 Ry cutoffs. Twelve converged SCFs covered
the dimer and both monomers at their original positions for RVV10 and its
semilocal base, in both programs. The RVV10 binding energies were
-0.0182223856070962 Ry (ABACUS) and -0.0182210818350868 Ry (QE).
The difference in the RVV10-minus-base binding-energy change was
-7.16959008e-8 Ry; this includes the self-consistent density response.
This single hydrogen-bonded geometry is a software comparison, not an S22-wide
accuracy test or a cutoff/cell-converged isolated-dimer prediction. The archived
binary predates only the subsequent XC-name ownership fix; the rVV10 numerical
sources are unchanged. Archive inputs, raw outputs and binary hashes with these
results rather than assigning them to a later executable without rerunning.

Two further 12-SCF matrices used the post-fix executable, the same 24-bohr cell
and an actual 120-cubed grid at 40/160 and 60/240 Ry cutoffs. RVV10 binding-energy
differences (ABACUS minus QE) were -5.4020e-7 and -2.8618e-7 Ry, respectively.
At 40 Ry, increasing the grid from 100 to 120 changed the RVV10 binding energy
by -0.211 meV in ABACUS and -0.222 meV in QE. At fixed grid 120, increasing the
cutoff from 40 to 60 Ry changed it by +0.00302 and -0.000431 meV, but changed
the semilocal-base binding energy by -0.342 and -0.344 meV. The latter also
changes the RVV10-minus-base effect; the small RVV10 total change alone is not
proof that every contribution is converged. All 24 final energies, actual grids
and QE basis cutoffs were reread from the archived logs/XML. These are two-level
parameter sensitivity checks, not a cell-convergence or complete S22 result.

Pointwise potentials differ from this original QE reference. Graphene vacuum
regions had maximum differences of 0.657--1.603 Ry (8.94--21.82 eV), with
density-weighted RMS differences of 0.00631--0.00771 Ry. Small total-energy or
`vtxc` differences therefore do not establish pointwise potential agreement.

Separate, explicitly modified QE diagnostic modules isolated two numerical
choices: projection before the divergence and the small-gradient derivative
cutoff. Matching these choices reduced the
graphene maximum potential difference to at most \(3.31\times10^{-8}\) Ry;
the dominant difference in these cells was QE's suppression of the gradient
derivative for \(|\nabla n|^2\le10^{-12}\). Central finite differences of the
*original* QE energy supported the ABACUS potential at the sampled discrepant
points. The diagnostic agreement is distinct from the original-QE comparison;
it is not a result obtained with unchanged QE.

A further 2026-09-15 check linked evidence-only wrappers against the unchanged
QE module and the current ABACUS native FFT objects. Seven prescribed densities
on an 8x10x12 grid in a 10x12.5x15-bohr cell covered ordinary density with fixed
core, ordinary uniform density, the lower q clip, gradients below and across
the QE cutoff, ordinary nonzero gradients, and inactive density. The density
basis contained 245 G vectors at a 4-Ry cutoff. Perturbations contained only
the constant or lowest Fourier modes within that basis, so this check separates
local derivative conventions from the high-G divergence projection difference.

Baseline nonlocal energies agreed within 8.1e-16 Ry. Nine dimensionless central
difference steps from 0.1 to 1e-5 were compared with the potential contraction
sum(v*delta_n*dV), with fixed core density. All seven ABACUS contractions matched
both programs' energy differences within the preset 1e-6 Ry normalized tolerance
for at least three consecutive steps. Normalization divides by
sum(abs(delta_n)*dV); without it, dilute perturbations can produce misleadingly
small absolute energy errors. The archived original-QE residual for the
lower-clip row predates the present QE-parity correction and must not be reused
as current validation. The remaining archived residuals concern the independent
small-gradient mask and hard density cutoff. The other four cases passed the
same check. All nine steps and all grid potentials are retained in the external
`directional-fd-round6-20260915` evidence, but the lower-clip case must be
rerun before making a current cross-code potential claim.

These checks establish energy agreement for the tested discretizations and
support the derivative convention at sampled points. Material-property
convergence, broader vacuum/NLCC cases, force, stress and GPU require further
work. Spin-polarized and distributed runs still need fresh CI evidence before
the expanded scope is presented as fully validated. The graphene energy still changed by about 7.10 meV/cell from
30 to 40 bohr in both programs, so the cross-code comparison does not establish
vacuum convergence. External inputs, raw results and diagnostic changes should
accompany any PR claim based on these figures.

A subsequent 21-evaluation check held each of those three graphene densities,
cells and grids fixed and applied seven isolated kernel variants. The current
default reproduced every archived ABACUS potential value exactly and agreed
with the original QE nonlocal energies within 3.5e-15 Ry. Increasing radial
intervals from 1024 to 4096 at R=100 changed energies by at most 6.16e-13 Ry.
Increasing R from 100 to 400 at fixed radial step changed them by up to
2.90e-6 Ry; refining the nested q mesh from 20 to 40 channels changed them by
up to 4.68e-6 Ry. The latter two comparisons reached maximum pointwise
potential changes of 0.0166 and 0.0289 Ry, respectively. Their density-weighted
RMS changes were at most 6.33e-5 and 7.50e-5 Ry; no constant offset was fitted.

These are fixed-density sensitivities, not new SCFs, binding energies or
converged-reference errors. Changing R also changes G-spline spacing, and the
20-to-30-to-40 q sequence is not monotonically converged. Low-q electron
populations alone do not bound their nonlocal energy contribution because
channel amplitudes and the pair kernel, not electron counts, set those weights.
The full outputs and two deliberately permuted-data rejection checks are in
the external `material-kernel-round10-20260915` evidence. No production kernel
constant, reference energy or tolerance was changed for these experiments.

Native IntelLLVM 2022.2.1/oneMKL 2022.2.1 builds with OpenMP enabled were also
checked on 2026-09-15. Non-MPI and Intel-MPI builds passed their 7 and 6 focused
CTest targets, respectively, with one OpenMP thread. Eleven repeated SCF
energies differed from the frozen GNU references by at most 3.35e-13 Ry.
On the same three archived graphene densities, native Intel/MKL versus
GNU/FFTW evaluation gave maximum nonlocal-energy and unshifted pointwise
potential differences of 1.39e-17 and 6.31e-8 Ry; the largest density-weighted
potential RMS difference was 8.57e-11 Ry. These satisfy the predeclared 1e-10 Ry
energy and 1e-7 Ry maximum-potential criteria for these fixtures. They do not
establish GPU, multiple-thread or distributed rVV10 support. Runtime library
identities and full outputs are retained in `intel-round11-evidence-20260915`.
The build reuses `abacus::linalg_libs` and the existing Intel strict-FP options
from `abacus_compile_requirements`, rather than specifying private MKL libraries
or overriding floating-point policy in the rVV10 targets.
An isolated core object compiled with only `-fp-model=strict` removed failed
the existing invalid-input test; the native object passed. This negative
control confirms that the test detects lost floating-point requirements.

A CUDA 11.8 CMake build with MPI and LCAO disabled also compiled
`abacus_pw_gpu`; the rVV10 core test passed 1/1. The archived driver then used
that binary with `device=cpu`, so its five-SCF and six pre-existing rejection
cases are CPU-fallback checks only. No GPU was visible during this run, and no
GPU rVV10 evaluation, CUDA-aware MPI evaluation or ROCm evaluation is claimed.
The configure log records a mixed system/conda runtime-library search path; the
library identities and warning must accompany any PR claim based on this build.

## References and provenance

The numerical model is independently implemented from Sabatini, Gorni and
de Gironcoli, [Phys. Rev. B 87, 041108(R) (2013)](https://doi.org/10.1103/PhysRevB.87.041108),
building on Vydrov and Van Voorhis,
[J. Chem. Phys. 133, 244103 (2010)](https://doi.org/10.1063/1.3521275).
The cardinal convolution strategy follows Román-Pérez and Soler,
[Phys. Rev. Lett. 103, 096102 (2009)](https://doi.org/10.1103/PhysRevLett.103.096102).

The numerical mesh and quadrature conventions follow QE commit
`b3af0d0e810a611c455fd80fca6ee523dfaaee95`, `Modules/xc_rVV10.f90`, SHA-256
`7e694bf652b8128ae7baf28e075383ab3e5fa83308ae389d71afe3a71c76c24c`.
QE is an external reference; its GPL routines and spline implementation were
not copied into this source tree. Diagnostic modifications are kept separate
from the unchanged reference implementation.

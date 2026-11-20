# Analytic auxiliary-TDA derivative couplings

The real Gamma-point `ks-lr` solver can calculate analytic nuclear derivative
couplings between the ground state and a TDA excited state, or between two TDA
excited states. These couplings differentiate the spin-adapted or collinear
unrestricted auxiliary determinant expansion built from the canonical
Kohn–Sham orbitals and normalized
TDA amplitudes. They are not full quadratic-response TDDFT couplings.

## Input and scope

Enable the force-response infrastructure with `cal_force 1`, solve the required
roots with `lr_nstates`, and request ascending surface pairs, for example:

```text
esolver_type    ks-lr
gamma_only     1
cal_force      1
lr_nstates     3
lr_nac_pairs   0 1 1 2
lr_nac_spin    singlet
lr_grad_solver cg
lr_grad_thr    1e-10
lr_grad_maxiter 800
```

Surface zero denotes the ground state. Excited-state indices are one-based.
The default empty `lr_nac_pairs` disables NAC calculations. When `lr_nstates`
is automatic (nonpositive), requested surfaces are checked against the resolved
number of LR roots before any state is accessed. The configuration
above illustrates a tighter response solve; the defaults remain `1e-6` Ry and
100 iterations. These two controls apply only to iterative CG orbital response for NAC;
direct response solvers do not use them. A failed certified response solve is an
error, rather than a silently accepted coupling. See the generated INPUT help
for the authoritative parameter defaults, allowed values and restrictions.
If the recursive CG residual converges before the explicitly reapplied residual,
the solver retains its accumulated solution and restarts a correction equation
with the true residual. All correction updates count toward the same limit;
this does not relax the requested accuracy.

Use one real Gamma point, `kpar 1` and `vh_in_h 1`. Hartree-disabled debugging
configurations are rejected because NAC orbital response includes the Hartree kernel. For a nonmagnetic closed-shell reference,
select `lr_nac_spin singlet` or `triplet` (the latter requires `nspin 2`).
For an unrestricted reference, use `nspin 2` and `lr_nac_spin updown`. The latter
uses both alpha and beta orbital spaces at fixed reference M_S; it does not
project onto a pure total-spin eigenstate or enable spin-flip couplings.
Both spin channels must have nonempty occupied and virtual spaces.
An explicitly unrestricted nonmagnetic reference (`lr_unrestricted true`)
uses the same `updown` channel.
The occupied bands must be fully occupied and the virtual bands empty, up to
floating-point roundoff. Fractional occupations, including a thermally smeared
frontier or an odd-electron spin-summed reference, are rejected: the required
occupation-response terms are not part of this auxiliary-determinant model.
The occupied and virtual LR windows must include every occupied and virtual
orbital of the basis: `nocc + nvirt = nbasis` separately in each spin channel.
Zero-padding a truncated amplitude window does not restore the omitted orbital response. Noncollinear/SOC,
spin-flip couplings, general complex k points, DFT+U and pseudopotentials with
nonlinear
core corrections are outside this implementation's scope. Ground–triplet
couplings vanish for a closed-shell spin-adapted reference with a spin-free
Hamiltonian; this zero rule does not apply to the unrestricted `updown` channel.
Ordinary ground-state support
for a functional does not by itself imply that its LR kernel, third XC
derivatives and analytic NAC response are available or verified.
The NAC evaluator rejects unavailable ground-state and LR response combinations,
including a meta-GGA reference with an unrelated available LR kernel. The
resolved reference and response names are checked once at the NAC batch entry,
after the KS/LR reference is available. Ground-state aliases
are resolved before selecting the ground-state response
kernel; use the corresponding canonical name for `xc_kernel`.

## Definition, units and output

For normalized real auxiliary electronic states, the raw coupling is
`d_AB = <Phi_A | d/dR | Phi_B>`. The ETF form removes the antisymmetric
moving-atomic-orbital connection, `d_ETF = d_raw - d_AO`. Reversing the state pair
changes the sign of all three vectors. It leaves the gap-weighted numerators
unchanged because the energy gap also changes sign.

Records are written to the existing running log on its writer rank; this
feature does not create a separate output file for every pair. The existing log
keeps its enclosing solver's lifecycle. The separate-file controls
`out_app_flag` and `out_freq_ion` do not create or select NAC files: the requested
records appear in the running log at every solved LR geometry.
Each `NACV_BEGIN` record contains the one-based geometry index, bra and ket
surfaces, spin label, resolved-gap flag and signed energy gap in Hartree.
Each `NACV_COMPONENT` contains the zero-based atom and Cartesian-axis indices,
then raw, ETF and AO vectors in bohr^-1, followed by raw and ETF numerators in
Hartree/bohr. `NACV_END` terminates the pair record. The public C++ result stores
energies in Ry and numerators in Ry/bohr.

`LR_EIGENPAIR` records give the geometry, spin, excited-state index, excitation
energy and independently reapplied operator-residual norm, in Hartree. A gap
is resolved only when its magnitude exceeds the sum of the two certified
residual norms and a floating-point roundoff bound. For an unresolved gap,
the finite numerators remain available and the divided vectors are zero with
the resolved flag false. These zeros must not be interpreted as physical zero
couplings. No artificial gap clipping is applied.

When NAC output is enabled, existing excitation-energy and amplitude files retain
all nonzero coefficients at full double precision, including magnitudes at or below
`1e-10`. Successful writes are confirmed with their absolute paths in the running
log. An open, write, flush or close failure aborts the output operation.

## Numerical checks and interpretation

Independent finite differences should use the same pseudopotentials, orbital
basis, functional, full LR orbital window and converged electronic states.
Check at least two displacement sizes and verify their convergence. Electronic
phases and state identities must be consistent between geometries. Individual
couplings inside a degenerate manifold depend on the chosen state gauge;
state-by-state comparison without subspace alignment is not unique.

Raw coupling finite differences require electronic overlaps between displaced
geometries. An independent ETF check additionally differentiates the moving-AO
connection numerically and subtracts that numerical connection from the raw
finite difference. Subtracting the tested analytic AO term would not provide
an independent ETF reference. A rigid-translation sum is a useful additional
diagnostic but does not replace these componentwise checks.

Agreement with finite differences tests the implementation at the chosen
numerical settings. Convergence with respect to real-space grid cutoff,
orbital basis and excited-state representation is a separate accuracy question.
Near small gaps, tighter LR and orbital-response residuals may be needed before
reducing the finite-difference step size.

## Review snapshot: independent HSE fixes excluded

This review branch intentionally excludes the separate HSE fixes for the GS
orbital-response sigma floor and the sigma-clipping chain rule in LR gradients.
It retains the upstream HSE numerical model. Results obtained from a snapshot
containing those fixes must not be attributed to this branch. Availability of
an LR kernel does not certify its NAC accuracy; check the same-definition
finite difference before quantitative use.

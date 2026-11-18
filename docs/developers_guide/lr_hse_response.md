# HSE ground-state response fallback in LR gradients

The SCF potential and the LR excitation kernel have different numerical roles.
The HSE LR kernel retains its existing sigma floor and the corresponding force
chain rule. The ground-state orbital-response kernel is constructed separately,
requests no KXC, and first evaluates raw sigma. Its coefficients are shared by
all GS-response terms, including the Z-vector equation and the force potentials.

For HSE06 only, the raw `e_sigma`, `e_rhorho`, `e_rhosigma`, and
`e_sigmasigma` arrays are checked with `isfinite` before spin masking or
accumulation. If any MPI rank reports NaN or infinity, every rank recomputes the
complete HSE response at the same sigma floor as LR (currently 1e-6 bohr^-8).
Finite raw components are replaced together with non-finite ones, preserving a
common evaluation point for the first- and second-order coefficients. Empty
grid ranks participate in both collective checks. If the retry is still
non-finite, construction fails on every rank rather than propagating invalid
coefficients into the linear solve.

A successful retry emits a standard warning with the global non-finite counts,
threshold and units. `KernelXC::gs_sigma_fallback_used()` records this
construction outcome. The fallback uses an approximate GS orbital Hessian:
it does not modify the SCF potential and does not guarantee exact agreement
between analytic excitation-energy gradients and finite differences. Finite
but inaccurate coefficients are outside the detection criterion. A trajectory
can switch between raw and fallback evaluation, so this mechanism is numerical
recovery, not a smooth regularization of the SCF energy functional.

No new INPUT settings or environment overrides are introduced. Non-HSE
functionals and finite raw HSE GS coefficients keep their existing behavior.
The LR sigma Jacobian is unchanged and is not applied to the approximate GS
Hessian as an extra prefactor.

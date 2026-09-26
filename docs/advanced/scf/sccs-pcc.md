# SCCS and PCC

The solvent model and electrostatic boundary correction are selected independently:

| INPUT | Values |
| --- | --- |
| `imp_sol` | `0`: no solvent; `1`: original ABACUS solvent; `2`: SCCS |
| `assume_isolated` | `none`: periodic; `pcc_0d`: cubic isolated molecule; `pcc_2d`: slab periodic in x-z and open along y |
| `sccs_debug` | `0`: no per-step SCCS/PCC output; `1`: iteration/time/energy summary; `2`: all detailed diagnostics |

The existing `makov-payne` option and its aliases remain available through
`assume_isolated`. They cannot be combined with PCC or SCCS.
PCC with the original solvent (`imp_sol 1`) is rejected. Use either vacuum
(`imp_sol 0`) or SCCS (`imp_sol 2`).

## Example

For a slab in SCCS water:

```text
imp_sol          2
assume_isolated  pcc_2d
sccs_preset      water-neutral
sccs_debug       1
```

Choose the solvent preset appropriate to the intended parametrization; the
example does not infer a preset from the system charge. Keep the existing
SCCS cavity and solver parameters when migrating an established calculation.
To remove the solvent while retaining PCC, set `imp_sol 0` and remove SCCS-only
parameters. To retain the SCCS code path at dielectric constant one, use
`imp_sol 2` and `sccs_preset vacuum`.

With `sccs_start_drho > 0`, SCCS polarization is delayed, but PCC remains active
from the initial potential onward. Once activated, SCCS remains enabled in
later electronic and ionic steps. Debug level zero also suppresses deferred
and activation messages. Standalone PCC summaries use `PCC_TIME/s` and
`E_PCC/Ry`; there is no polarization iteration in this path.

## Preset values and solver choices

Only `custom` uses the five physical parameters supplied in INPUT. The other
presets replace them with the effective values below, even if INPUT explicitly
specifies different values. Input range checks still apply to supplied values.
To change a preset's permittivity or cavity parameters, select `custom` and copy
the desired row before modifying it.

| sccs_preset | sccs_epsilon | sccs_rho_min (bohr^-3) | sccs_rho_max (bohr^-3) | sccs_gamma (dyn/cm) | sccs_pressure (GPa) |
| --- | ---: | ---: | ---: | ---: | ---: |
| custom (defaults; INPUT overrides) | 78.3 | 1.0e-4 | 5.0e-3 | 0 | 0 |
| vacuum | 1 | 1.0e-4 | 5.0e-3 | 0 | 0 |
| water-neutral | 78.3 | 1.0e-4 | 5.0e-3 | 47.9 | -0.36 |
| water-cation | 78.3 | 2.0e-4 | 3.5e-3 | 5.0 | 0.125 |
| water-anion | 78.3 | 2.4e-3 | 1.55e-2 | 0 | 0.45 |

The accepted `sccs_mixing_type` strings are exactly `linear`, `pulay`, and
`anderson`. They select damped fixed-point, Pulay DIIS and Anderson acceleration,
respectively. `broyden` and `andersonb` are not accepted INPUT values. This inner
solver setting is independent of the outer electronic `mixing_type`.

All presets retain user control over mixing, tolerances, maximum iteration
count, surface regularization, delayed activation and diagnostic output.
Defaults are `sccs_mixing_type linear`, `sccs_mixing 0.5`,
`sccs_mixing_ndim 8`, `sccs_mixing_adaptive 0`, `sccs_mixing_min 0.1`,
`sccs_mixing_max 0.8`, `sccs_tol_rms 1e-10`, `sccs_tol_max 1e-8`,
`sccs_maxiter 200`, `sccs_surface_eta 1e-8`, `sccs_start_drho 0`,
`sccs_start_nmax 30`, and `sccs_debug 0`.

## Migration

`solvation_model`, `sccs_boundary`, and `pcc_boundary` have been removed.
Old INPUT files containing these names are rejected rather than silently
changing their electrostatics.

- Replace `imp_sol 1` plus `solvation_model sccs` with `imp_sol 2`.
- Replace `sccs_boundary pcc_0d`/`pcc_2d`, or the corresponding
  `pcc_boundary`, with `assume_isolated pcc_0d`/`pcc_2d`.
- Replace `sccs_boundary periodic` with `assume_isolated none`.
- Set `sccs_debug 2` to retain the former detailed boolean debug output.
  Set `sccs_debug 1` for compact per-step summaries.

PCC supports CPU KS-DFT PW/LCAO SCF and fixed-cell relaxation with
`nspin 1` or `2`. SCCS/PCC with ultrasoft pseudopotentials emits a warning and
continues; its numerical compatibility has not been validated. The host ABACUS
USPP restrictions still apply, including the current PW-only restriction.
Stress and simultaneous
external electric/gate fields are unsupported. The 2D open direction is y;
charged-slab absolute energies at different y cell lengths are not directly
comparable. Small FFT grids must be distributed so that each MPI rank has
local real-space grid points; increasing MPI ranks beyond the number of z
blocks is not supported by the current SCCS kernels.

## Implementation

The implementations live under `source/source_hamilt/module_surchem/sccs/`
and `source/source_hamilt/module_surchem/pcc/`. The common PCC host correction
in `pcc/h_corr_pcc.cpp` supplies point-ion/electron moments, energy, potential,
and ionic force to both vacuum and SCCS calculations. SCCS additionally uses
the PCC-aware Coulomb operator during polarization iterations and retains the
2D ionic-shape correction. Its solvent response is therefore solved with the
selected boundary condition.

The `surchem_input.cpp` adapter translates explicit INPUT, cell and parallel
configuration into module settings and validates the open-direction k points.
The electronic solver only invokes this adapter and the correction hooks.
`uses_surchem_correction()` is an internal predicate, not an INPUT parameter:
it is true for either an enabled solvent or a PCC boundary, including
`imp_sol 0` with `assume_isolated pcc_0d`/`pcc_2d`.

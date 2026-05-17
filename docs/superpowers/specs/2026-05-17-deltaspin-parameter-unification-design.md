# DeltaSpin Parameter Unification Design

## Problem Statement

The DeltaSpin module has several parameter-related issues:

1. **`sc_scf_thr` default value contradiction**: `input_parameter.h` declares `1e-3` but `read_input_item_other.cpp` documents `1e-4`. The comment in `input_parameter.h` says "minimum number of outer scf loop" which is semantically wrong (it's a density error threshold, not an iteration count).

2. **Magic number `10.0`**: PW basis uses `sc_scf_thr=10.0` to mean "always pass the drho gate" (enter lambda loop from the first iteration with valid wavefunctions). The `reset_value` logic checks `sc_scf_thr != 10.0` to distinguish this mode, which is implicit and undocumented.

3. **`mixing_restart` interaction with `direction_only`**: When `direction_only=true`, `mixing_restart` is still auto-set to `sc_scf_thr` (~1e-3). But the `direction_only` path enters BFGS from iter=1 without the `sc_scf_thr` gate. If drho naturally drops below 1e-3 during Phase 1, `mixing_restart` triggers `init_mixing()` unexpectedly, polluting `mixing_restart_count` and affecting `mixing_dmr` logic.

4. **Hardcoded Phase 1 duration**: The `direction_only` two-phase strategy hardcodes `iter <= 5` / `iter == 6` with no user control.

5. **Unfinished strategy code**: `lambda_strategy_integration.cpp`, `lambda_update_strategies.h/cpp` and their tests are incomplete and reference undeclared members, creating dead code and maintenance burden.

## Design

### 1. `sc_scf_thr` Default Value and Comment Fix

- Unify default value to `1e-3` across all files.
- Fix `input_parameter.h` comment: `"density error threshold for activating the lambda loop in spin-constrained DFT"`.
- Fix `read_input_item_other.cpp` `item.default_value` from `"1.0e-4"` to `"1.0e-3"`.
- Update `item.description` to reference `sc_scf_thr_mode` for immediate activation.

### 2. New Parameter: `sc_scf_thr_mode`

Enum-type string parameter:

| Value | Semantics | Behavior |
|---|---|---|
| `"threshold"` (default) | Standard drho gate | Lambda loop activates when `drho > 0 && drho < sc_scf_thr` |
| `"immediate"` | Immediate activation | Lambda loop activates from iter>=2 (first iteration with valid wavefunctions). Replaces the old `sc_scf_thr=10.0` convention. |

**Decision logic in `hamilt2rho_single`**:

```
standard DeltaSpin (no direction_only):
  if sc_scf_thr_mode == "immediate":
      if iter > 1: run lambda loop, skip_solve = true
  else:  // "threshold"
      if !converged && drho > 0 && drho < sc_scf_thr: first lambda loop
      elif converged: refine lambda loop
```

**Parameter declaration**:
- `input_parameter.h`: `std::string sc_scf_thr_mode = "threshold";`
- `item.description`: `"Controls when the DeltaSpin lambda loop is activated. \"threshold\" (default): activate when drho < sc_scf_thr. \"immediate\": activate from the first iteration with valid wavefunctions (iter>=2), used for PW basis where the first iteration cannot compute initial magnetic moments."`
- `item.check_value`: validate value is `"threshold"` or `"immediate"`

### 3. `mixing_restart` Auto-Setting Logic Revision

Replace the magic-number-based logic:

```cpp
// Old (read_input_item_elec_stru.cpp reset_value):
if (sc_mag_switch) {
    if (sc_scf_thr != 10.0) {
        mixing_restart = sc_scf_thr;
    } else {
        mixing_restart = scf_thr / 10.0;
    }
}

// New:
if (sc_mag_switch) {
    if (sc_direction_only) {
        mixing_restart = 0.0;            // direction_only: no auto mixing_restart
    }
    else if (sc_scf_thr_mode == "threshold") {
        mixing_restart = sc_scf_thr;     // standard: clean history before lambda loop
    }
    else {  // "immediate"
        mixing_restart = scf_thr / 10.0; // PW: restart on oscillation
    }
}
```

Rationale for `direction_only` path:
- Phase 1 BFGS updates DM directly, not through charge mixing. `mixing_restart` triggering during Phase 1 is harmful.
- Phase 2 transition uses explicit `mix_reset()` (see Section 4), so `mixing_restart=0` is safe.
- Setting `mixing_restart=0` avoids polluting `mixing_restart_count`, which affects `mixing_dmr` logic.

### 4. New Parameter: `sc_dir_phase1_steps`

Integer parameter controlling the Phase 1 duration in `direction_only` two-phase strategy for collinear (nspin=2) calculations.

| Attribute | Value |
|---|---|
| Default | `5` |
| Semantics | Number of SCF iterations for Phase 1 (magnitude constraint via BFGS with direction_only temporarily disabled) |
| Availability | `sc_mag_switch=true && sc_direction_only=true && nspin=2` |
| Validation | `sc_dir_phase1_steps >= 2` |

**Code usage in `esolver_ks_lcao.cpp`**:

```cpp
if (sc_direction_only && nspin == 2) {
    // collinear direction_only: two-phase strategy
    if (iter <= PARAM.inp.sc_dir_phase1_steps) {
        // Phase 1: BFGS with direction_only disabled
        sc.set_direction_only(false);
        sc.run_lambda_loop(iter - 1);
        sc.set_direction_only(true);
        skip_solve = true;
    }
    else {
        // Phase 2: lambda decay + normal SCF
        if (iter == PARAM.inp.sc_dir_phase1_steps + 1) {
            this->p_chgmix->mix_reset();
            this->p_chgmix->mixing_restart_count = 0;
            this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
        }
        // lambda decay: factor = 0.5^(1/3) per step
        ...
    }
}
else if (sc_direction_only && nspin == 4) {
    // non-collinear direction_only: standard sc_scf_thr_mode gate
    // (direction_only works correctly for nspin=4, no two-phase needed)
    if (sc_scf_thr_mode == "immediate") {
        if (iter > 1) { sc.run_lambda_loop(iter - 1); skip_solve = true; }
    } else {
        if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr) {
            sc.run_lambda_loop(iter - 1); sc.set_mag_converged(true); skip_solve = true;
        }
        else if (sc.mag_converged()) {
            sc.run_lambda_loop(iter - 1); skip_solve = true;
        }
    }
}
else {
    // standard DeltaSpin (no direction_only)
    if (PARAM.inp.sc_scf_thr_mode == "immediate") {
        if (iter > 1) { sc.run_lambda_loop(iter - 1); sc.set_mag_converged(true); skip_solve = true; }
        else if (sc.mag_converged()) { sc.run_lambda_loop(iter - 1); skip_solve = true; }
    } else {
        // existing "threshold" logic (unchanged)
        ...
    }
}
```

**Parameter declaration**:
- `input_parameter.h`: `int sc_dir_phase1_steps = 5;`
- `item.description`: `"Number of SCF iterations for Phase 1 (magnitude constraint) in the direction_only two-phase strategy for collinear (nspin=2) calculations. During Phase 1, direction_only projection is temporarily disabled so BFGS can constrain the magnetic moment magnitude. After Phase 1, lambda decays and the system relaxes naturally. Minimum: 2."`

### 5. `drho > 0` Condition Comment

Add comment explaining the condition:

```cpp
// drho > 0: excludes iter=1 where drho has not been computed yet.
// At iter=1, the initial charge density is arbitrary and drho=0,
// so we cannot assess convergence. The lambda loop should only start
// after at least one SCF step has produced a meaningful drho.
// For "immediate" mode, this condition is bypassed by checking iter > 1 directly.
```

### 6. PW/LCAO Oscillation Detection Difference

Document the difference without changing behavior:

```cpp
// NOTE: Unlike the PW path (deltaspin_pw.cpp), the LCAO path does not
// include SCF oscillation detection. This is because:
// - LCAO's BFGS inner loop updates the density matrix directly without
//   going through the charge mixing loop, so oscillation in Broyden
//   mixing history is not a concern during Phase 1.
// - Phase 2 uses normal SCF with charge mixing, and the mix_reset()
//   at the Phase transition already provides a fresh mixing state.
```

### 7. Delete Unfinished Strategy Code

Remove the following files:
- `source/source_lcao/module_deltaspin/lambda_strategy_integration.cpp`
- `source/source_lcao/module_deltaspin/lambda_update_strategies.h`
- `source/source_lcao/module_deltaspin/lambda_update_strategies.cpp`
- `source/source_lcao/module_deltaspin/test/lambda_update_strategies_test.cpp`

Remove from `spin_constrain.h`:
- Any references to `LambdaUpdateStrategy`, `HybridDelayedUpdate`, `AugmentedLagrangianUpdate`, `LinearResponseUpdate`, `LambdaStrategyType`, `strategy_`, `strategy_type_`

Remove from CMakeLists.txt if listed.

## Control Flow Summary

```
hamilt2rho_single:
  if linear_scan: ...

  elif direction_only && nspin==2:
      if iter <= sc_dir_phase1_steps: Phase1 (BFGS, dir_only=false)
      else: Phase2 (lambda decay, mix_reset at transition)

  elif direction_only && nspin==4:
      standard sc_scf_thr_mode gate (direction_only works natively)

  else (standard DeltaSpin):
      if sc_scf_thr_mode == "immediate":
          if iter > 1: run lambda loop
      else:  // "threshold"
          if !converged && drho>0 && drho<sc_scf_thr: first lambda loop
          elif converged: refine lambda loop

mixing_restart auto-setting:
  if direction_only: mixing_restart = 0
  elif threshold:    mixing_restart = sc_scf_thr
  elif immediate:    mixing_restart = scf_thr / 10
```

## Files to Modify

| File | Change |
|---|---|
| `source/source_io/module_parameter/input_parameter.h` | Fix sc_scf_thr comment; add sc_scf_thr_mode, sc_dir_phase1_steps |
| `source/source_io/module_parameter/read_input_item_other.cpp` | Fix sc_scf_thr default_value; add sc_scf_thr_mode, sc_dir_phase1_steps items |
| `source/source_io/module_parameter/read_input_item_elec_stru.cpp` | Rewrite mixing_restart reset_value logic |
| `source/source_esolver/esolver_ks_lcao.cpp` | Restructure direction_only branching; use sc_dir_phase1_steps; add drho>0 comments |
| `source/source_lcao/module_deltaspin/spin_constrain.h` | Remove strategy-related declarations |
| `source/source_lcao/module_deltaspin/lambda_loop.cpp` | Update comments for direction_only behavior |
| `source/source_lcao/module_deltaspin/init_sc.cpp` | Update comments |
| `source/source_pw/module_pwdft/deltaspin_pw.cpp` | Update to use sc_scf_thr_mode instead of sc_scf_thr=10.0 |
| `source/source_lcao/module_deltaspin/lambda_strategy_integration.cpp` | DELETE |
| `source/source_lcao/module_deltaspin/lambda_update_strategies.h` | DELETE |
| `source/source_lcao/module_deltaspin/lambda_update_strategies.cpp` | DELETE |
| `source/source_lcao/module_deltaspin/test/lambda_update_strategies_test.cpp` | DELETE |

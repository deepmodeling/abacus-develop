# DeltaSpin Parameter Unification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify the DeltaSpin parameter scheme by fixing contradictions, replacing magic numbers with named parameters, and correcting the mixing_restart interaction with direction_only.

**Architecture:** Add two new input parameters (`sc_scf_thr_mode`, `sc_dir_phase1_steps`), rewrite the mixing_restart auto-setting logic, restructure the esolver control flow for direction_only, and delete unfinished strategy code.

**Tech Stack:** C++17, ABACUS input parameter system, CMake

---

### Task 1: Delete Unfinished Lambda Strategy Code

**Files:**
- Delete: `source/source_lcao/module_deltaspin/lambda_strategy_integration.cpp`
- Delete: `source/source_lcao/module_deltaspin/lambda_update_strategies.h`
- Delete: `source/source_lcao/module_deltaspin/lambda_update_strategies.cpp`
- Delete: `source/source_lcao/module_deltaspin/test/lambda_update_strategies_test.cpp`
- Modify: `source/source_lcao/module_deltaspin/spin_constrain.h`
- Modify: `source/source_io/module_parameter/read_input_item_other.cpp`

- [ ] **Step 1: Delete the four files**

```bash
rm source/source_lcao/module_deltaspin/lambda_strategy_integration.cpp
rm source/source_lcao/module_deltaspin/lambda_update_strategies.h
rm source/source_lcao/module_deltaspin/lambda_update_strategies.cpp
rm source/source_lcao/module_deltaspin/test/lambda_update_strategies_test.cpp
```

- [ ] **Step 2: Remove strategy-related includes and forward declarations from spin_constrain.h**

The current `spin_constrain.h` has no references to `LambdaUpdateStrategy`, `HybridDelayedUpdate`, etc. (confirmed by grep). However, verify and remove any if present. Also check if `lambda_update_strategies.h` is included anywhere:

```bash
rg "lambda_update_strategies" source/ -l
rg "lambda_strategy_integration" source/ -l
```

Remove any `#include` of these deleted files from any source file.

- [ ] **Step 3: Update sc_lambda_strategy parameter description**

In `source/source_io/module_parameter/read_input_item_other.cpp`, the `sc_lambda_strategy` item description lists `linear_response`, `augmented_lagrangian`, and `hybrid_delayed` as valid strategies. Since these are deleted, update the description to only list the strategies that actually work.

Find the `sc_lambda_strategy` item (around line 200-214) and replace:

```cpp
item.description = R"(Lambda update strategy for spin-constrained DFT:
* bfgs: BFGS quasi-Newton method
* linear_scan: linear sweep of lambda for testing magnetic moment response)";
```

Also update the `item.check_value` if it validates against the deleted strategy names. Add a check:

```cpp
item.check_value = [](const Input_Item& item, const Parameter& para) {
    if (para.input.sc_lambda_strategy != "bfgs" && para.input.sc_lambda_strategy != "linear_scan")
    {
        ModuleBase::WARNING_QUIT("ReadInput", "sc_lambda_strategy must be bfgs or linear_scan");
    }
};
```

- [ ] **Step 4: Update input_parameter.h comment for sc_lambda_strategy**

In `source/source_io/module_parameter/input_parameter.h`, line 604, the comment references deleted strategies:

```cpp
std::string sc_lambda_strategy = "bfgs";  ///< lambda update strategy: bfgs or linear_scan
```

- [ ] **Step 5: Build to verify no broken includes**

```bash
cd build && cmake .. -DBUILD_TESTING=ON -DENABLE_LCAO=ON && make -j4 2>&1 | head -50
```

Expected: no errors about missing `lambda_update_strategies.h` or `lambda_strategy_integration.cpp`.

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "chore(deltaspin): delete unfinished lambda strategy code

Remove lambda_strategy_integration.cpp, lambda_update_strategies.h/cpp,
and their test. Only bfgs and linear_scan strategies remain."
```

---

### Task 2: Fix sc_scf_thr Default Value and Comment

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h:602`
- Modify: `source/source_io/module_parameter/read_input_item_other.cpp:169-184`

- [ ] **Step 1: Fix input_parameter.h comment**

In `source/source_io/module_parameter/input_parameter.h`, line 602, change:

```cpp
double sc_scf_thr = 1e-3;       ///< minimum number of outer scf loop before initial lambda loop
```

to:

```cpp
double sc_scf_thr = 1e-3;       ///< density error threshold for activating the lambda loop in spin-constrained DFT
```

- [ ] **Step 2: Fix read_input_item_other.cpp default_value and description**

In `source/source_io/module_parameter/read_input_item_other.cpp`, find the `sc_scf_thr` item (lines 168-184). Change:

```cpp
item.description = "Density error threshold for inner loop of spin-constrained SCF";
item.default_value = "1.0e-4";
```

to:

```cpp
item.description = "When the charge density error drho falls below sc_scf_thr, the DeltaSpin lambda optimization loop is activated. Should be 10-100x larger than scf_thr. For immediate activation (PW basis), use sc_scf_thr_mode=\"immediate\" instead of setting sc_scf_thr to a large value.";
item.default_value = "1.0e-3";
```

- [ ] **Step 3: Commit**

```bash
git add source/source_io/module_parameter/input_parameter.h source/source_io/module_parameter/read_input_item_other.cpp
git commit -m "fix(deltaspin): correct sc_scf_thr default value and comment

- Unify default value to 1e-3 (was 1e-4 in read_input_item_other.cpp)
- Fix comment from 'minimum number of outer scf loop' to 'density error threshold'
- Update description to reference sc_scf_thr_mode for immediate activation"
```

---

### Task 3: Add sc_scf_thr_mode Parameter

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_other.cpp`

- [ ] **Step 1: Add sc_scf_thr_mode declaration to input_parameter.h**

After the `sc_scf_thr` line (line 602), add:

```cpp
std::string sc_scf_thr_mode = "threshold"; ///< controls when the lambda loop activates: "threshold" (drho<sc_scf_thr) or "immediate" (from iter>=2)
```

- [ ] **Step 2: Add sc_scf_thr_mode input item to read_input_item_other.cpp**

After the `sc_scf_thr` item block (after line 184), add a new item:

```cpp
    {
        Input_Item item("sc_scf_thr_mode");
        item.annotation = "controls when the DeltaSpin lambda loop is activated";
        item.category = "Spin-Constrained DFT";
        item.type = "String";
        item.description = R"(Controls when the DeltaSpin lambda loop is activated.
* threshold (default): activate when drho < sc_scf_thr. The lambda loop starts once the charge density is reasonably stable.
* immediate: activate from the first iteration with valid wavefunctions (iter>=2). Used for PW basis where the first iteration cannot compute initial magnetic moments. Replaces the old convention of setting sc_scf_thr=10.0.)";
        item.default_value = "threshold";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_string(input.sc_scf_thr_mode);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_scf_thr_mode != "threshold" && para.input.sc_scf_thr_mode != "immediate")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_scf_thr_mode must be threshold or immediate");
            }
        };
        this->add_item(item);
    }
```

- [ ] **Step 3: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

Expected: clean build.

- [ ] **Step 4: Commit**

```bash
git add source/source_io/module_parameter/input_parameter.h source/source_io/module_parameter/read_input_item_other.cpp
git commit -m "feat(deltaspin): add sc_scf_thr_mode parameter

New parameter controls when the lambda loop activates:
- threshold (default): activate when drho < sc_scf_thr
- immediate: activate from iter>=2 (replaces sc_scf_thr=10.0 convention)"
```

---

### Task 4: Add sc_dir_phase1_steps Parameter

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_other.cpp`

- [ ] **Step 1: Add sc_dir_phase1_steps declaration to input_parameter.h**

After the `sc_direction_only` line (line 605), add:

```cpp
int sc_dir_phase1_steps = 5;    ///< number of Phase 1 iterations in direction_only two-phase strategy for collinear (nspin=2)
```

- [ ] **Step 2: Add sc_dir_phase1_steps input item to read_input_item_other.cpp**

After the `sc_direction_only` item block (after line 199), add:

```cpp
    {
        Input_Item item("sc_dir_phase1_steps");
        item.annotation = "Phase 1 iterations for direction_only collinear strategy";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = R"(Number of SCF iterations for Phase 1 (magnitude constraint) in the direction_only two-phase strategy for collinear (nspin=2) calculations. During Phase 1, direction_only projection is temporarily disabled so BFGS can constrain the magnetic moment magnitude. After Phase 1, lambda decays and the system relaxes naturally. Minimum: 2.)";
        item.default_value = "5";
        item.unit = "";
        item.availability = "sc_mag_switch is true and sc_direction_only is true and nspin=2";
        read_sync_int(input.sc_dir_phase1_steps);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_dir_phase1_steps < 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_dir_phase1_steps must >= 2");
            }
        };
        this->add_item(item);
    }
```

- [ ] **Step 3: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 4: Commit**

```bash
git add source/source_io/module_parameter/input_parameter.h source/source_io/module_parameter/read_input_item_other.cpp
git commit -m "feat(deltaspin): add sc_dir_phase1_steps parameter

Controls Phase 1 duration in direction_only two-phase strategy
for collinear (nspin=2) calculations. Default: 5."
```

---

### Task 5: Rewrite mixing_restart Auto-Setting Logic

**Files:**
- Modify: `source/source_io/module_parameter/read_input_item_elec_stru.cpp:661-672`

- [ ] **Step 1: Replace the reset_value lambda**

In `source/source_io/module_parameter/read_input_item_elec_stru.cpp`, find the `mixing_restart` item's `reset_value` (lines 661-672). Replace the entire lambda:

Old:
```cpp
item.reset_value = [](const Input_Item& item, Parameter& para) {
    if (para.input.sc_mag_switch == 1)
    {// for DeltaSpin calculation, the mixing_restart should be same as sc_scf_thr
        if(para.input.sc_scf_thr != 10.0)
        {
            para.input.mixing_restart = para.input.sc_scf_thr;
        }
        else
        {// no mixing_restart until oscillation happen in PW base
            para.input.mixing_restart = para.input.scf_thr / 10.0;
        }
    }
};
```

New:
```cpp
item.reset_value = [](const Input_Item& item, Parameter& para) {
    if (para.input.sc_mag_switch)
    {
        if (para.input.sc_direction_only)
        {
            // direction_only path: Phase 2 uses explicit mix_reset(),
            // no auto mixing_restart needed (would interfere with Phase 1)
            para.input.mixing_restart = 0.0;
        }
        else if (para.input.sc_scf_thr_mode == "threshold")
        {
            // standard path: clean mixing history before lambda loop starts
            para.input.mixing_restart = para.input.sc_scf_thr;
        }
        else // "immediate"
        {
            // PW path: wait for oscillation before restarting
            para.input.mixing_restart = para.input.scf_thr / 10.0;
        }
    }
};
```

- [ ] **Step 2: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 3: Commit**

```bash
git add source/source_io/module_parameter/read_input_item_elec_stru.cpp
git commit -m "fix(deltaspin): rewrite mixing_restart auto-setting logic

Replace magic number 10.0 check with sc_scf_thr_mode.
direction_only path: mixing_restart=0 (Phase 2 uses explicit mix_reset).
threshold mode: mixing_restart=sc_scf_thr (standard gate).
immediate mode: mixing_restart=scf_thr/10 (PW oscillation detection)."
```

---

### Task 6: Restructure esolver_ks_lcao.cpp Control Flow

**Files:**
- Modify: `source/source_esolver/esolver_ks_lcao.cpp:438-574`

This is the most complex task. The existing block from line 438 to 574 needs to be restructured.

- [ ] **Step 1: Replace the entire DeltaSpin block in hamilt2rho_single**

Find the block starting at `bool skip_solve = false;` (line 438) and the following `if (PARAM.inp.sc_mag_switch)` block. Replace from `bool skip_solve = false;` through the end of the `if (PARAM.inp.sc_mag_switch)` block (up to but not including `// 3) run Hsolver`).

Replace with:

```cpp
    bool skip_solve = false;
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();

        if (PARAM.inp.sc_lambda_strategy == "linear_scan")
        {
            sc.run_lambda_linear_scan(iter - 1);
            skip_solve = true;
        }
        else if (PARAM.inp.sc_direction_only && PARAM.inp.nspin == 2)
        {
            // ================================================================
            // Collinear direction_only: two-phase strategy
            // ================================================================
            // For nspin=2, direction_only projection zeroes lambda entirely
            // (see lambda_loop.cpp). The two-phase strategy works around this:
            //
            // Phase 1 (iter 1..sc_dir_phase1_steps): BFGS with direction_only
            //   temporarily disabled, constraining moment MAGNITUDE to target.
            //   skip_solve=true: BFGS inner loop handles diagonalization.
            //
            // Phase 2 (iter > sc_dir_phase1_steps): Lambda decays gradually,
            //   normal SCF runs, system relaxes to magnetic ground state.
            // ================================================================
            if (iter <= PARAM.inp.sc_dir_phase1_steps)
            {
                sc.set_direction_only(false);
                sc.run_lambda_loop(iter - 1);
                sc.set_direction_only(true);
                skip_solve = true;
            }
            else
            {
                if (iter == PARAM.inp.sc_dir_phase1_steps + 1)
                {
                    // Reset mixing at Phase 1->2 transition.
                    // Phase 1 BFGS updates DM directly without charge mixing,
                    // so Broyden history is incompatible with Phase 2 SCF.
                    // Also reset mixing_restart_count and mixing_restart_step
                    // to avoid polluting mixing_dmr logic.
                    this->p_chgmix->mix_reset();
                    this->p_chgmix->mixing_restart_count = 0;
                    this->p_chgmix->mixing_restart_step = PARAM.inp.scf_nmax + 1;
                }

                // Gradual lambda decay: factor = 0.5^(1/3) per step
                // (~halves every 3 steps). Gradual decay avoids discontinuous
                // Hamiltonian change that would cause charge density oscillations.
                int nat = sc.get_nat();
                auto lambda = sc.get_sc_lambda();
                const double DECAY = std::pow(0.5, 1.0 / 3.0);
                for (int ia = 0; ia < nat; ++ia)
                    for (int ic = 0; ic < 3; ++ic)
                        lambda[ia][ic] *= DECAY;
                sc.set_lambda(lambda);
            }
        }
        else if (PARAM.inp.sc_direction_only && PARAM.inp.nspin == 4)
        {
            // Non-collinear direction_only: direction_only projection works
            // correctly for nspin=4 (only removes parallel component, leaving
            // transverse constraint). Use standard sc_scf_thr_mode gate.
            if (PARAM.inp.sc_scf_thr_mode == "immediate")
            {
                if (iter > 1)
                {
                    sc.run_lambda_loop(iter - 1);
                    if (!sc.mag_converged()) { sc.set_mag_converged(true); }
                    skip_solve = true;
                }
            }
            else
            {
                if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
                {
                    sc.run_lambda_loop(iter - 1);
                    sc.set_mag_converged(true);
                    skip_solve = true;
                }
                else if (sc.mag_converged())
                {
                    sc.run_lambda_loop(iter - 1);
                    skip_solve = true;
                }
            }
        }
        else
        {
            // Standard DeltaSpin (no direction_only)
            if (PARAM.inp.sc_scf_thr_mode == "immediate")
            {
                // "immediate" mode: activate lambda loop from iter>=2.
                // iter=1 is skipped because initial wavefunctions are not
                // available to compute initial magnetic moments.
                if (iter > 1)
                {
                    sc.run_lambda_loop(iter - 1);
                    if (!sc.mag_converged()) { sc.set_mag_converged(true); }
                    skip_solve = true;
                }
            }
            else
            {
                // "threshold" mode: activate when drho < sc_scf_thr.
                // drho > 0 excludes iter=1 where drho has not been computed yet.
                if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
                {
                    sc.run_lambda_loop(iter - 1);
                    sc.set_mag_converged(true);
                    skip_solve = true;
                }
                else if (sc.mag_converged())
                {
                    sc.run_lambda_loop(iter - 1);
                    skip_solve = true;
                }
            }
        }
    }
```

- [ ] **Step 2: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 3: Commit**

```bash
git add source/source_esolver/esolver_ks_lcao.cpp
git commit -m "refactor(deltaspin): restructure esolver control flow

- Separate nspin=2 and nspin=4 direction_only paths
- Use sc_scf_thr_mode instead of hardcoded sc_scf_thr=10.0
- Use sc_dir_phase1_steps instead of hardcoded 5
- Reset mixing_restart_count and mixing_restart_step at Phase transition
- Add drho>0 condition comment"
```

---

### Task 7: Update deltaspin_pw.cpp for sc_scf_thr_mode

**Files:**
- Modify: `source/source_pw/module_pwdft/deltaspin_pw.cpp`

- [ ] **Step 1: Update run_deltaspin_lambda_loop**

Replace the function body in `source/source_pw/module_pwdft/deltaspin_pw.cpp`. The current logic checks `drho < inp.sc_scf_thr` directly. Update to use `sc_scf_thr_mode`:

```cpp
bool run_deltaspin_lambda_loop(const int iter,
                               const double drho,
                               const Input_para& inp)
{
    if (!inp.sc_mag_switch)
    {
        return false;
    }

    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();

    if (inp.sc_lambda_strategy == "linear_scan")
    {
        sc.run_lambda_linear_scan(iter);
        return true;
    }

    if (inp.sc_scf_thr_mode == "immediate")
    {
        // "immediate" mode: activate from iter>=2.
        // iter=0/1 is skipped (no valid wavefunctions for initial Mi).
        if (iter >= 1) // iter is 0-indexed in PW (passed as iter-1 from caller)
        {
            sc.run_lambda_loop(iter);
            if (!sc.mag_converged()) { sc.set_mag_converged(true); }
            return true;
        }
        return false;
    }

    // "threshold" mode: activate when drho < sc_scf_thr
    // drho > 0 excludes iterations where drho has not been computed.
    if (!sc.mag_converged() && drho > 0 && drho < inp.sc_scf_thr)
    {
        sc.run_lambda_loop(iter);
        sc.set_mag_converged(true);
        return true;
    }
    else if (sc.mag_converged())
    {
        sc.run_lambda_loop(iter);
        return true;
    }

    return false;
}
```

**Note:** The PW esolver calls this with `iter - 1` (see `esolver_ks_pw.cpp:211`), so `iter >= 1` here corresponds to outer SCF iter >= 2. Verify by reading `esolver_ks_pw.cpp` to confirm the call convention.

- [ ] **Step 2: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 3: Commit**

```bash
git add source/source_pw/module_pwdft/deltaspin_pw.cpp
git commit -m "refactor(deltaspin): update PW path for sc_scf_thr_mode

Replace direct sc_scf_thr comparison with sc_scf_thr_mode logic.
immediate mode: activate from iter>=1 (outer iter>=2).
threshold mode: activate when drho < sc_scf_thr."
```

---

### Task 8: Update Comments in lambda_loop.cpp, init_sc.cpp, spin_constrain.h

**Files:**
- Modify: `source/source_lcao/module_deltaspin/lambda_loop.cpp`
- Modify: `source/source_lcao/module_deltaspin/init_sc.cpp`
- Modify: `source/source_lcao/module_deltaspin/spin_constrain.h`

- [ ] **Step 1: Update spin_constrain.h class documentation**

The current class doc comment (lines 26-55) references deleted strategies and outdated parameter descriptions. Replace the `@par` sections:

Find the existing `@par` sections and update them:

Replace `@par direction_only Mode and Two-Phase Strategy` with:

```
@par direction_only Mode and Two-Phase Strategy (collinear only)
The direction_only flag was designed for non-collinear calculations to constrain
only the spin DIRECTION (not magnitude). For collinear (nspin=2) mode,
the direction_only projection mathematically zeroes lambda (see lambda_loop.cpp).

The two-phase strategy in esolver_ks_lcao.cpp solves this for nspin=2:
  Phase 1 (iter 1..sc_dir_phase1_steps): BFGS with direction_only=FALSE constrains magnitude
  Phase 2 (iter > sc_dir_phase1_steps): Lambda decays to zero, system relaxes naturally

For nspin=4, direction_only works correctly and uses the standard sc_scf_thr_mode gate.
```

Replace `@par Parameter Recommendations` with:

```
@par Parameter Recommendations
- sc_scf_thr: Density error threshold for activating lambda loop. Default: 1e-3.
  Should be 10-100x larger than scf_thr.
- sc_scf_thr_mode: "threshold" (default, activate when drho<sc_scf_thr) or
  "immediate" (activate from iter>=2, for PW basis).
- mixing_restart: Auto-set based on sc_scf_thr_mode. See read_input_item_elec_stru.cpp.
- sc_dir_phase1_steps: Phase 1 duration for collinear direction_only. Default: 5.
- sc_direction_only: Set to true for direction-only constraint.
```

Remove any references to `HybridDelayed`, `AugmentedLagrangian`, `LinearResponse` strategies.

- [ ] **Step 2: Update init_sc.cpp comments**

In `source/source_lcao/module_deltaspin/init_sc.cpp`, find the comment block added at lines 82-85 and 96-98 (from recent commits). Update the `IMPORTANT` comment:

Old (around line 82-85):
```cpp
    // IMPORTANT: This masking is the root cause of why direction_only projection
    // zeroes lambda for nspin=2. Since only constrain.z is non-zero, the
    // direction_only projection (which removes the component parallel to target)
    // removes the entire lambda.z component, leaving lambda = 0.
```

New:
```cpp
    // This masking causes direction_only projection to zero lambda for nspin=2:
    // since only constrain.z is non-zero, the direction_only projection (which
    // removes the component parallel to target direction z) removes lambda.z entirely.
    // The two-phase strategy in esolver_ks_lcao.cpp handles this by temporarily
    // disabling direction_only during Phase 1 for collinear calculations.
```

Update the comment near line 96-98:

Old:
```cpp
    // NOTE: direction_only_ is designed for non-collinear (nspin=4) mode.
    // For nspin=2, it MUST be temporarily disabled during Phase 1 BFGS
    // (see esolver_ks_lcao.cpp) to allow magnitude constraint.
```

New:
```cpp
    // For nspin=2 with direction_only, the two-phase strategy in
    // esolver_ks_lcao.cpp temporarily disables direction_only during Phase 1
    // (controlled by sc_dir_phase1_steps). For nspin=4, direction_only
    // works correctly without two-phase.
```

- [ ] **Step 3: Update lambda_loop.cpp direction_only comments**

The direction_only comments in lambda_loop.cpp are already extensive but need minor updates to reference `sc_dir_phase1_steps` instead of "Phase 1 BFGS" as a fixed concept.

In the first direction_only comment block (around line 151-178), change the last line from:
```cpp
            // Therefore: direction_only MUST be disabled during Phase 1 BFGS for
            // collinear calculations. See esolver_ks_lcao.cpp for details.
```
to:
```cpp
            // Therefore: direction_only is temporarily disabled during Phase 1
            // (sc_dir_phase1_steps iterations) for collinear calculations.
            // See esolver_ks_lcao.cpp for details.
```

Similarly update the second direction_only comment block (around line 235-260), change:
```cpp
            //   This is another reason why direction_only must be disabled for
            //   collinear Phase 1 BFGS.
```
to:
```cpp
            //   This is another reason why direction_only is disabled during
            //   Phase 1 for collinear calculations.
```

And the third block (around line 383-390), change:
```cpp
            // This is why direction_only must be disabled for Phase 1.
```
to:
```cpp
            // This is why direction_only is disabled during Phase 1 for collinear calculations.
```

Same for the fourth block (around line 423-427):
```cpp
            // Same as above: after the optimal step correction, remove any parallel
            // component that may have been introduced.
            // For nspin=2: again zeroes the only non-zero component.
```

Add after:
```cpp
            // For nspin=2 collinear: this zeroes the only non-zero component (dnu.z).
```

- [ ] **Step 4: Build to verify compilation**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 5: Commit**

```bash
git add source/source_lcao/module_deltaspin/spin_constrain.h source/source_lcao/module_deltaspin/init_sc.cpp source/source_lcao/module_deltaspin/lambda_loop.cpp
git commit -m "docs(deltaspin): update comments for unified parameter scheme

- Remove references to deleted strategies (HybridDelayed, etc.)
- Reference sc_scf_thr_mode and sc_dir_phase1_steps in comments
- Clarify nspin=2 vs nspin=4 direction_only behavior"
```

---

### Task 9: Update Existing Input Parameter Tests

**Files:**
- Modify: `source/source_io/test/read_input_ptest.cpp`
- Modify: `source/source_io/test_serial/read_input_item_test.cpp`

- [ ] **Step 1: Update read_input_ptest.cpp**

Find the line that checks `sc_scf_thr` default (line 436):
```cpp
EXPECT_EQ(param.inp.sc_scf_thr, 1e-3);
```

This is already correct. Add checks for new parameters:
```cpp
EXPECT_EQ(param.inp.sc_scf_thr_mode, "threshold");
EXPECT_EQ(param.inp.sc_dir_phase1_steps, 5);
```

- [ ] **Step 2: Update read_input_item_test.cpp**

Find the `sc_scf_thr` test section (around line 1678) and verify it tests the correct default value. Add test sections for `sc_scf_thr_mode` and `sc_dir_phase1_steps`:

```cpp
{ // sc_scf_thr_mode
    auto it = find_label("sc_scf_thr_mode", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());
    param.input.sc_scf_thr_mode = "invalid";
    it->check_value(*it, param); // should WARNING_QUIT for invalid value
}
{ // sc_dir_phase1_steps
    auto it = find_label("sc_dir_phase1_steps", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());
    param.input.sc_dir_phase1_steps = 0;
    // check_value should reject < 2
}
```

Note: The exact test pattern depends on how `WARNING_QUIT` is handled in the test framework. Follow the existing pattern in the file.

- [ ] **Step 3: Build and run tests**

```bash
cd build && make -j4 2>&1 | tail -20
```

- [ ] **Step 4: Commit**

```bash
git add source/source_io/test/read_input_ptest.cpp source/source_io/test_serial/read_input_item_test.cpp
git commit -m "test(deltaspin): add tests for sc_scf_thr_mode and sc_dir_phase1_steps"
```

---

### Task 10: Full Build and Integration Test

**Files:** None (verification only)

- [ ] **Step 1: Full clean build**

```bash
cd build && cmake .. -DBUILD_TESTING=ON -DENABLE_LCAO=ON && make -j4 2>&1 | tail -30
```

Expected: clean build with no errors.

- [ ] **Step 2: Run unit tests for parameter reading**

```bash
cd build && ctest -R "read_input" -V 2>&1 | tail -30
```

Expected: all parameter tests pass.

- [ ] **Step 3: Verify the direction_only test case still works**

Run the existing DeltaSpin test case (e.g., case 24 or 26 from `tests/17_DS_DFTU/`) to verify the direction_only two-phase strategy still produces reasonable results.

- [ ] **Step 4: Final commit (if any fixes needed)**

If any fixes were needed during verification, commit them with appropriate messages.

# Weak Schur Spectral Preconditioner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a metric-compatible analytic spectral preconditioner to the fine-integral/coarse-response weak Schur solver, validate `1e-6` production accuracy, and continue the hBN adsorption campaign only after fixed-input gates pass.

**Architecture:** A new FFT preconditioner approximates the inverse analytic coarse kinetic shift in the physical regular-grid representation of the normalized Schur coordinates. `SternheimerWeakAugmented::Worker` exposes it directly to normalized-coordinate GMRES, while original-equation residuals remain the convergence authority.

**A/B correction:** The initial implementation wrapped the physical preconditioner in two complement square-root maps. The first remote hBN+H2O batch did not converge in 300 iterations, while the no-preconditioner control converged. A dedicated remote RED test proved that the physical callback must act directly in normalized coordinates; the complement square-root eigensolver separately passed a non-diagonal dense-reference test.

**Tech Stack:** C++17, FFTW3, BLAS/LAPACK, GoogleTest, ABACUS, Slurm on `df_iopcas_ghj`, LibRPA reader-v1 analysis.

---

### Task 1: Metric-Compatible Worker Hook

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_weak_augmented.h`
- Modify: `source/source_lcao/module_ri/sternheimer_weak_augmented.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_weak_augmented_test.cpp`

- [ ] **Step 1: Write failing tests for the complement square root and normalized preconditioner hook**

Add tests that verify `S_W^(1/2) S_W^(-1/2) x = x`, count preconditioner invocations during GMRES, compare the preconditioned solve with the independent dense solution, and reject wrong-sized or non-finite callback output.

- [ ] **Step 2: Run the remote test and verify RED**

Run `MODULE_RI_sternheimer_weak_augmented_test` in a clean remote build on `df_iopcas_ghj`. Expected: compilation fails because `apply_complement_sqrt` and the optional worker preconditioner API do not exist.

- [ ] **Step 3: Implement the low-rank square-root map and worker callback**

Expose `apply_complement_sqrt` for coordinate diagnostics and implement it as `S_W * S_W^(-1/2)` without a coarse square matrix. Accept an optional physical preconditioner callback in `Worker`; install it directly as the normalized right-preconditioner and validate each output. Do not add complement metric factors around the callback.

- [ ] **Step 4: Run the remote test and verify GREEN**

Run `MODULE_RI_sternheimer_weak_augmented_test`. Expected: all tests pass, including dense-solution equality and explicit invalid-output failures.

- [ ] **Step 5: Commit the worker hook**

Commit only the augmented worker and its tests with author `Codex <codex@openai.com>` and committer `AroundPeking <gonghuanjing@iphy.ac.cn>`.

### Task 2: Analytic Coarse-Grid Spectral Preconditioner

**Files:**
- Create: `source/source_lcao/module_ri/sternheimer_weak_preconditioner.h`
- Create: `source/source_lcao/module_ri/sternheimer_weak_preconditioner.cpp`
- Create: `source/source_lcao/module_ri/test/sternheimer_weak_preconditioner_test.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Write failing Fourier-mode tests**

Define tests for an orthogonal grid and a skew grid. Apply the requested preconditioner to individual Bloch Fourier modes and compare with `1/(alpha*|G+k|^2-epsilon+i*omega+regularization)`. Add signed-frequency conjugacy, alias-safe application, invalid grid, invalid shift, wrong-size, non-finite input, and singular-mode tests.

- [ ] **Step 2: Run the remote test and verify RED**

Run `MODULE_RI_sternheimer_weak_preconditioner_test`. Expected: compilation fails because the new class is absent.

- [ ] **Step 3: Implement the FFT preconditioner**

Create one FFTW forward/backward context per object, cache exact analytic inverse denominators, preserve the explicit Bloch phase, normalize the inverse FFT, and validate all public inputs and outputs. Do not use the finite-difference stencil symbol.

- [ ] **Step 4: Run focused and transfer regression tests**

Run `MODULE_RI_sternheimer_weak_preconditioner_test`, `MODULE_RI_sternheimer_grid_transfer_test`, and `MODULE_RI_sternheimer_fd_preconditioner_test`. Expected: all pass.

- [ ] **Step 5: Commit the analytic preconditioner**

Commit the class, build registration, and tests with the required attribution.

### Task 3: Weak-Q Runtime Integration

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] **Step 1: Write failing runtime-contract tests**

Add parser-level tests for default `spectral`, explicit `none`, default residual `1e-6`, explicit finite positive tolerance, and rejection of unknown modes or invalid tolerances. Verify the audit labels are stable.

- [ ] **Step 2: Run the remote smoke test and verify RED**

Run `MODULE_RI_sternheimer_abacus_st_smoke_test`. Expected: the new weak-q runtime helpers or labels are missing.

- [ ] **Step 3: Connect the preconditioner and tolerance**

Construct one analytic spectral preconditioner per weak worker and signed frequency, pass its physical regular-grid callback directly to `Worker`, and write mode, regularization, tolerance, iteration, Schur residual, and original residual to the audit. Keep `none` as a strict identity path.

- [ ] **Step 4: Run focused remote tests**

Run `MODULE_RI_sternheimer_abacus_st_smoke_test`, `MODULE_RI_sternheimer_weak_augmented_test`, `MODULE_RI_sternheimer_weak_preconditioner_test`, `MODULE_RI_sternheimer_weak_q_unit_test`, `MODULE_RI_sternheimer_weak_grid_test`, `MODULE_RI_sternheimer_delta_test`, and `MODULE_RI_sternheimer_periodic_solver_test`. Expected: all pass.

- [ ] **Step 5: Commit runtime integration**

Commit the integration and tests with the required attribution.

### Task 4: Remote Build and Fixed-Input Accuracy/Speed Gates

**Files:**
- Create remotely: isolated source archive, build directory, Slurm scripts, provenance files, and A/B result summaries under a new `df_iopcas_ghj` campaign directory
- Modify locally after evidence: `docs/superpowers/specs/2026-09-15-hbn-weak-schur-preconditioner-design.md`

- [ ] **Step 1: Verify source before upload**

Run `git diff --check`, inspect the three implementation commits and attribution, and create a source archive whose commit and SHA256 are recorded.

- [ ] **Step 2: Build once on `df_iopcas_ghj`**

Use an isolated remote build with `ENABLE_MPI=ON`, `ENABLE_LCAO=ON`, `ENABLE_LIBRI=ON`, `ENABLE_LIBCOMM=ON`, `DEBUG_INFO=ON`, and `ENABLE_MLALGO=OFF`. Record the binary SHA256 and CMake gates.

- [ ] **Step 3: Run the regression suite**

Run all focused tests from Task 3 plus the existing ABACUS-to-LibRPA smoke gates. Preserve the first failure and repair through a new unique build if required.

- [ ] **Step 4: Run hBN+H2O fixed-input A/B at `1e-8`**

Use the same physical q/frequency/unit inputs, binary, worker layout, and resources. Compare `none` with `spectral`; require finite matching reader-v1 metadata, all original residuals `<=1e-8`, matrix relative Frobenius difference `<=1e-8`, and reduced mean iterations and wall time.

- [ ] **Step 5: Run the `1e-6` tolerance comparison**

Compare `spectral,1e-6` with `spectral,1e-8`. Require every original residual `<=1e-6`, finite matching metadata, and a measured weighted single-q/frequency LibRPA contribution change `<=1 meV`.

- [ ] **Step 6: Record the measured gate**

Append exact job IDs, resources, commit, binary SHA256, equation counts, residual maxima, iteration distributions, wall time, MaxRSS, matrix errors, and RPA contribution difference to the design evidence. Do not label a partial response as a physical adsorption result.

### Task 5: Integrate and Resume the Adsorption Campaign

**Files:**
- Modify: `/Users/ghj/Downloads/同步空间/AITP_project/h2o_bn_rpa_adsorption/main.tex`
- Modify remotely: guarded remaining-q submission manifests and final LibRPA dependency chain

- [ ] **Step 1: Integrate the tested commits into `master_ghj`**

Push only after all gates pass. Verify remote `master_ghj`, author, committer, commit message, and binary provenance.

- [ ] **Step 2: Preserve current q001 evidence**

Verify its finalizer, equation completeness, residuals, finite reader-v1 fragments, and metadata. Do not recompute q001 solely to use the new solver.

- [ ] **Step 3: Submit missing representative q points once**

Use strict active/completed-output duplicate guards, the tested binary, `spectral` mode, residual `1e-6`, and resource layout informed by measured A/B wall time. Do not use `df_dcu` for execution.

- [ ] **Step 4: Complete LibRPA and adsorption subtraction**

Require full adsorbed and clean-system q/frequency coverage, finite matching reader-v1 data, accepted strict-2D head/wing handling, successful LibRPA outputs, and an independently checked adsorption-energy subtraction.

- [ ] **Step 5: Update and compile the research note**

Document implementation, producer, numerical-equivalence, tolerance, LibRPA, and physical-result gates separately. Compile the PDF and inspect the affected pages.

## Self-Review

- The plan covers the metric transform, analytic symbol, runtime switch, `1e-6` request, remote-only build, numerical/energy A/B, attribution, guarded production release, and final adsorption-energy gate.
- No production code is changed before a failing test is observed remotely.
- Existing q001 and unrelated local dirty worktrees are preserved.

# Sternheimer Operator Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the validated spectral preconditioner the no-configuration default and remove exact repeated work from the FD Hamiltonian and fixed-subspace projector without changing Delta-Sternheimer responses or RPA energies.

**Architecture:** Add one small runtime-option parser for strict default-on environment handling. Cache immutable FD stencil geometry inside each Hamiltonian, iterate only active mixed-derivative pairs, and add a reusable projector that caches fixed-vector norms while retaining the original sequential projection order.

**Tech Stack:** C++17, GoogleTest, CMake/CTest, OpenMP, Slurm on df_dcu, ABACUS Delta-Sternheimer, LibRPA.

---

### Task 1: Default-on spectral preconditioner

**Files:**
- Create: `source/source_lcao/module_ri/sternheimer_runtime_options.h`
- Create: `source/source_lcao/module_ri/sternheimer_runtime_options.cpp`
- Create: `source/source_lcao/module_ri/test/sternheimer_runtime_options_test.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.h`
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`

- [ ] **Step 1: Write failing runtime-option tests**

Add tests that unset a unique environment variable and expect the caller-supplied default, accept `true/1/on/yes` and `false/0/off/no` case-insensitively, and reject empty or unknown nonempty values. Add a solver-options test expecting the spectral preconditioner to default to enabled with zero regularization.

```cpp
TEST(SternheimerRuntimeOptions, UsesDefaultWhenEnvironmentIsMissing)
{
    unsetenv("ABACUS_STERNHEIMER_TEST_FLAG");
    EXPECT_TRUE(ModuleRI::sternheimer_environment_flag("ABACUS_STERNHEIMER_TEST_FLAG", true));
    EXPECT_FALSE(ModuleRI::sternheimer_environment_flag("ABACUS_STERNHEIMER_TEST_FLAG", false));
}

TEST(SternheimerRuntimeOptions, RejectsInvalidBooleanText)
{
    setenv("ABACUS_STERNHEIMER_TEST_FLAG", "maybe", 1);
    EXPECT_THROW(ModuleRI::sternheimer_environment_flag("ABACUS_STERNHEIMER_TEST_FLAG", true),
                 std::invalid_argument);
    unsetenv("ABACUS_STERNHEIMER_TEST_FLAG");
}

TEST(SternheimerRPA, SpectralPreconditionerIsDefault)
{
    const ModuleRI::SternheimerRPA::SolverOptions options;
    EXPECT_TRUE(options.use_fd_spectral_preconditioner);
    EXPECT_DOUBLE_EQ(options.fd_spectral_preconditioner_regularization, 0.0);
}
```

- [ ] **Step 2: Build the two focused tests and verify RED**

Run:

```bash
cmake --build build --target MODULE_RI_sternheimer_runtime_options_test MODULE_RI_sternheimer_rpa_test -j8
```

Expected: compilation fails because `sternheimer_environment_flag` does not exist and the solver-options default assertions fail after the parser test is made buildable.

- [ ] **Step 3: Implement strict default-on behavior**

Declare and define:

```cpp
bool sternheimer_environment_flag(const char* name, bool default_value);
```

Missing variables return `default_value`; recognized true and false strings return the corresponding value; empty or unknown values throw `std::invalid_argument` naming the variable. Change `SolverOptions` defaults to `true` and `0.0`. Replace the three production assignments with:

```cpp
solver_options.use_fd_spectral_preconditioner
    = sternheimer_environment_flag(kSpectralPreconditionerEnv, true);
solver_options.fd_spectral_preconditioner_regularization
    = nonnegative_double_from_env(kSpectralPreconditionerRegularizationEnv, 0.0);
```

- [ ] **Step 4: Verify GREEN and neighboring tests**

Run:

```bash
ctest --test-dir build -R '^(MODULE_RI_sternheimer_runtime_options_test|MODULE_RI_sternheimer_rpa_test|MODULE_RI_sternheimer_fd_preconditioner_test|MODULE_RI_sternheimer_periodic_solver_test|MODULE_RI_sternheimer_abacus_st_smoke_test)$' --output-on-failure
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit Task 1**

```bash
git add source/source_lcao/module_ri
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'perf(sternheimer): enable spectral preconditioner by default'
```

### Task 2: Cache FD stencil geometry and skip exact zero work

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp`

- [ ] **Step 1: Write failing stencil-plan tests**

Add public read-only diagnostic accessors and tests proving that an explicit orthogonal cell has zero active mixed pairs and zero cached first-derivative directions, while a skew cell retains its nonzero pairs. Also compare FD8 kinetic output against the existing manual/reference stencil for periodic Gamma, periodic twisted, and nonperiodic grids.

```cpp
EXPECT_EQ(orthogonal.active_mixed_derivative_pair_count(), 0);
EXPECT_EQ(orthogonal.cached_first_derivative_direction_count(), 0);
EXPECT_GT(skew.active_mixed_derivative_pair_count(), 0);
```

- [ ] **Step 2: Build and verify RED**

Run:

```bash
cmake --build build --target MODULE_RI_sternheimer_fd_hamiltonian_test -j8
```

Expected: compilation fails because the cache diagnostics do not exist.

- [ ] **Step 3: Implement immutable constructor cache**

Store FD radius/weights, Laplacian coefficients, center coefficient, active `(left,right)` pairs, required right-derivative directions, and a Gamma flag in the Hamiltonian. Initialize the cache after constructor validation. Define an active pair only by exact comparison:

```cpp
if (laplacian_coefficients_[left][right] != 0.0)
{
    active_mixed_derivative_pairs_[active_mixed_derivative_pair_count_++] = {left, right};
    required_first_derivative_directions_[right] = true;
}
```

Use the cache in `apply_grid_terms`. Iterate only active pairs and allocate first-derivative buffers only for required directions. Preserve the existing arithmetic order for each retained term.

- [ ] **Step 4: Add exact phase shortcuts**

In `shifted_grid_point`, return `1+0i` when the lattice translation is zero or when the cached k point is Gamma. Continue calling `sternheimer_bloch_phase` for non-Gamma boundary translations.

- [ ] **Step 5: Verify GREEN and the FD solver neighborhood**

Run:

```bash
ctest --test-dir build -R '^(MODULE_RI_sternheimer_fd_hamiltonian_test|MODULE_RI_sternheimer_fd_preconditioner_test|MODULE_RI_sternheimer_fd_solver_test|MODULE_RI_sternheimer_delta_test|MODULE_RI_sternheimer_periodic_solver_test)$' --output-on-failure
```

Expected: all selected tests pass with numerical tolerances unchanged.

- [ ] **Step 6: Commit Task 2**

```bash
git add source/source_lcao/module_ri/sternheimer_fd_hamiltonian.* \
        source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'perf(sternheimer): cache finite-difference operator geometry'
```

### Task 3: Cache fixed-subspace projection norms

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.h`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_delta.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_delta_test.cpp`

- [ ] **Step 1: Write a failing dot-count test**

Construct `SternheimerSubspaceProjector` with two non-unit vectors and a dot callback that increments a counter. Verify two norm dots occur at construction and exactly one dot per basis vector occurs for each later projection. Verify zero-norm vectors are skipped and cached projection matches `project_out_subspace` to `1e-13`.

```cpp
int dot_calls = 0;
auto counted_dot = [&dot_calls](const Vector& lhs, const Vector& rhs) {
    ++dot_calls;
    return ModuleRI::SternheimerRPA::local_grid_dot(lhs, rhs, 1.0);
};
ModuleRI::SternheimerSubspaceProjector projector(subspace, counted_dot);
EXPECT_EQ(dot_calls, 2);
projector.project(vec);
EXPECT_EQ(dot_calls, 4);
```

- [ ] **Step 2: Build and verify RED**

Run:

```bash
cmake --build build --target MODULE_RI_sternheimer_rpa_test -j8
```

Expected: compilation fails because `SternheimerSubspaceProjector` does not exist.

- [ ] **Step 3: Implement the cached projector**

Add a class storing a reference to the immutable subspace, the dot callback, and precomputed norms. Keep sequential projection and the current OpenMP `axpy` loop. Make `SternheimerRPA::project_out_subspace` construct a temporary projector for compatibility.

- [ ] **Step 4: Use one cached projector per Delta linear solve**

Create `std::shared_ptr<const SternheimerSubspaceProjector>` after the fixed subspace and dot callback are available. Use it for projected RHS, operator input/output, preconditioner output, and final response projection. Do not alter virtual-residual low-rank correction order.

- [ ] **Step 5: Verify GREEN and solver tests**

Run:

```bash
ctest --test-dir build -R '^(MODULE_RI_sternheimer_rpa_test|MODULE_RI_sternheimer_delta_test|MODULE_RI_sternheimer_fd_solver_test|MODULE_RI_sternheimer_periodic_solver_test)$' --output-on-failure
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit Task 3**

```bash
git add source/source_lcao/module_ri/sternheimer_rpa.* \
        source/source_lcao/module_ri/sternheimer_delta.cpp \
        source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp \
        source/source_lcao/module_ri/test/sternheimer_delta_test.cpp
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'perf(sternheimer): cache fixed-subspace projection norms'
```

### Task 4: Full local regression and source audit

**Files:**
- Modify only if a test exposes a defect in Tasks 1-3.

- [ ] **Step 1: Build all directly affected targets**

```bash
cmake --build build --target \
  MODULE_RI_sternheimer_runtime_options_test \
  MODULE_RI_sternheimer_rpa_test \
  MODULE_RI_sternheimer_fd_hamiltonian_test \
  MODULE_RI_sternheimer_fd_preconditioner_test \
  MODULE_RI_sternheimer_fd_solver_test \
  MODULE_RI_sternheimer_delta_test \
  MODULE_RI_sternheimer_periodic_solver_test \
  MODULE_RI_sternheimer_abacus_st_smoke_test -j8
```

- [ ] **Step 2: Run all directly affected tests**

```bash
ctest --test-dir build -R '^MODULE_RI_sternheimer_(runtime_options|rpa|fd_hamiltonian|fd_preconditioner|fd_solver|delta|periodic_solver|abacus_st_smoke)_test$' --output-on-failure
```

Expected: eight tests pass.

- [ ] **Step 3: Audit defaults and exact-zero semantics**

```bash
rg -n 'use_fd_spectral_preconditioner = true|regularization = 0.0|!= 0.0|sternheimer_environment_flag' source/source_lcao/module_ri
git diff --check origin/master_ghj...HEAD
```

Expected: one default definition, three production default-on assignments, exact-zero mixed-pair selection, and no whitespace errors.

### Task 5: Build and validate on df_dcu

**Files:**
- Create: `server_jobs/sternheimer_performance_20260830/build_operator_acceleration_dfdcu.slurm`
- Create: `server_jobs/sternheimer_performance_20260830/benchmark_hf_operator_acceleration_dfdcu.slurm`
- Create: `server_jobs/sternheimer_performance_20260830/benchmark_si_q7_operator_acceleration_dfdcu.slurm`
- Create: `server_jobs/sternheimer_performance_20260830/compare_sternheimer_v1_matrices.py`

- [ ] **Step 1: Package the exact feature commit and stage immutable inputs**

Create a git archive, SHA256 manifest, and scripts under local `server_jobs/sternheimer_performance_20260830`. Transfer them through the existing df_dcu ControlMaster connection to `/work1/ghj/sternheimer-performance-20260830`.

- [ ] **Step 2: Submit one normal-partition build**

Use one node, one MPI task, 30 build threads, `RelWithDebInfo`, LibRI/LibComm, GreenX, and `BUILD_TESTING=ON`. Run the eight tests from Task 4, build `abacus_3p`, and create one immutable artifact with source, executable, dependency, test, and SHA provenance.

- [ ] **Step 3: Run HF precision and timing regression**

Reuse `/work1/ghj/sternheimer_abacus_tests/molecular_convergence_20260818/response_frequency/job21654986/hf_mol_90`, four nodes, one MPI rank and 30 OpenMP threads per node, one accepted frequency, FD8, tolerance `1e-6`, and no explicit preconditioner environment variable. Compare against the accepted response matrix with threshold `1e-8`.

- [ ] **Step 4: Run Si q=7 precision and timing regression**

Reuse `/work1/ghj/si-solid-grid-fd8-final-20260823/materials/si/solid`, grid `30^3`, q index 7, accepted frequency 4, one occupied-band timing limit, four nodes, one MPI rank and 30 OpenMP threads per node, FD8, tolerance `1e-8`, and no explicit preconditioner environment variable. Compare with the accepted preconditioned matrix at threshold `1e-8`.

- [ ] **Step 5: Record incremental speedups**

For HF and Si record application status, convergence, equations, iterations, response differences, maximum residuals, total process wall, maximum RSS, nodes, ranks, threads, and node-hours. Compare against the accepted `4bee9e3b1` spectral-preconditioner timings.

### Task 6: Commit, integrate, and document accepted acceleration

**Files:**
- Create: `server_jobs/sternheimer_performance_20260830/STERNHEIMER_OPERATOR_ACCELERATION_BENCHMARK_20260830.md`
- Modify: `/Users/ghj/同步空间/AITP_project/delta_st_rpa_project/development_notes/sections/sternheimer_fd_spectral_preconditioner.tex`
- Modify: `/Users/ghj/.codex/skills/abacus-delta-st-efficient-production/SKILL.md`

- [ ] **Step 1: Record exact implementation and A/B evidence**

Document the mathematical equivalence, test results, feature commit, executable hash, HF/Si matrix differences, residual gates, timing, resources, and limits of the benchmark.

- [ ] **Step 2: Review and merge the feature branch**

Run final spec and code-quality review, verify commit attribution, merge the feature branch into local `master_ghj` with a merge commit using Codex author and AroundPeking committer, and push `master_ghj` only after confirming the remote head has not moved unexpectedly.

- [ ] **Step 3: Rebuild the merged commit if its tree differs**

If the merge commit tree differs from the validated feature artifact, build and rerun focused tests. Otherwise record tree identity and retain the accepted executable artifact.

- [ ] **Step 4: Update and verify documentation**

Update the development TeX note and production skill. Run `latexmk -xelatex`, render the changed pages, inspect them, and record the final PDF page count and SHA256.

### Task 7: Final full-material benchmark

**Files:**
- Create: `server_jobs/sternheimer_performance_20260830/run_si_fd8_full_operator_acceleration_dfdcu.slurm`
- Append: `server_jobs/sternheimer_performance_20260830/STERNHEIMER_OPERATOR_ACCELERATION_BENCHMARK_20260830.md`
- Append: `/Users/ghj/同步空间/AITP_project/delta_st_rpa_project/development_notes/sections/sternheimer_fd_spectral_preconditioner.tex`

- [ ] **Step 1: Audit the accepted Si FD8 baseline**

Confirm the accepted grid, all canonical q indices, 4x4x4 k/q mesh, 12 frequencies, PCA `1e-6`, full Coulomb, analytic q-average head/wing with body start 1, PP/NAO/ABFS hashes, and prior response/LibRPA success. Do not rerun any baseline physics solely for timing.

- [ ] **Step 2: Submit a duplicate-safe optimized response array**

Use only the final validated merged executable, `normal` partition, `/work1`, the accepted one-rank-per-node and 30-thread layout, and a wall limit below 24 hours. Each array task gets a unique case directory and checks for an existing success marker before submission.

- [ ] **Step 3: Validate every q and run LibRPA head/wing postprocessing**

Require scheduler success, application `status success`, `all_converged yes`, complete frequency/response counts, partial and symmetry manifests, matching provenance, and finite output. Reuse the validated LibRPA head/wing artifact; do not rerun ABACUS for postprocessing failures.

- [ ] **Step 4: Report final end-to-end acceleration**

Report solid total energy and comparison with the accepted value, producer wall time by q and aggregate, LibRPA wall time, critical-path time, node-hours, peak RSS, node/rank/thread layout, and speedup against the same-input accepted baseline. Distinguish first-order solver speedup from full ABACUS and complete ABACUS-plus-LibRPA speedup.

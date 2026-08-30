# RPA On-Site Coulomb Preorthogonalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in ABACUS RPA auxiliary-basis on-site Coulomb preorthogonalization switch, use the same transformed basis in reader-v1 producer and Sternheimer response paths, and make it the validated default in the molecular Delta-ST skill.

**Architecture:** A focused `rpa_abfs_preorthogonalization` module owns the radial Coulomb metric, repeated modified Gram-Schmidt, channel diagnostics, and orbital reconstruction. RPA producer and Sternheimer code call that one module only after their existing product-PCA or explicit-ABFS construction; ordinary EXX/GW construction remains unchanged. ABACUS defaults to the legacy `none` mode, while the installed molecular skill explicitly selects `onsite_coulomb`, threshold `1e-2`, LibRPA hard cutoff `1e-4`, and the already validated spectral Sternheimer preconditioner.

**Tech Stack:** C++17, ABACUS `Numerical_Orbital_Lm`, GoogleTest/GoogleMock, CMake/CTest, reader-v1 binary matrix comparison, Markdown Codex skill files, Slurm validation on `df_iopcas_ghj`.

---

## File Map

- `source/source_io/module_parameter/input_parameter.h`: stores the two new runtime inputs and backward-compatible defaults.
- `source/source_io/module_parameter/read_input_item_exx_dftu.cpp`: parses and validates mode and threshold.
- `source/source_io/test_serial/read_input_item_test.cpp`: proves parsing, defaults, and invalid-value failures.
- `source/source_lcao/module_ri/rpa_abfs_preorthogonalization.h`: public types and the single transformation API shared by producer and response.
- `source/source_lcao/module_ri/rpa_abfs_preorthogonalization.cpp`: Coulomb metric, repeated MGS, reconstruction, diagnostics, and report formatting.
- `source/source_lcao/module_ri/test/rpa_abfs_preorthogonalization_test.cpp`: numerical, legacy no-op, and provenance tests.
- `source/source_lcao/module_ri/CMakeLists.txt`: builds the shared implementation.
- `source/source_lcao/module_ri/test/CMakeLists.txt`: registers the focused test target.
- `source/source_lcao/module_ri/RPA_LRI.hpp`: applies the shared transform only in RPA producer ABFS preparation.
- `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: applies the same transform to Sternheimer ABFS and writes matching provenance.
- `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`: verifies Sternheimer provenance formatting.
- `docs/advanced/input_files/input-main.md`: documents both public input keys and their scope.
- `/Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence/SKILL.md`: makes the validated stable route the default for new molecular work.
- `/Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence/references/parameter-contract.md`: gives exact producer/response and LibRPA settings.

### Task 1: Parse And Validate The Runtime Switch

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_exx_dftu.cpp`
- Test: `source/source_io/test_serial/read_input_item_test.cpp`

- [ ] **Step 1: Write the failing input tests**

Add assertions beside `rpa_ccp_rmesh_times` that exercise the intended public contract:

```cpp
EXPECT_EQ(param.input.rpa_abfs_preorth, "none");
EXPECT_DOUBLE_EQ(param.input.rpa_abfs_preorth_threshold, 1.0e-2);

auto mode = find_label("rpa_abfs_preorth", readinput.input_lists);
ASSERT_NE(mode, readinput.input_lists.end());
mode->second.str_values = {"onsite_coulomb"};
mode->second.read_value(mode->second, param);
EXPECT_EQ(param.input.rpa_abfs_preorth, "onsite_coulomb");

param.input.rpa_abfs_preorth = "invalid";
EXPECT_EXIT(mode->second.check_value(mode->second, param),
            ::testing::ExitedWithCode(1), "");

auto threshold = find_label("rpa_abfs_preorth_threshold", readinput.input_lists);
ASSERT_NE(threshold, readinput.input_lists.end());
for (const double invalid : {0.0, 1.0, -1.0,
                             std::numeric_limits<double>::infinity()})
{
    param.input.rpa_abfs_preorth_threshold = invalid;
    EXPECT_EXIT(threshold->second.check_value(threshold->second, param),
                ::testing::ExitedWithCode(1), "");
}
```

- [ ] **Step 2: Build and run the input test to verify RED**

Run on the configured ABACUS build host:

```bash
cmake --build "$build" --parallel 40 --target MODULE_IO_read_item_serial
ctest --test-dir "$build" --output-on-failure -R '^MODULE_IO_read_item_serial$'
```

Expected: compilation fails because the two fields do not exist.

- [ ] **Step 3: Add the minimal input implementation**

Add these defaults to `Input_para`:

```cpp
std::string rpa_abfs_preorth = "none";
double rpa_abfs_preorth_threshold = 1.0e-2;
```

Register both keys in `read_input_item_exx_dftu.cpp`. Accept exactly `none` and `onsite_coulomb`; require a finite threshold strictly inside `(0, 1)` and issue `ModuleBase::WARNING_QUIT` with the rejected key and allowed range.

- [ ] **Step 4: Rebuild and verify GREEN**

Run the same build and CTest commands. Expected: `MODULE_IO_read_item_serial` passes.

- [ ] **Step 5: Commit the input contract**

```bash
git add source/source_io/module_parameter/input_parameter.h \
        source/source_io/module_parameter/read_input_item_exx_dftu.cpp \
        source/source_io/test_serial/read_input_item_test.cpp
git commit -m 'feat(rpa): add auxiliary preorthogonalization inputs'
```

### Task 2: Implement The Shared Coulomb-Metric Transformation

**Files:**
- Create: `source/source_lcao/module_ri/rpa_abfs_preorthogonalization.h`
- Create: `source/source_lcao/module_ri/rpa_abfs_preorthogonalization.cpp`
- Create: `source/source_lcao/module_ri/test/rpa_abfs_preorthogonalization_test.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Define the tested public API in the new test**

The header must expose one shared orbital-set type and one mutating entry point:

```cpp
namespace ModuleRI
{
using RpaAbfsOrbitalSet
    = std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>;

struct RpaAbfsPreorthChannelReport
{
    int atom_type = -1;
    int angular_momentum = -1;
    int input_count = 0;
    int output_count = 0;
    std::vector<int> rejected_indices;
    double minimum_residual_norm = 1.0;
    double maximum_identity_error = 0.0;
};

struct RpaAbfsPreorthReport
{
    std::string mode = "none";
    double threshold = 1.0e-2;
    int input_expanded_size = 0;
    int output_expanded_size = 0;
    std::vector<RpaAbfsPreorthChannelReport> channels;
};

RpaAbfsPreorthReport apply_rpa_abfs_preorthogonalization(
    RpaAbfsOrbitalSet& abfs,
    const std::string& mode,
    double threshold,
    bool force_flag);

std::string format_rpa_abfs_preorth_report(
    const RpaAbfsPreorthReport& report);
}
```

Construct synthetic odd uniform radial meshes in the test and prove separately that:

1. `none` leaves count, `psi`, and metadata byte-for-byte unchanged;
2. two identical `l=0` radials retain one radial at threshold `1e-2`;
3. two independent radials remain and their recomputed Coulomb metric differs from identity by at most `1e-10`;
4. rejecting one `l=2` radial removes exactly five expanded functions;
5. zero/negative norms and invalid modes throw `std::invalid_argument`;
6. the report contains input/output counts, rejected indices, minimum residual norm, and identity error.

- [ ] **Step 2: Register the test and verify RED**

Add a `MODULE_RI_rpa_abfs_preorthogonalization_test` target with `parameter`, `base`, `${math_libs}`, and the new `.cpp`. Build and run:

```bash
cmake --build "$build" --parallel 40 \
  --target MODULE_RI_rpa_abfs_preorthogonalization_test
ctest --test-dir "$build" --output-on-failure \
  -R '^MODULE_RI_rpa_abfs_preorthogonalization_test$'
```

Expected: compilation fails because the shared header and implementation are absent.

- [ ] **Step 3: Implement the discrete on-site Coulomb metric**

For every `(type,l)` channel, build Simpson quadrature weights from the channel's odd radial mesh and `rab`. Exclude `r=0` from singular powers. Evaluate the separable discrete form

```text
sum_i w_i r_i^2 f_a(i) [
  r_i^(-l-1) sum_(j<=i) w_j r_j^(l+2) f_b(j)
  + r_i^l sum_(j>i) w_j r_j^(1-l) f_b(j)
]
```

and multiply by `4*pi/(2*l+1)`. Use prefix/suffix sums so cost is `O(nradial^2 * nmesh)`, symmetrize the matrix, and reject inconsistent channel grids before orthogonalization.

- [ ] **Step 4: Implement repeated modified Gram-Schmidt**

Traverse input radials in their original order. Coulomb-normalize each candidate, compute its squared residual against retained vectors, reject when `residual_norm <= threshold`, and otherwise reorthogonalize for at most four passes until summed squared overlap is at most `1e-20`. Throw on non-finite/non-positive norms. Rebuild each output `Numerical_Orbital_Lm` with sequential `chi` indices and the Coulomb-unit radial values; do not L2-normalize after reconstruction.

- [ ] **Step 5: Verify GREEN and numerical diagnostics**

Run the focused test. Expected: all six behaviors pass and the identity error is no larger than `1e-10`.

- [ ] **Step 6: Commit the numerical core**

```bash
git add source/source_lcao/module_ri/rpa_abfs_preorthogonalization.* \
        source/source_lcao/module_ri/test/rpa_abfs_preorthogonalization_test.cpp \
        source/source_lcao/module_ri/CMakeLists.txt \
        source/source_lcao/module_ri/test/CMakeLists.txt
git commit -m 'feat(rpa): preorthogonalize ABFS in onsite Coulomb metric'
```

### Task 3: Apply The Shared Transform Only To RPA Producer Bases

**Files:**
- Modify: `source/source_lcao/module_ri/RPA_LRI.hpp`
- Test: `source/source_lcao/module_ri/test/rpa_abfs_preorthogonalization_test.cpp`

- [ ] **Step 1: Add a failing producer-finalization test**

Add a small helper to the shared module:

```cpp
RpaAbfsPreorthReport finalize_rpa_abfs_from_input(
    RpaAbfsOrbitalSet& abfs,
    const Input_para& input,
    bool force_flag);
```

The test supplies `Input_para` first with `none`, then `onsite_coulomb`, and verifies the same synthetic basis is unchanged or reduced respectively. Expected RED: missing helper.

- [ ] **Step 2: Implement the input adapter and verify GREEN**

The adapter delegates directly to `apply_rpa_abfs_preorthogonalization` using only `input.rpa_abfs_preorth` and `input.rpa_abfs_preorth_threshold`.

- [ ] **Step 3: Call the adapter at every RPA-only preparation site**

In `RPA_LRI.hpp`, call the adapter immediately after `ExxLriDetail::prepare_abfs` for:

- the normal `abfs` branch in `cal_postSCF_exx`;
- the optional `abfs_shrink` branch in `cal_postSCF_exx`;
- the `cal_large_Cs` basis;
- any additional RPA-only `prepare_abfs` call found by
  `rg -n 'prepare_abfs\(' source/source_lcao/module_ri/RPA_LRI.hpp`.

Write the formatted report to `GlobalV::ofs_running` on rank zero. Do not modify `ExxLriDetail::prepare_abfs` or any `Exx_LRI::init*` path, because those are shared with ordinary EXX/GW.

- [ ] **Step 4: Build the RPA objects and focused test**

```bash
cmake --build "$build" --parallel 40 \
  --target ri MODULE_RI_rpa_abfs_preorthogonalization_test
ctest --test-dir "$build" --output-on-failure \
  -R '^MODULE_RI_rpa_abfs_preorthogonalization_test$'
```

Expected: compile and test pass; `rg` shows each RPA preparation followed by the shared adapter and no new call under ordinary `Exx_LRI` initialization.

- [ ] **Step 5: Commit producer integration**

```bash
git add source/source_lcao/module_ri/RPA_LRI.hpp \
        source/source_lcao/module_ri/rpa_abfs_preorthogonalization.* \
        source/source_lcao/module_ri/test/rpa_abfs_preorthogonalization_test.cpp
git commit -m 'feat(rpa): use preorthogonalized basis in producer output'
```

### Task 4: Apply The Same Transform To Sternheimer Response Bases

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] **Step 1: Add a failing Sternheimer provenance test**

Extend the existing report test so a result carrying a preorth report must emit:

```text
rpa_abfs_preorth onsite_coulomb
rpa_abfs_preorth_threshold 0.01
rpa_abfs_preorth_expanded_size <input> <output>
rpa_abfs_preorth_channel <type> <l> <input> <output> ...
```

Expected RED: the result/report fields do not exist.

- [ ] **Step 2: Route both Sternheimer ABFS builders through the shared adapter**

After product-PCA generation or explicit `ABFS_ORBITAL` loading and `filter_empty_orbs`, invoke `finalize_rpa_abfs_from_input` in both `build_sternheimer_abfs` and `build_abfs_ccp_data`. Preserve the returned report in `SternheimerABFBuildData` and in the final status/report object. Do not transform `abfs_ccp` independently; it must be generated from the already transformed ABFS.

- [ ] **Step 3: Write complete Sternheimer provenance**

Append the shared formatted report to `STERNHEIMER_CHI0.dat` and the running log. Keep the existing `basis_aux_out` and reader-v1 dimensions authoritative. The `none` report must explicitly state `none` but must not change any matrix or dimension.

- [ ] **Step 4: Build and verify the Sternheimer tests**

```bash
cmake --build "$build" --parallel 40 --target \
  MODULE_RI_rpa_abfs_preorthogonalization_test \
  MODULE_RI_sternheimer_abacus_st_smoke_test
ctest --test-dir "$build" --output-on-failure -R \
  '^(MODULE_RI_rpa_abfs_preorthogonalization_test|MODULE_RI_sternheimer_abacus_st_smoke_test)$'
```

Expected: both pass; report text includes mode, threshold, dimensions, and per-channel diagnostics.

- [ ] **Step 5: Commit Sternheimer integration**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp \
        source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp
git commit -m 'feat(sternheimer): share RPA ABFS preorthogonalization'
```

### Task 5: Document The Switch And Update The Molecular Skill

**Files:**
- Modify: `docs/advanced/input_files/input-main.md`
- Modify: `/Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence/SKILL.md`
- Modify: `/Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence/references/parameter-contract.md`

- [ ] **Step 1: Record the pre-edit skill behavior**

Confirm the installed skill currently lacks the new keys and still recommends the old PCA ladder:

```bash
rg -n 'rpa_abfs_preorth|sqrt_coulomb_threshold|fd_spectral|1e-4/1e-6/1e-8' \
  /Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence
```

Expected RED evidence: no `rpa_abfs_preorth` or `fd_spectral` production contract, and the old ladder remains.

- [ ] **Step 2: Document the ABACUS input keys**

Add both entries to the Exact Exchange/RPA table of contents and parameter sections. State:

- global default `rpa_abfs_preorth none`;
- opt-in `onsite_coulomb` affects only RPA producer and Sternheimer LibRPA auxiliary bases;
- threshold must lie inside `(0,1)` and defaults to `1e-2`;
- ordinary EXX/GW remains unchanged;
- the global LibRPA Coulomb cutoff is a separate input.

- [ ] **Step 3: Update the molecular skill's stable production contract**

Add these mandatory settings for new molecular producer and response inputs:

```text
rpa_abfs_preorth              onsite_coulomb
rpa_abfs_preorth_threshold    1e-2
```

Set the auxiliary ladder to `exx_pca_threshold 1e-4, 1e-5, 1e-6`, require p5-to-p6 changes below `0.1 kcal/mol` independently for every absolute molecule/atom `Ec` and `Dc`, and use `sqrt_coulomb_threshold = 1e-4` as the primary LibRPA definition with `1e-5/1e-6/0` archived. State that `none` is for historical reproduction or explicit A/B diagnostics only.

Require the status fields:

```text
rpa_abfs_preorth=onsite_coulomb
rpa_abfs_preorth_threshold=1e-2
sternheimer_preconditioner=fd_spectral
sternheimer_preconditioner_regularization_Ry=0
```

Clarify that the algorithm is general but automatic skill use remains molecule-only until a periodic-solid regression passes.

- [ ] **Step 4: Validate the skill and inspect the diff**

```bash
python /Users/ghj/.codex/skills/.system/skill-creator/scripts/quick_validate.py \
  /Users/ghj/.codex/skills/abacus-delta-st-molecular-convergence
git diff --check
```

Expected: validator succeeds, the parameter contract contains the exact keys, and no unrelated skill text changes.

- [ ] **Step 5: Commit repository documentation**

The installed skill is outside the ABACUS repository and is recorded separately in the final report. Commit the repository input documentation:

```bash
git add docs/advanced/input_files/input-main.md
git commit -m 'docs(rpa): document onsite Coulomb preorthogonalization'
```

### Task 6: Run Regression Gates, Commit Attribution Checks, And Push

**Files:**
- Verify: all files changed in Tasks 1-5
- Produce: remote N2 p5 one-frequency regression artifacts under `/data/home/df_iopcas_ghj/sternheimer_tests/n2_rpa_preorth_runtime_regression_20260830`

- [ ] **Step 1: Run the focused local/remote unit suite**

Build and require:

```bash
cmake --build "$build" --parallel 40 --target \
  MODULE_IO_read_item_serial \
  MODULE_RI_rpa_abfs_preorthogonalization_test \
  MODULE_RI_sternheimer_abacus_st_smoke_test \
  MODULE_RI_sternheimer_fd_preconditioner_test \
  MODULE_RI_sternheimer_fd_solver_test \
  MODULE_RI_sternheimer_delta_test
ctest --test-dir "$build" --output-on-failure -R \
  '^(MODULE_IO_read_item_serial|MODULE_RI_rpa_abfs_preorthogonalization_test|MODULE_RI_sternheimer_abacus_st_smoke_test|MODULE_RI_sternheimer_fd_preconditioner_test|MODULE_RI_sternheimer_fd_solver_test|MODULE_RI_sternheimer_delta_test)$'
```

Expected: all selected tests pass with no warnings or crashes.

- [ ] **Step 2: Build and fingerprint the ABACUS executable on df**

Stage the clean branch, build with the accepted `master_ghj` Intel/MPI profile, and record source commit, executable path, `stat`, and SHA256. Verify the binary contains both input labels before any calculation.

- [ ] **Step 3: Run the N2 p5 lowest-frequency regression once**

Reuse the archived N2 p5 QZTP8 geometry, SG15 PP, FD8 grid, frequency, and product-PCA definition. Use:

```text
rpa_abfs_preorth              onsite_coulomb
rpa_abfs_preorth_threshold    1e-2
```

and the default-on `fd_spectral` preconditioner with regularization zero. Before submission, require script syntax, `sbatch --test-only`, active-job/fingerprint, and destination-directory duplicate checks. Submit exactly one normal `p1` calculation; do not rerun the archived external-preorth reference.

- [ ] **Step 4: Compare runtime and archived reader-v1 artifacts**

Require scheduler `COMPLETED/0:0`, `status=success`, `all_converged=yes`, one response matrix, expected transformed `naux`, solved equations `5*naux`, maximum residual at most `1e-6`, and no NaN/OOM/fatal/unconverged marker. Compare runtime versus archived external-preorth artifacts with the v1 parser:

```text
metadata equality                  required
relative Frobenius(V)              <= 1e-6
relative Frobenius(B)              <= 1e-6
max Sternheimer residual           <= 1e-6
fd_spectral regularization         0
```

Record wall time, binary hash, producer/response basis hashes, transformed dimension, and report identity error.

- [ ] **Step 5: Run repository hygiene and attribution checks**

```bash
git diff --check origin/master_ghj...HEAD
git status --short
git log --format='%h Author:%an<%ae> Committer:%cn<%ce> %s' origin/master_ghj..HEAD
```

Every implementation commit must show author `Codex <codex@openai.com>` and committer `AroundPeking <gonghuanjing@iphy.ac.cn>`. The worktree must contain no generated build or regression artifacts.

- [ ] **Step 6: Rebase safely and push only after all gates pass**

Fetch `origin/master_ghj`, rebase the feature branch if the remote moved, rerun the focused tests after conflict resolution, then push the verified branch head to `refs/heads/master_ghj`. Confirm `git ls-remote origin refs/heads/master_ghj` equals local `HEAD`. If SSH port 22 fails, use the configured GitHub key through `ssh.github.com:443`; do not weaken host verification or overwrite an unexpected remote head.

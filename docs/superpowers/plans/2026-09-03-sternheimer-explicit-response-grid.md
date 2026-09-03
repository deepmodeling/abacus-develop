# Explicit Sternheimer Response Grid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three Delta-Sternheimer-only response-grid dimensions and use the accepted implementation to converge the 13 Angstrom graphene plus water RPA adsorption energy to 10 meV.

**Architecture:** Preserve the PBE `PW_Basis` and the existing scalar response-cutoff path.  Pass an explicit integer dimension triplet into the response-grid factory, validate it before allocation, construct a serial response basis with those exact dimensions, and transfer the PBE potential through a full rectangular FFT restriction.  The physical campaign first screens complete response spaces at representative q/frequency points, then runs full adsorbed/slab RPA pairs for energy acceptance.

**Tech Stack:** C++17, ABACUS input parameter framework, `PW_Basis`, GoogleTest, CMake/CTest, Slurm, LibRPA reader-v1.

---

### Task 1: Add And Validate The Input Contract

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_output.cpp`
- Test: `source/source_io/test_serial/read_input_item_test.cpp`

- [ ] **Step 1: Write failing parser tests**

Add assertions that `sternheimer_response_nx`, `sternheimer_response_ny`, and
`sternheimer_response_nz` default to zero, read positive integers, and reject
negative integers through each item's `check_value` callback.

- [ ] **Step 2: Verify the parser target fails**

Build and run `MODULE_IO_read_item_serial`; expect compilation to fail because
the three `Input_para` members and input labels do not exist.

- [ ] **Step 3: Add the three input fields and items**

Add integer members with zero defaults.  Register each item beside
`sternheimer_response_ecutwfc`, document that it affects only the
Delta-Sternheimer response grid, and reject negative values.

- [ ] **Step 4: Verify the parser target passes**

Build and run `MODULE_IO_read_item_serial`; expect all tests to pass.

### Task 2: Construct Explicit Anisotropic Response Grids

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_response_grid.h`
- Modify: `source/source_lcao/module_ri/sternheimer_response_grid.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_response_grid_test.cpp`

- [ ] **Step 1: Write failing grid-selection tests**

Cover exact explicit dimensions, all-zero compatibility, same-as-PBE identity,
mixed zero/positive rejection, negative values, dimensions larger than PBE,
conflict with positive `sternheimer_response_ecutwfc`, and FD8 dimensions
smaller than nine points.

- [ ] **Step 2: Write failing anisotropic restriction tests**

Use a fine `12x10x14` basis and coarse `8x6x10` basis.  Verify preservation of
a common low-G mode and removal without aliasing of a mode excluded only by
the coarse z dimension.  Add controls proving that sparse z does not remove a
retained in-plane mode and that positive/negative Nyquist coefficients combine
correctly on an even target grid.

- [ ] **Step 3: Run the response-grid target and observe the expected failures**

Build and run `MODULE_RI_sternheimer_response_grid_test`; expect failures or
compile errors caused only by the missing explicit-grid API.

- [ ] **Step 4: Implement the minimal explicit-grid API**

Pass `std::array<int, 3>` and FD order to the factory.  Validate the complete
triplet and scalar-cutoff conflict, create the response basis with exact
dimensions, retain `4*pbe_ecutwfc` as the reciprocal information ceiling, use
a complete rectangular FFT restriction for the explicit route, and record
whether the selected source is `pbe`, `cutoff`, or `explicit`.

- [ ] **Step 5: Run the response-grid target**

Expect every response-grid test to pass.

### Task 3: Wire The Producer And Provenance

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] **Step 1: Write a failing producer-level test where practical**

Exercise the public response-grid selection used by the producer and assert
that explicit dimensions are selected without changing the supplied PBE
basis.  If the smoke fixture cannot construct a full `Parameter` cheaply,
retain the behavior test in the focused response-grid target and test the
producer through the integration smoke run in Task 4.

- [ ] **Step 2: Pass the three dimensions into the grid factory**

Reject explicit dimensions when `sternheimer_delta=false`.  Replace duplicated
source-label expressions with the factory's recorded source and emit requested
dimensions, PBE dimensions, actual response dimensions, and point counts in
both diagnostic-only and normal status files.

- [ ] **Step 3: Re-run focused unit tests**

Build and run the input, response-grid, adapter, Delta, periodic-solver, k/q,
and smoke targets; expect zero failures.

### Task 4: Verify And Commit The Implementation

**Files:**
- Verify all modified source and test files
- Update: `docs/superpowers/specs/2026-09-03-sternheimer-explicit-response-grid-design.md` only if implementation details require clarification

- [ ] **Step 1: Run formatting and whitespace checks**

Run the repository formatter on changed C++ files and `git diff --check`;
expect no diagnostics.

- [ ] **Step 2: Run the seven required focused tests**

Run `MODULE_IO_read_item_serial`, `MODULE_RI_sternheimer_response_grid_test`,
`MODULE_RI_sternheimer_abacus_fd_adapter_test`,
`MODULE_RI_sternheimer_delta_test`,
`MODULE_RI_sternheimer_periodic_solver_test`,
`MODULE_RI_sternheimer_kq_test`, and
`MODULE_RI_sternheimer_abacus_st_smoke_test`; expect all tests to pass.

- [ ] **Step 3: Build `abacus_3p` remotely**

Build from an immutable source archive on the selected HPC host.  Record source
commit, archive hash, executable hash, CMake gates, and test results.

- [ ] **Step 4: Commit and integrate into `master_ghj`**

Create implementation commits with Author `Codex <codex@openai.com>` and
Committer `AroundPeking <gonghuanjing@iphy.ac.cn>`.  Integrate the verified
commit into `master_ghj` without touching unrelated dirty work and verify the
final identities and branch head.

### Task 5: Run The 13 Angstrom Response-Grid Campaign

**Files:**
- Create remotely: immutable inputs, Slurm scripts, manifests, diagnostics, and validation summaries under a new campaign directory
- Update after accepted evidence: `/Users/ghj/Downloads/同步空间/AITP_project/h2o_bn_rpa_adsorption/main.tex`

- [ ] **Step 1: Establish duplicate and provenance guards**

Confirm no active job or completed case already has the same source hash,
geometry, q/frequency definition, and response dimensions.  Verify PBE remains
80 Ry, vacuum is 13 Angstrom, and no H2O producer is submitted.

- [ ] **Step 2: Run the same-grid control**

Use explicit dimensions equal to the 80 Ry PBE grid and compare reader-v1
metadata and matrices against the scalar/PBE path.  Require finite complete
output and the established same-grid numerical tolerance.

- [ ] **Step 3: Screen grids in increasing point-count order**

Use full occupied, virtual, and auxiliary spaces.  Increase in-plane and z
dimensions separately at fixed representative q/frequency points, stop failed
or clearly unconverged branches, and record equations, residuals, wall time,
node-hours, and MaxRSS.  Do not use `sternheimer_delta_max_states` or channel
truncation.

- [ ] **Step 4: Run complete adsorbed/slab pairs**

For the lowest candidate that survives screening and a strictly denser
reference, produce complete canonical q and 16-frequency reader-v1 data and
run LibRPA with the same strict-2D Coulomb settings.  Reuse the validated H2O
energy.

- [ ] **Step 5: Apply the 10 meV gate and document**

Accept only when the complete RPA adsorption energy changes by at most 10 meV
against the denser reference and separate x/y/z refinements satisfy the same
limit.  Keep scheduler, producer, numerical, LibRPA, and physical-result gates
separate, update the TeX evidence, and compile the PDF.

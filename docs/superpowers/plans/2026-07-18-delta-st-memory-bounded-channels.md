# Delta-ST Memory-Bounded Channel Execution Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the peak memory of molecular Delta-Sternheimer response calculations enough to complete the 20 Angstrom N2/N calculation without changing the response equations.

**Architecture:** Keep frequency MPI unchanged. Build the occupied-plus-Delta fixed subspace once per spin and share it read-only across auxiliary-channel workers. Add a function-local maximum worker count to the OpenMP channel scheduler, controlled by `ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS`; an unset variable preserves the current all-thread behavior.

**Tech Stack:** C++17, OpenMP, GoogleTest, ABACUS Delta-ST, Slurm on df_dcu.

---

### Task 1: Bound Auxiliary-Channel Concurrency

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_channel_parallel.h`
- Test: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] Add a test that requests two workers while the OpenMP runtime exposes four and asserts that the observed peak concurrency is exactly two.
- [ ] Run `MODULE_RI_sternheimer_abacus_st_smoke_test` on df_dcu and verify the new test fails before the scheduler API accepts a worker limit.
- [ ] Extend `run_sternheimer_channel_tasks` with a `max_workers` argument; reject negative values, interpret zero as the existing runtime default, and use `num_threads` for positive limits.
- [ ] Re-run the targeted test and verify ordered results, exception propagation, and the two-worker cap.

### Task 2: Share the Delta Fixed Subspace

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_delta.h`
- Modify: `source/source_lcao/module_ri/sternheimer_delta.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_delta_test.cpp`

- [ ] Add a test that constructs one fixed subspace, invokes the shared-subspace solver overload, and compares all response components and residuals with the compatibility API.
- [ ] Run `MODULE_RI_sternheimer_delta_test` on df_dcu and verify the new overload is missing.
- [ ] Expose one fixed-subspace builder and one solver overload taking a const reference to that shared subspace.
- [ ] Keep the old solver signature as a compatibility wrapper that constructs one subspace and delegates.
- [ ] In the ABACUS production path, build the fixed subspace once after `delta_subspace_ready` and capture it read-only in the channel worker.
- [ ] Re-run both targeted unit-test executables.

### Task 3: Runtime Control and Diagnostics

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] Parse `ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS` as a positive integer, defaulting to zero/unbounded.
- [ ] Pass the value to the channel scheduler and include `max_workers=<N>` in the `channels_ready` progress event.
- [ ] Verify invalid values fail explicitly rather than silently changing concurrency.

### Task 4: Remote Build and Numerical Regression

**Files:**
- Modify on server: `/public/home/ghj/app/src/abacus-delta-st-postprocess-20260705`
- Build on server: `/public/home/ghj/app/src/abacus-delta-st-postprocess-20260705/build_delta_st_test`

- [ ] Transfer only the changed source and test files while preserving relative paths.
- [ ] Configure a test-enabled build if the production tree still reports zero registered tests.
- [ ] Run the two targeted test executables, then the complete Sternheimer test set.
- [ ] Rebuild `abacus_3p` and record source/binary timestamps.
- [ ] Run a completed smaller N2 case twice with worker limits 1 and 4; require matching reader-v1 response matrices within floating-point roundoff and matching LibRPA EcRPA.

### Task 5: 20 Angstrom Production Retry

**Files:**
- Create on server: new N2 and N retry directories under `/work1/ghj/sternheimer_abacus_tests/n2_delta_st_box_convergence_20260716`
- Modify: `sternheimer_siab_project/main.tex`

- [ ] Stage clean N2/N response and full-Coulomb producer inputs with `exx_ccp_rmesh_times=1` and `rpa_ccp_rmesh_times=1`.
- [ ] Start with `ABACUS_STERNHEIMER_CHANNEL_MAX_WORKERS=16`, six frequency MPI ranks, one rank per normal node, 30 allocated CPUs and 110610 MiB per node.
- [ ] Confirm progress reaches at least one `equation` event and check live memory before allowing the jobs to continue.
- [ ] If memory remains above the node limit, reduce only the worker cap to 12 and preserve every physics/input parameter.
- [ ] Run matching full-Coulomb producers and LibRPA only after all six reader-v1 response files are complete.
- [ ] Record job IDs, peak RSS, wall time, response convergence, EcRPA and N2 binding energy in the TeX document and rebuild the PDF.

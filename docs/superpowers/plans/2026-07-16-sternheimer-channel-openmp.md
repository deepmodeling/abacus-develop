# Sternheimer Auxiliary-Channel OpenMP Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Use all cores of each frequency MPI node by solving independent Sternheimer auxiliary-basis channels with OpenMP while preserving reader-v1 results and deterministic diagnostics.

**Architecture:** Keep MPI ownership over frequencies and the serial spin, occupied-band, and output order. Add a small indexed channel executor that runs a thread-safe callback with OpenMP dynamic scheduling, stores results and exceptions by channel index, and returns to a serial ordered reduction in the response driver.

**Tech Stack:** C++17, OpenMP, ABACUS finite-difference Delta-Sternheimer, GoogleTest, MPI frequency decomposition, LibRPA reader-v1, Slurm on `df_dcu`.

---

## File Map

- Create `source/source_lcao/module_ri/sternheimer_channel_parallel.h`: indexed OpenMP channel executor with ordered results and exception propagation.
- Modify `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`: test ordered concurrent execution and deterministic exception propagation.
- Modify `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: call the executor for each occupied-band channel batch and reduce diagnostics serially.
- Modify `sternheimer_siab_project/main.tex` in the project workspace after remote verification: record the design, equivalence evidence, CPU utilization, and timing.

### Task 1: Test The Indexed Channel Executor

- [ ] Add tests that request four OpenMP threads, execute at least 32 indexed tasks, record peak concurrent workers with atomics, and require result `i` to remain in slot `i`.
- [ ] Add a test with failures at two indices and require the lower-index exception to be re-thrown after all workers leave the parallel region.
- [ ] Sync the test-only change to the remote feature source and build `MODULE_RI_sternheimer_abacus_st_smoke_test`.
- [ ] Verify RED because `sternheimer_channel_parallel.h` and `run_sternheimer_channel_tasks` do not exist.

### Task 2: Implement The Executor

- [ ] Add `run_sternheimer_channel_tasks<Result>(int, Worker)` with a result vector and an exception vector indexed by channel.
- [ ] Use `#pragma omp parallel for schedule(dynamic)` when `_OPENMP` is defined and a serial loop otherwise.
- [ ] Reject negative task counts and re-throw the first indexed exception after the parallel region.
- [ ] Rebuild the focused test remotely and verify GREEN.

### Task 3: Parallelize The Physical Channel Loop

- [ ] Add a private per-equation result containing `delta_wavefunction`, `SolverResult`, and the explicit equation residual.
- [ ] Move RHS construction, Delta or standard solve, and distinct-column response accumulation into the indexed worker.
- [ ] Keep `all_converged`, `solved_equations`, maximum residuals, and `append_chi0_progress_event` in a serial channel-order pass.
- [ ] Format changed C++ files and run `git diff --check` locally without building ABACUS locally.

### Task 4: Remote Build And Focused Regression

- [ ] Sync the changed source files to `/public/home/ghj/app/src/abacus-delta-st-postprocess-20260705`.
- [ ] Build focused Sternheimer tests and `abacus_3p` on a `normal` node with the current Intel/LibRI configuration.
- [ ] Run the focused tests and retain exact pass/fail output.

### Task 5: Physical OMP Equivalence And Timing

- [ ] Prepare one small fixed-frequency Delta-Sternheimer input and duplicate it into separate `omp1` and `omp30` directories.
- [ ] Run both on `normal`, one MPI rank per node, with `MKL_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.
- [ ] Compare every reader-v1 matrix element, frequency, weight, equation count, convergence flag, and maximum residual.
- [ ] Compare wall time and Slurm CPU utilization; require visible multi-core use before using the build for the N2 convergence campaign.

### Task 6: Documentation And Commit

- [ ] Record implementation and measured evidence in `sternheimer_siab_project/main.tex`, rebuild with XeLaTeX, and inspect the final pages for overflow.
- [ ] Run final focused verification and inspect the complete diff.
- [ ] Commit with `Author: Codex <codex@openai.com>` and `Committer: AroundPeking <gonghuanjing@iphy.ac.cn>`.


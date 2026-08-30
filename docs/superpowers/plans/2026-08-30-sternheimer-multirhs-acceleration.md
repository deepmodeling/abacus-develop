# Sternheimer Multi-RHS Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Batch four independent auxiliary-channel right-hand sides through shared Delta-Sternheimer operator traversals while preserving one scalar GMRES history and convergence decision per column.

**Architecture:** Add column-equivalent batch methods to the FD Hamiltonian, nonlocal projector, fixed-subspace projector, and spectral preconditioner.  Build an independent lockstep batch GMRES and a shared-operator Delta-ST batch entry point, then group stable `owned_channels` into memory-accounted batches in the ABACUS producer.

**Tech Stack:** C++17, GoogleTest, CMake, OpenMP, FFTW, Slurm on df_dcu, ABACUS Delta-Sternheimer, LibRPA.

---

### Task 1: Batch FD Hamiltonian, nonlocal projector, and fixed-subspace projection

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_nonlocal_projector.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.h`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`

- [ ] Add failing tests that call `apply_batch` and `project_batch` and compare every column with the existing scalar calls for empty, one-column, and four-column inputs.
- [ ] Build the focused targets and confirm compilation fails because the batch APIs do not exist.
- [ ] Implement the smallest column-equivalent batch APIs, including explicit size validation and a fused nonlocal-projector block traversal.
- [ ] Run the focused tests and the FD Hamiltonian/RPA neighbors until all pass.
- [ ] Commit the accepted primitive batch APIs.

### Task 2: Independent lockstep batch GMRES

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.h`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`

- [ ] Add a failing test with three diagonal linear systems: one initially converged, one converging before restart, and one requiring a restart.  Compare batch and scalar solutions, convergence flags, iteration counts, and residuals.
- [ ] Build the focused target and confirm compilation fails because `BatchLinearProblem` and `solve_gmres_batch` do not exist.
- [ ] Implement lockstep independent GMRES with active-column compaction at every operator call and unchanged per-column algebra order.
- [ ] Run the focused test, then all RPA and FD solver tests.
- [ ] Commit the accepted batch GMRES.

### Task 3: Shared-operator Delta-ST batch solve

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_fd_preconditioner.h`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_preconditioner.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_delta.h`
- Modify: `source/source_lcao/module_ri/sternheimer_delta.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_delta_test.cpp`

- [ ] Add a failing two-channel Delta-ST test comparing scalar and batch response branches, coefficients, solver records, and physical residuals.
- [ ] Build the focused target and confirm compilation fails because the batch Delta-ST entry point does not exist.
- [ ] Add scalar-equivalent spectral-preconditioner batching and implement one shared projector/operator/denominator path for all columns.
- [ ] Batch the final Hamiltonian residual application while retaining scalar reconstruction formulas.
- [ ] Run Delta, preconditioner, RPA, and periodic solver tests.
- [ ] Commit the accepted batch Delta-ST solve.

### Task 4: Memory-safe production channel grouping

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_channel_resources.h`
- Modify: `source/source_lcao/module_ri/sternheimer_channel_resources.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_runtime_options.h`
- Modify: `source/source_lcao/module_ri/sternheimer_runtime_options.cpp`
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_channel_resources_test.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_runtime_options_test.cpp`
- Modify: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] Add failing tests for default width four, rollback width one, invalid widths, stable channel grouping, and a batch-worker plan whose total estimated storage never exceeds the scalar plan.
- [ ] Build the focused targets and verify the expected missing-API failures.
- [ ] Implement strict diagnostic width parsing, batch grouping, batch-aware worker planning, stable result flattening, and output provenance.
- [ ] Route width-one through the existing scalar implementation and widths above one through the batch solve.
- [ ] Run all directly affected tests and commit production integration.

### Task 5: Local regression and source audit

**Files:**
- Modify only when a test exposes a defect in Tasks 1-4.

- [ ] Configure a clean Debug build with `BUILD_TESTING=ON`, MPI, LibRI, and LibComm.
- [ ] Build the Sternheimer runtime-options, resource, FD Hamiltonian, preconditioner, RPA, Delta, periodic solver, and ABACUS smoke targets.
- [ ] Run every built target directly and require zero failures.
- [ ] Run `git diff --check origin/master_ghj...HEAD` and audit that no solver tolerance, FD order, physical parameter, or output ordering changed.

### Task 6: Remote precision and performance A/B

**Files:**
- Create: `server_jobs/sternheimer_performance_20260830/build_multirhs_dfdcu.slurm`
- Create: `server_jobs/sternheimer_performance_20260830/benchmark_si_multirhs_dfdcu.slurm`
- Create: `server_jobs/sternheimer_performance_20260830/compare_multirhs_response.py`

- [ ] Archive the exact feature commit and stage it under a new immutable `/work1/ghj/sternheimer-performance-20260830` directory without overwriting prior artifacts.
- [ ] Build one RelWithDebInfo artifact on the df_dcu `normal` partition and run the full affected test set.
- [ ] Reuse the accepted Si single-frequency inputs for width-one and width-four A/B jobs with identical MPI/OMP resources and physical parameters.
- [ ] Require complete responses, `all_converged=yes`, residual and matrix gates, then compare process wall time and node-hours excluding queue time.
- [ ] Run the complete accepted Si canonical-q workload only when the single-frequency A/B is both accurate and faster.

### Task 7: Integration and documentation decision

**Files:**
- Create: `server_jobs/sternheimer_performance_20260830/STERNHEIMER_MULTIRHS_BENCHMARK_20260830.md`
- Modify: `/Users/ghj/Downloads/同步空间/AITP_project/delta_st_rpa_project/development_notes/sections/sternheimer_fd_spectral_preconditioner.tex`
- Modify: `/Users/ghj/.codex/skills/abacus-delta-st-efficient-production/SKILL.md`

- [ ] Record commits, executable hashes, test results, response differences, residuals, wall times, resources, node-hours, and speedups.
- [ ] If width four passes accuracy but not performance, keep production default at width one and document the measured result.
- [ ] If width four passes both gates, set it as the no-configuration default, retain width one only for rollback, and merge the attributed commits into `master_ghj`.
- [ ] Rebuild the development-note PDF with XeLaTeX, render the changed pages, and verify the final PDF hash and page count.

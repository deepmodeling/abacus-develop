# Sternheimer SIAB Global-Equation MPI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a deterministic global response-equation MPI layout that mixes all imaginary frequencies across ranks, validate it against the existing grouped layout, and record its molecular speedup separately from solver-tolerance effects.

**Architecture:** Preserve the current grouped implementation as the default. A new input selects global ownership by a modulo mapping over `(occupied state, frequency, auxiliary channel)`; every rank traverses the same loops and solves only its owned equations, while the existing SIAB row gather and canonical writer reconstruct one unchanged target. Solver tolerance and preconditioning remain unchanged during the MPI A/B.

**Tech Stack:** C++14, ABACUS `Input_Item`, MPI, OpenMP, GoogleTest, Slurm on `df_dcu/normal`, XeLaTeX project notes.

---

## File Map

- `source/source_io/module_parameter/input_parameter.h`: store the MPI layout string.
- `source/source_io/module_parameter/read_input_item_output.cpp`: register and validate the input keyword.
- `source/source_io/test_serial/read_input_item_test.cpp`: parser RED/GREEN coverage.
- `source/source_lcao/module_ri/sternheimer_rpa.h`: declare layout parsing and global owner APIs.
- `source/source_lcao/module_ri/sternheimer_rpa.cpp`: implement deterministic ownership.
- `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`: ownership uniqueness and balance tests.
- `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: integrate global ownership, diagnostics, and reductions into the SIAB producer.
- `source/source_lcao/module_ri/test/test_sternheimer_siab_mpi.cpp`: two-rank row assembly regression for globally distributed frequencies.
- `server_jobs/sternheimer_global_equation_mpi_20260727/`: immutable remote build, MPI regression, production A/B, and tolerance-scan Slurm scripts.
- `sternheimer_siab_project/main.tex`: molecular performance ledger and overall acceleration summary.

### Task 1: Add the MPI layout input

- [ ] **Step 1: Write the parser RED test**

Add to `read_input_item_test.cpp` next to `sternheimer_channel_mpi`:

```cpp
{ // sternheimer_mpi_layout
    auto it = find_label("sternheimer_mpi_layout", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());
    EXPECT_EQ(param.input.sternheimer_mpi_layout, "frequency_grouped");
    it->second.str_values = {"global_equation"};
    it->second.read_value(it->second, param);
    EXPECT_EQ(param.input.sternheimer_mpi_layout, "global_equation");
    param.input.sternheimer_mpi_layout = "invalid";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("frequency_grouped or global_equation"));
}
```

- [ ] **Step 2: Run RED remotely**

Submit one `normal` node with `--cpus-per-task=30 --mem=110610M --time=1-00:00:00`. Configure the existing Sternheimer unit-test build and run the serial input test. Expected result: compile failure because `Input_para::sternheimer_mpi_layout` does not exist.

- [ ] **Step 3: Implement the input**

Add to `input_parameter.h`:

```cpp
std::string sternheimer_mpi_layout = "frequency_grouped";
```

Register a string `Input_Item` after `sternheimer_channel_mpi`. Its checker accepts only `frequency_grouped` and `global_equation` and reports both accepted values.

- [ ] **Step 4: Run GREEN remotely**

Re-run the same input test. Expected result: the new parser case passes and existing Sternheimer input cases remain green.

- [ ] **Step 5: Commit**

Commit only the three input files with message:

```text
feat(sternheimer): add MPI layout input
```

### Task 2: Implement deterministic global equation ownership

- [ ] **Step 1: Write the ownership RED tests**

Add tests requiring:

```cpp
EXPECT_EQ(SternheimerRPA::global_equation_owner(0, 0, 0, 16, 428, 32, 0), 0);
EXPECT_EQ(SternheimerRPA::global_equation_owner(0, 0, 31, 16, 428, 32, 0), 31);
EXPECT_EQ(SternheimerRPA::global_equation_owner(0, 0, 32, 16, 428, 32, 0), 0);
EXPECT_EQ(SternheimerRPA::global_equation_owner(0, 0, 0, 16, 428, 32, 1), 1);
```

Loop over `nocc=2`, `nfreq=3`, `nchannel=5`, and `nranks=4`: count every returned owner, assert every equation has one valid owner, and assert `max(count)-min(count) <= 1`. Add invalid-dimension and negative-index cases.

- [ ] **Step 2: Run RED remotely**

Expected result: compile failure because `global_equation_owner` is absent.

- [ ] **Step 3: Implement the owner API**

Use checked 64-bit indexing:

```cpp
const std::int64_t task =
    (static_cast<std::int64_t>(occupied_state) * frequency_count + frequency_index)
        * auxiliary_channel_count
    + auxiliary_channel;
return static_cast<int>((task + normalized_shift) % mpi_ranks);
```

Reject invalid occupied, frequency, channel, rank-count, and dimension values before arithmetic.

- [ ] **Step 4: Run GREEN and grouped regressions remotely**

Run the focused `sternheimer_rpa_test` binary. Expected result: new ownership tests pass and existing `frequency_mpi_assignment`/`channel_group_owner` tests remain unchanged.

- [ ] **Step 5: Commit**

```text
feat(sternheimer): map global response equations to MPI ranks
```

### Task 3: Integrate global ownership into the SIAB producer

- [ ] **Step 1: Add a producer validation RED test**

Extend the existing testable layout validation so that `global_equation` is rejected unless frequency MPI, channel MPI, and SIAB output are all enabled, and accepted when all three are enabled even when `nranks % nfreq != 0`.

- [ ] **Step 2: Run RED remotely**

Expected result: the new validation symbol or overload is missing.

- [ ] **Step 3: Implement layout selection and equation filtering**

Parse `PARAM.inp.sternheimer_mpi_layout` before building grid data. In `global_equation` mode:

```cpp
const int equation_owner = SternheimerRPA::global_equation_owner(
    occupied_state_offset + ib,
    ifrequency,
    ichannel,
    nfreq,
    num_channels,
    GlobalV::NPROC,
    frequency_rank_shift);
if (equation_owner != GlobalV::MY_RANK)
{
    continue;
}
```

All ranks own every frequency for loop traversal and progress timing. The grouped branch continues to use `frequency_mpi_assignment` and `channel_group_owner` without changed semantics.

- [ ] **Step 4: Add load diagnostics**

Accumulate `local_solved_equations` and `local_iteration_sum` before global reductions. Reduce minimum and maximum values across ranks and write:

```text
sternheimer_mpi_layout global_equation
equation_owner_formula occupied_frequency_channel_modulo
rank_local_equations_min <value>
rank_local_equations_max <value>
rank_local_iterations_min <value>
rank_local_iterations_max <value>
```

Retain `solved_equations=6848`, convergence, and residual reductions. Rely on the existing SIAB gather/writer canonical ordering; do not add a second sort.

- [ ] **Step 5: Run unit GREEN remotely**

Run all module_ri Sternheimer tests. Expected result: all tests pass with no changed grouped-layout expectation.

- [ ] **Step 6: Commit**

```text
feat(sternheimer): distribute SIAB equations across all MPI ranks
```

### Task 4: Verify MPI row completeness and numerical identity

- [ ] **Step 1: Extend the two-rank MPI regression**

Construct two frequencies and at least four channels. Assign each row with `global_equation_owner`, gather through `gather_reference_rows_to_root`, and compare the written SIAB file byte-for-byte against a serial row list in a deliberately different order.

- [ ] **Step 2: Run the MPI RED/GREEN cycle remotely**

Run with exactly two MPI ranks. First verify the test fails before producer integration is complete; after integration, require pass. Also run the existing empty-local-row MPI test.

- [ ] **Step 3: Build an immutable ABACUS binary**

On `df_dcu/normal`, build from the committed source in a new directory. Record source commit, compiler/module environment, binary SHA256, and test counts. Do not reuse a path writable by another job.

- [ ] **Step 4: Run small H2 1/2/4-rank gates**

Use one fixed 20 A, 10 Ry H2 input, fixed frequency file, explicit ABFS, identical solver tolerance, and the immutable binary. Compare each global layout to the grouped/serial reference:

- row keys and count exactly equal;
- all equations converged;
- response target and overlap relative Frobenius errors below `1e-10`;
- maximum solver and explicit equation residuals unchanged within floating-point reduction order.

- [ ] **Step 5: Commit regression scripts**

```text
test(sternheimer): validate global equation MPI output
```

### Task 5: Run the full production A/B

- [ ] **Step 1: Freeze the comparison manifest**

Use the completed H2 baseline parameters: 20 A, 50 Ry, 16 fixed frequencies, 428 explicit ABFS channels, 1250 SIAB primitives, full-Coulomb whitening, tolerance `1e-8`, 32 `normal` nodes, one rank/node, 30 OpenMP threads/rank, 110610 MiB/node, and 24 h limit.

- [ ] **Step 2: Run `frequency_grouped` and `global_equation`**

Use the same immutable executable and all identical physical files. Do not submit to `debug`. Save Slurm accounting, progress files, status, `sternheimer_matrix.dat`, and SHA256 manifests.

- [ ] **Step 3: Evaluate correctness and speed**

Require 6848/6848 converged equations, matching provenance, and target/overlap relative Frobenius errors below `1e-10`. Report wall speedup and node-hour ratio. The performance gate is at least `1.5x` lower wall time than 14:06:59 with the same 32 nodes.

- [ ] **Step 4: Commit only reusable runner/comparator changes**

```text
test(sternheimer): benchmark global equation MPI layout
```

### Task 6: Measure solver-tolerance sensitivity separately

- [ ] **Step 1: Run fixed-layout tolerance points**

With `global_equation`, the same immutable binary, and unchanged resources, run `1e-8`, `1e-7`, and `1e-6`. Do not combine this scan with grouped/global speedup numbers.

- [ ] **Step 2: Compare numerical targets**

Use `1e-8` as reference. Record wall time, node-hours, mean/max iterations by frequency, solver/equation residuals, response-target relative Frobenius error, maximum element difference, and SIAB selection loss/rank.

- [ ] **Step 3: Apply the physical gate**

Do not change the default tolerance from `1e-8` based on iteration count alone. A looser value becomes a candidate only after its selected basis is used in the same SOS RPA H/H2 calculation and changes the binding energy by less than `0.05 kcal/mol`, reserving half of the project's `0.1 kcal/mol` convergence budget for other errors.

### Task 7: Update the performance documentation

- [ ] **Step 1: Update the molecular ledger**

Add separate rows for global-equation MPI and each tolerance point, including system, Ecut, grid, frequencies, channels, equations, nodes/ranks/threads, wall time, node-hours, MaxRSS, speedup, and matrix error.

- [ ] **Step 2: Update the overall acceleration summary near the solid table**

Add the molecular global-equation result with an explicit `molecular H2 SIAB` label. Do not report it as a solid speedup. Keep solid OpenMP, symmetry, and k-point results in their existing rows.

- [ ] **Step 3: Compile and inspect the PDF**

Run `latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex`, render the affected pages, and inspect table width, page breaks, labels, absolute times, resources, and numerical gates.

- [ ] **Step 4: Commit documentation changes in their owning repository**

Use the established Codex author and AroundPeking committer identities and verify the final commit attribution.

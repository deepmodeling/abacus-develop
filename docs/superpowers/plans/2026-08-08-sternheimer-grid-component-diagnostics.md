# Sternheimer Grid-Component Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in Delta-ST grid diagnostics that decompose the response and Hamiltonian, then use PCA `1e-4` N2/N 40/50 Ry runs to quantify each component's effect on LibRPA correlation energies.

**Architecture:** Preserve the existing total reader-v1 path byte-for-byte when diagnostics are off. Extend the Delta grid-matrix assembly to retain `T`, `Vloc`, and `Vnl`; add a focused diagnostics module for versioned matrix/tensor files and response reconstruction checks; accumulate SOS, Pulay, and Q-space response branches through the existing MPI reduction path. A standalone Python tool builds 40-to-50 Ry hybrid reader-v1 matrices for normal LibRPA energy evaluation.

**Tech Stack:** C++17, GoogleTest, ABACUS input framework, Intel MPI, reader-v1 binary matrices, Python 3/NumPy, LibRPA `sternheimer_rpa`.

---

### Task 1: Add the opt-in input parameter

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_output.cpp`
- Test: `source/source_io/test_serial/read_input_item_test.cpp`

- [ ] **Step 1: Write the failing parser test**

Add after the `sternheimer_delta` test:

```cpp
{ // sternheimer_grid_diagnostics
    auto it = find_label("sternheimer_grid_diagnostics", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());
    EXPECT_FALSE(param.input.sternheimer_grid_diagnostics);
    it->second.str_values = {"true"};
    it->second.read_value(it->second, param);
    EXPECT_TRUE(param.input.sternheimer_grid_diagnostics);
}
```

- [ ] **Step 2: Run the remote parser test and verify RED**

Run the serial input test target on the configured server.  Expected result:
failure because `Input_para` and the input list do not yet contain
`sternheimer_grid_diagnostics`.

- [ ] **Step 3: Add the parameter and input item**

Add to `Input_para`:

```cpp
bool sternheimer_grid_diagnostics = false;
```

Register a Boolean `Input_Item` beside `sternheimer_delta`, default `False`,
available for LCAO Delta-ST LibRPA output.  The description must state that it
writes component response and grid-matrix diagnostics without changing the
normal total output.

- [ ] **Step 4: Run the remote parser test and verify GREEN**

Expected result: parser target passes and the default remains false.

- [ ] **Step 5: Commit**

```bash
git add source/source_io/module_parameter/input_parameter.h \
        source/source_io/module_parameter/read_input_item_output.cpp \
        source/source_io/test_serial/read_input_item_test.cpp
git commit -m "feat(sternheimer): add grid diagnostics input"
```

### Task 2: Decompose the Delta grid Hamiltonian

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_delta.h`
- Modify: `source/source_lcao/module_ri/sternheimer_delta.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_delta_test.cpp`

- [ ] **Step 1: Extend the existing test and verify RED**

In `AssemblesReferenceGridHamiltonianWithAnalyticGradientsAndNonlocalProjector`,
require:

```cpp
EXPECT_NEAR(matrices.kinetic[0].real(), 1.25, 1.0e-14);
EXPECT_NEAR(matrices.local_potential[0].real(), 0.5, 1.0e-14);
EXPECT_NEAR(matrices.nonlocal[0].real(), 1.5, 1.0e-14);
EXPECT_NEAR((matrices.kinetic[0] + matrices.local_potential[0]
             + matrices.nonlocal[0] - matrices.hamiltonian[0]).real(),
            0.0,
            1.0e-14);
```

Run `MODULE_RI_sternheimer_delta_test` remotely.  Expected result: compile
failure because the component fields do not exist.

- [ ] **Step 2: Extend `SternheimerDeltaGridMatrices`**

Use LAPACK column-major vectors:

```cpp
std::vector<Complex> overlap;
std::vector<Complex> kinetic;
std::vector<Complex> local_potential;
std::vector<Complex> nonlocal;
std::vector<Complex> hamiltonian;
```

- [ ] **Step 3: Accumulate each term independently**

Inside `assemble_delta_sternheimer_grid_matrices`, accumulate `T`, `Vloc`, and
`Vnl` separately, assign each component vector, and form `hamiltonian` by their
sum.  Do not change the analytic-gradient kinetic expression or the nonlocal
projector implementation.

- [ ] **Step 4: Verify GREEN and existing Delta tests**

Run `MODULE_RI_sternheimer_delta_test`.  Expected result: all tests pass,
including the pre-existing total-Hamiltonian value.

- [ ] **Step 5: Commit**

```bash
git add source/source_lcao/module_ri/sternheimer_delta.h \
        source/source_lcao/module_ri/sternheimer_delta.cpp \
        source/source_lcao/module_ri/test/sternheimer_delta_test.cpp
git commit -m "feat(sternheimer): expose delta grid operator components"
```

### Task 3: Add reusable diagnostic data and writers

**Files:**
- Create: `source/source_lcao/module_ri/sternheimer_grid_diagnostics.h`
- Create: `source/source_lcao/module_ri/sternheimer_grid_diagnostics.cpp`
- Create: `source/source_lcao/module_ri/test/sternheimer_grid_diagnostics_test.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Write failing tests for reconstruction and file headers**

Test deterministic complex matrices with

```cpp
const double error = relative_component_reconstruction_error(total, sos, pulay, qspace);
EXPECT_NEAR(error, 0.0, 1.0e-15);
```

Perturb one entry and require a nonzero error.  Write a two-dimensional
operator diagnostic and require the text to contain version `1`, grid
dimensions, spin, matrix labels `overlap`, `kinetic`, `local_potential`,
`nonlocal`, and `hamiltonian`, and indexed complex entries.

Run the new target remotely.  Expected result: target or API is missing.

- [ ] **Step 2: Define focused diagnostic types**

Add types for:

```cpp
struct SternheimerGridDiagnosticMetadata {
    int nx, ny, nz, spin, occupied, virtuals, auxiliaries;
    double volume_element;
};
struct SternheimerPerturbationTensor {
    int occupied, virtuals, auxiliaries;
    std::vector<std::complex<double>> values;
};
```

Provide checked row-major indexing for the perturbation tensor and checked
column-major writing for operator matrices.

- [ ] **Step 3: Implement invariant checks and versioned text writers**

Implement:

```cpp
double relative_component_reconstruction_error(...);
double relative_operator_reconstruction_error(const SternheimerDeltaGridMatrices&);
void write_delta_grid_matrices(...);
void write_delta_perturbation_tensor(...);
```

Reject dimension mismatches, nonpositive volume elements, and operator or
response reconstruction errors above their caller-provided thresholds.

- [ ] **Step 4: Verify GREEN**

Run the new diagnostic test and all existing module-RI Sternheimer unit tests.
Expected result: all pass.

- [ ] **Step 5: Commit**

```bash
git add source/source_lcao/module_ri/sternheimer_grid_diagnostics.* \
        source/source_lcao/module_ri/test/sternheimer_grid_diagnostics_test.cpp \
        source/source_lcao/module_ri/CMakeLists.txt \
        source/source_lcao/module_ri/test/CMakeLists.txt
git commit -m "feat(sternheimer): add grid diagnostic writers"
```

### Task 4: Integrate diagnostics into Delta-ST response output

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Test: `source/source_lcao/module_ri/test/test_sternheimer_siab_mpi.cpp`

- [ ] **Step 1: Add a failing MPI/component regression**

Extend the deterministic MPI fixture to form three local component branches,
reduce them with `sternheimer_chi0::reduce_branch_to_root`, symmetrize on rank
zero, and require their sum to equal the reduced total within `1e-12`.  Run the
MPI test with one and two ranks remotely.  Expected result: component helper or
integration API is missing.

- [ ] **Step 2: Validate the input combination**

At response entry, when `sternheimer_grid_diagnostics=true`, require
`out_sternheimer_librpa=true`, `sternheimer_delta=true`, and the LCAO zero-order
route.  Record `sternheimer_grid_diagnostics yes/no` in the status file.

- [ ] **Step 3: Write stable-basis operator data**

Before Delta grid Gram-Schmidt, assemble `S/T/Vloc/Vnl/H` in the ordered raw KS
virtual basis and compute occupied-virtual grid overlaps.  Rank zero writes
`STERNHEIMER_DELTA_GRID_MATRICES_spin_<spin>.dat` after the `1e-12` operator
invariant passes.

- [ ] **Step 4: Compute the frequency-independent perturbation tensor once**

For each spin, distribute `(occupied, auxiliary)` pairs independently of
frequency, compute all ordered raw-KS-virtual overlaps
`<eta_a|P_mu,h|psi_i>`, gather them to rank zero, check for missing or duplicate
rows, and write `STERNHEIMER_DELTA_PERTURBATION_spin_<spin>.dat`.

- [ ] **Step 5: Accumulate component response branches**

Allocate three extra `nfreq * naux * naux` branch arrays only when diagnostics
are enabled.  Within each owned Delta equation, accumulate
`in_sos_wavefunction`, `in_pulay_wavefunction`, and `out_wavefunction` through
the same `accumulate_chi0_branch_column` call used by the total response.

- [ ] **Step 6: Reduce, validate, and write component reader-v1 files**

Reduce each component with the existing frequency communicator.  On each
frequency leader, symmetrize all four matrices, require relative total-minus-
sum error below `1e-10`, and write `sos`, `pulay`, and `qspace` files with the
normal metadata and auxiliary ordering.  Keep the existing total filename and
index unchanged.

- [ ] **Step 7: Verify diagnostic-off identity and MPI equality**

Run the existing H2 fixture once with the frozen parent executable and once
with the new executable with diagnostics off; strip timing-only status lines
and require byte identity for all reader-v1 matrices.  Run diagnostics on with
one rank and global-equation MPI and require relative Frobenius differences
below `1e-12` for every component.

- [ ] **Step 8: Commit**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp \
        source/source_lcao/module_ri/test/test_sternheimer_siab_mpi.cpp
git commit -m "feat(sternheimer): output response grid components"
```

### Task 5: Build hybrid reader-v1 matrices

**Files:**
- Create: `tools/sternheimer_component_hybrid.py`
- Create: `tools/test_sternheimer_component_hybrid.py`

- [ ] **Step 1: Write failing Python tests**

Generate small synthetic reader-v1 matrices and require the tool to build

```text
hybrid_sos = total40 - sos40 + sos50
hybrid_pulay = total40 - pulay40 + pulay50
hybrid_qspace = total40 - qspace40 + qspace50
```

while preserving marker, metadata, channels, frequencies, and weights.  Require
rejection of mismatched basis order, dimensions, or frequency metadata.

- [ ] **Step 2: Verify RED**

Run `python3 tools/test_sternheimer_component_hybrid.py`.  Expected result:
import or executable missing.

- [ ] **Step 3: Implement checked reader-v1 parsing and writing**

Use `struct` and NumPy, preserve all non-matrix bytes from the 40 Ry total file,
replace only the complex matrix payload, and emit a CSV manifest with input
hashes and relative component differences.

- [ ] **Step 4: Verify GREEN**

Run the Python tests and compare a no-op replacement against the original file
byte-for-byte.

- [ ] **Step 5: Commit**

```bash
git add tools/sternheimer_component_hybrid.py tools/test_sternheimer_component_hybrid.py
git commit -m "feat(sternheimer): build component hybrid matrices"
```

### Task 6: Remote validation and PCA 1e-4 production diagnosis

**Files:**
- Create locally in the TeX project: `scripts/run_n2_n_grid_components_pca1em4_df_20260808.slurm`
- Update locally in the TeX project: `sections/n2_grid_convergence_20260808.tex`

- [ ] **Step 1: Build and run all relevant tests remotely**

Build the feature branch on `df_iopcas_ghj` using the established ABACUS Intel
MPI profile.  Run parser, Delta, RPA, diagnostics, and MPI tests.  Record source
commit, source archive hash, executable hash, commands, and complete results.

- [ ] **Step 2: Run the H2 diagnostic smoke**

Use one low-cost frequency first.  Require normal total output, three component
files, both spin-independent diagnostic text files, all invariant checks, and
successful LibRPA total energy.

- [ ] **Step 3: Submit N2/N 40 and 50 Ry jobs**

Use PCA `1e-4`, fixed 12-point GreenX frequencies, TZDP-8au, full-Ewald LibRPA,
`ks_bands`, `global_equation`, and the established node/OMP policy.  Set
`sternheimer_grid_diagnostics 1`.  N2 and N remain independent convergence
targets.

- [ ] **Step 4: Build hybrid matrices and run LibRPA**

For N2 and N separately, run total40, total50, hybrid-SOS, hybrid-Pulay, and
hybrid-Q LibRPA calculations.  Report component matrix changes, component
energy influences, the full 40-to-50 change, and the nonlinear remainder.

- [ ] **Step 5: Update and compile the research note**

Add tables for `S/T/Vloc/Vnl/B/M` changes and `Ec` influences, state the Ecut
decision, compile with XeLaTeX, and visually inspect the affected pages.

- [ ] **Step 6: Final commit**

Commit code and durable documentation with Codex as author and AroundPeking as
committer, then verify attribution and the clean worktree.

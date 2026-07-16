# Delta-Sternheimer Collinear `nspin=2` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the Gamma-point LCAO Delta-Sternheimer path so independent spin-up and spin-down responses are solved with their own Hamiltonians and projectors and summed into one LibRPA v1 response.

**Architecture:** Keep grid, AO samples, auxiliary perturbations, and frequency metadata shared. Process occupied spin channels sequentially so each channel owns its Hamiltonian, occupied projector, and Delta subspace; add each result to rank-local frequency accumulators and write each frequency only after all spins have contributed.

**Tech Stack:** C++17, ABACUS LCAO/finite-difference Sternheimer modules, GoogleTest, MPI frequency decomposition, LibRPA reader-v1 output, Slurm remote builds and smoke tests.

---

## File Map

- Modify `source/source_lcao/module_ri/sternheimer_rpa.h`: expose transition-window union.
- Modify `source/source_lcao/module_ri/sternheimer_rpa.cpp`: implement validated transition-window union.
- Modify `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`: verify spin-window union and additive occupation-weighted response.
- Modify `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h`: expose pure spin-channel diagnostics helpers.
- Modify `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`: require two occupied spin channels and per-spin band diagnostics.
- Modify `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: replace the one-spin restriction with sequential spin contexts and deferred v1 writing.
- Modify `sternheimer_siab_project/main.tex` in the project workspace after remote verification: record equations, implementation mapping, tests, and N-atom evidence.

### Task 1: Define And Test Spin-Channel Metadata

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h`
- Test: `source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp`

- [ ] **Step 1: Write the failing two-spin metadata test**

Extend `ValidatesSpinResolvedLCAOOccupiedChannels` with a spin-down channel containing one occupied band, while spin-up contains four. Require both channels to survive validation and require helpers that return `{0,1}` and `{4,1}`:

```cpp
ModuleRI::SternheimerLCAOOccupiedChannel spin_down;
spin_down.spin_index = 1;
spin_down.coefficients = {{Complex(0.0, 0.0), Complex(1.0, 0.0), Complex(0.0, 0.0)}};

spin_up.coefficients.resize(4, spin_up.coefficients.front());
const std::vector<ModuleRI::SternheimerLCAOOccupiedChannel> quartet = {spin_up, spin_down};

EXPECT_NO_THROW(ModuleRI::validate_sternheimer_lcao_occupied_channels(quartet, 2, 3));
EXPECT_EQ(ModuleRI::sternheimer_lcao_spin_indices(quartet), (std::vector<int>{0, 1}));
EXPECT_EQ(ModuleRI::sternheimer_lcao_occupied_bands_per_spin(quartet), (std::vector<int>{4, 1}));
EXPECT_EQ(ModuleRI::sternheimer_lcao_total_occupied_bands(quartet), 5);
```

- [ ] **Step 2: Run the focused test remotely and verify RED**

Sync only the test change to the feature-branch server source and run:

```bash
ctest --test-dir build_delta_st_test -R sternheimer_abacus_st_smoke --output-on-failure
```

Expected: compile failure because `sternheimer_lcao_spin_indices` and `sternheimer_lcao_occupied_bands_per_spin` do not exist.

- [ ] **Step 3: Implement the metadata helpers**

Add these inline helpers after validation:

```cpp
inline std::vector<int> sternheimer_lcao_spin_indices(
    const std::vector<SternheimerLCAOOccupiedChannel>& channels)
{
    std::vector<int> result;
    result.reserve(channels.size());
    for (const auto& channel: channels)
    {
        result.push_back(channel.spin_index);
    }
    return result;
}

inline std::vector<int> sternheimer_lcao_occupied_bands_per_spin(
    const std::vector<SternheimerLCAOOccupiedChannel>& channels)
{
    std::vector<int> result;
    result.reserve(channels.size());
    for (const auto& channel: channels)
    {
        result.push_back(static_cast<int>(channel.coefficients.size()));
    }
    return result;
}
```

- [ ] **Step 4: Rebuild and verify GREEN remotely**

Run the same focused CTest. Expected: one matching test target passes with zero failures.

- [ ] **Step 5: Commit the metadata slice**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h \
        source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp
git commit -m "test(sternheimer): define spin response metadata"
```

### Task 2: Union Spin-Resolved Minimax Windows

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.h`
- Modify: `source/source_lcao/module_ri/sternheimer_rpa.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`

- [ ] **Step 1: Write the failing union test**

```cpp
TEST(SternheimerRPA, MergeTransitionEnergyWindowsCoversAllSpinChannels)
{
    using Window = ModuleRI::SternheimerRPA::TransitionEnergyWindow;
    const auto merged = ModuleRI::SternheimerRPA::merge_transition_energy_windows(
        {Window{0.4, 2.0}, Window{0.2, 3.5}});
    EXPECT_DOUBLE_EQ(merged.emin_ha, 0.2);
    EXPECT_DOUBLE_EQ(merged.emax_ha, 3.5);
    EXPECT_THROW(ModuleRI::SternheimerRPA::merge_transition_energy_windows({}),
                 std::invalid_argument);
}
```

- [ ] **Step 2: Run the focused test remotely and verify RED**

```bash
ctest --test-dir build_delta_st_test -R sternheimer_rpa --output-on-failure
```

Expected: compile failure because `merge_transition_energy_windows` is absent.

- [ ] **Step 3: Implement validated union**

Declare in the public class and implement:

```cpp
SternheimerRPA::TransitionEnergyWindow SternheimerRPA::merge_transition_energy_windows(
    const std::vector<TransitionEnergyWindow>& windows)
{
    if (windows.empty())
    {
        throw std::invalid_argument("Sternheimer minimax window merge requires at least one spin channel.");
    }
    TransitionEnergyWindow result = windows.front();
    for (const TransitionEnergyWindow& window: windows)
    {
        if (window.emin_ha <= 0.0 || window.emax_ha < window.emin_ha)
        {
            throw std::invalid_argument("Sternheimer minimax spin window is invalid.");
        }
        result.emin_ha = std::min(result.emin_ha, window.emin_ha);
        result.emax_ha = std::max(result.emax_ha, window.emax_ha);
    }
    return result;
}
```

- [ ] **Step 4: Verify GREEN remotely**

Run the focused CTest and the complete `sternheimer_rpa` target. Expected: zero failures.

- [ ] **Step 5: Commit the minimax slice**

```bash
git add source/source_lcao/module_ri/sternheimer_rpa.h \
        source/source_lcao/module_ri/sternheimer_rpa.cpp \
        source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp
git commit -m "feat(sternheimer): merge spin transition windows"
```

### Task 3: Refactor The Response Driver To Sum Independent Spins

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`
- Test: `source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp`

- [ ] **Step 1: Add an additive spin-response regression test**

Call `accumulate_chi0_branch_column` once with occupation one for a synthetic
spin-up response and once with occupation one for a distinct spin-down
response. Require the final branch matrix to equal the elementwise sum of two
separately accumulated matrices, and require one final call to
`symmetrize_chi0_imaginary_frequency` to equal the sum of separately
symmetrized spin matrices.

- [ ] **Step 2: Verify the algebra test passes before orchestration changes**

Run `sternheimer_rpa`. Expected: pass, documenting that the existing matrix
primitive is additive and must not receive a new spin factor.

- [ ] **Step 3: Remove the one-occupied-spin restriction**

Delete the `lcao_occupied_channels->size() != 1` failure. Build a vector of
response channels: all validated LCAO channels for LCAO Delta-ST, or one
synthetic spin-zero descriptor for the existing FD path.

- [ ] **Step 4: Build a global transition window**

For every response spin, call
`transition_energy_window_from_eigenvalues_ry` with that spin's full
eigenvalue and occupation rows. Merge the resulting windows before GreenX grid
generation. External frequency files remain unchanged.

- [ ] **Step 5: Allocate owner-rank frequency accumulators**

Create `std::vector<std::vector<Complex>> chi0_branches(nfreq)`. Allocate
`naux*naux` zeros only for frequencies whose owner is `MY_RANK`. Record each
owned frequency start time before entering the spin loop.

- [ ] **Step 6: Move spin-specific state into a sequential loop**

For each response spin channel:

```cpp
const int spin_index = response_channel.spin_index;
const int occupied_count = occupied_band_count(elec_state, spin_index);
const SternheimerFDHamiltonian hamiltonian = use_frequency_mpi
    ? make_sternheimer_fd_full_hamiltonian(potential, pw_basis, ucell, spin_index, 1.0)
    : make_sternheimer_fd_hamiltonian(potential, pw_basis, ucell, spin_index, 1.0);
```

Construct that spin's LCAO occupied functions from the shared sampled AO
functions, normalize values and gradients together, orthonormalize only those
occupied functions, and build the reference Delta subspace with the same
spin-specific Hamiltonian and occupied projector.

- [ ] **Step 7: Accumulate every owned frequency without writing**

For each owned frequency, occupied band, and auxiliary perturbation, solve with
the current spin's Hamiltonian, projector, eigenvalue, and occupation. Pass the
same Hartree potentials and current `wg(spin_index, ib)` to
`accumulate_chi0_branch_column`, targeting `chi0_branches[ifrequency]`.

Do not multiply by `2/nspin`; do not symmetrize inside the spin loop.

- [ ] **Step 8: Release spin data and process the next spin**

Keep the Hamiltonian, states, projector functions, and Delta subspace scoped to
the loop body. Preserve only integer diagnostics per spin: spin index,
occupied-band count, projector dimension, virtual-state count, accepted
candidates, and discarded candidates.

- [ ] **Step 9: Write each owned frequency exactly once**

After the spin loop, symmetrize every owned `chi0_branches[ifrequency]`, write
one v1 file, and emit one `frequency_finish` event. Keep the existing
frequency-owner filename and index mapping unchanged.

- [ ] **Step 10: Replace scalar status fields with spin-resolved fields**

Write:

```text
sternheimer_response_spin_channels 2
sternheimer_response_spin_indices 1 2
occupied_bands_per_spin 4 1
occupied_bands 5
occupied_projector_dimensions_per_spin 4 1
sternheimer_delta_virtual_states_per_spin ...
sternheimer_delta_accepted_candidates_per_spin ...
sternheimer_delta_discarded_candidates_per_spin ...
```

For compatibility, retain `sternheimer_response_spin_channel <index>` only
when exactly one response channel is present.

- [ ] **Step 11: Run formatting and static checks locally**

```bash
clang-format -i source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp \
                source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h \
                source/source_lcao/module_ri/sternheimer_rpa.cpp \
                source/source_lcao/module_ri/sternheimer_rpa.h \
                source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp \
                source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp
git diff --check
```

- [ ] **Step 12: Commit the response-driver slice after remote verification**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp \
        source/source_lcao/module_ri/sternheimer_abacus_st_smoke.h \
        source/source_lcao/module_ri/sternheimer_rpa.cpp \
        source/source_lcao/module_ri/sternheimer_rpa.h \
        source/source_lcao/module_ri/test/sternheimer_abacus_st_smoke_test.cpp \
        source/source_lcao/module_ri/test/sternheimer_rpa_test.cpp
git commit -m "feat(sternheimer): sum collinear spin responses"
```

### Task 4: Remote Build And Focused Regression

**Files:**
- Verify: remote feature-branch source and `build_delta_st_test/abacus_3p`

- [ ] **Step 1: Verify source and executable provenance**

Record local and server branch, SHA, source path, executable path, CMake LibRI,
LibComm, and debug-info flags, and binary timestamp. Use the feature-branch
exception rather than claiming a `master_ghj` match.

- [ ] **Step 2: Build on the remote server**

Use the existing validated Delta-ST CMake cache and build `abacus_3p`; do not
compile ABACUS locally.

- [ ] **Step 3: Run focused tests**

```bash
ctest --test-dir build_delta_st_test \
  -R 'sternheimer_(rpa|abacus_st_smoke|delta|fd_hamiltonian)' \
  --output-on-failure
```

Expected: four matching test targets, zero failures.

- [ ] **Step 4: Re-run the spin-polarized H smoke**

Use the existing small H `nspin=2`, `nupdown=1` smoke input. Require success,
one response spin index, unchanged equation count, converged equations, and a
v1 matrix equal to the pre-change file within `1e-12` relative Frobenius norm.

### Task 5: Quartet N Atom End-To-End Validation

**Files:**
- Create remotely: a clean N-atom smoke directory derived from the validated H smoke input
- Use: SG15 N ONCV PBE pseudopotential and N TZDP-8au orbital

- [ ] **Step 1: Prepare the N quartet input**

Use one N atom in a modest smoke-test box, Gamma-only LCAO, `nspin=2`,
`nupdown=3`, fixed integer occupations, `out_sternheimer_librpa=1`,
`sternheimer_delta=true`, one frequency, and loose development-only auxiliary
settings. Do not use this smoke result as a converged N reference.

- [ ] **Step 2: Submit only to the normal partition**

Request one full node, explicit maximum normal-node memory and wall time, one
MPI rank per node, and all node cores as OpenMP threads. Check `squeue` after
submission and correct `BadConstraints` before waiting.

- [ ] **Step 3: Verify spin physics and output**

Require SCF convergence, total magnetization three, occupied-band diagnostics
`4 1`, two response spin indices, `N_ABFS*(4+1)` solved equations, all solvers
converged, and one valid v1 file.

- [ ] **Step 4: Verify frequency MPI equivalence**

Run a small six-frequency N case using six ranks/nodes and a one-rank serial or
single-owner comparison. Compare all v1 matrix elements and require relative
Frobenius difference below `1e-12`.

- [ ] **Step 5: Run LibRPA on the complete N output**

Use the same full-Coulomb v1 producer convention as the H2/H workflow. Confirm
LibRPA reads every frequency and reports a finite RPA correlation energy.

### Task 6: Documentation And Final Commit

**Files:**
- Modify: `/Users/ghj/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.tex`

- [ ] **Step 1: Record the implementation in the TeX project**

Add the spin-resolved equations, the absence of cross-spin terms, ABACUS
occupation convention, per-spin Hamiltonian/projector/Delta construction,
memory-preserving loop order, diagnostics, commits, remote build evidence, H
regression, and N quartet smoke results.

- [ ] **Step 2: Render and inspect the TeX document locally**

Use the project's existing render script. Check the new equations and tables
for overflow and verify the table of contents still includes the relevant
section.

- [ ] **Step 3: Verify git state and attribution**

Run focused tests again, `git diff --check`, and inspect every commit with:

```bash
git log -3 --format='%h Author:%an\ <%ae\> Committer:%cn\ <%ce\> %s'
```

Expected author is `Codex <codex@openai.com>` and committer is
`AroundPeking <gonghuanjing@iphy.ac.cn>`.

- [ ] **Step 4: Start the N2 campaign only after the N atom gate passes**

The next calculation stage is N2/N binding-energy convergence with the Thesis
bond length and frozen-core reference. Do not submit N2 production runs while
the quartet N atom response or MPI equivalence remains unverified.

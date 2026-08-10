# Efficient Eighth-Order Sternheimer Finite-Difference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an efficient eighth-order Cartesian finite-difference kinetic operator to the molecular Delta-Sternheimer route and determine, with cost-ordered gates, whether it should replace sixth order for production N/N2 calculations.

**Architecture:** Keep the existing Hamiltonian interface and select a radius-specialized local kernel once per `apply()` call. Precompute only one-dimensional shifted coordinate tables in the constructor, so the inner grid loop performs fixed-coefficient neighbor reads without order branches, full-grid work arrays, or additional MPI communication. Validate the operator independently before using the existing fixed-frequency, fixed-basis, fixed-full-Ewald N/N2 workflow.

**Tech Stack:** ABACUS-compatible C++11/14, GoogleTest, OpenMP, ABACUS input system, Slurm on `df_iopcas_ghj`, LibRPA Sternheimer task, Python/TeX result postprocessing.

---

## File Map

- Modify `source/source_io/module_parameter/input_parameter.h`: document the accepted finite-difference orders.
- Modify `source/source_io/module_parameter/read_input_item_output.cpp`: accept order 8 and reject all other unsupported values.
- Modify `source/source_io/test_serial/read_input_item_test.cpp`: cover accepted order 8 and rejected order 10.
- Modify `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h`: store compact shifted-coordinate tables and declare the radius-specialized kernel.
- Modify `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp`: build the tables, dispatch once by order, and apply the second-/fourth-/sixth-/eighth-order stencils.
- Modify `source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp`: test eighth-order accuracy, Hermiticity, OpenMP consistency, and unsupported values.
- Create `.codex_tmp/build_fd8_${short_commit}.slurm`: server-only build and regression gate; keep untracked.
- Create `.codex_tmp/run_fd8_h2_smoke_${short_commit}.slurm`: H2 order-2/6/8 numerical and timing gate; keep untracked.
- Create `.codex_tmp/run_n_fd8_grid_diag_40_50_${short_commit}.slurm`: N gate before any N2 submission; keep untracked.
- Create `.codex_tmp/run_n2_fd8_grid_diag_40_50_${short_commit}.slurm`: conditional N2 production gate; keep untracked.
- Modify `/Users/ghj/Downloads/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.tex`: record provenance, component convergence, timings, and the production-order decision after all gates finish.

### Task 1: Extend the Input Contract

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h:585`
- Modify: `source/source_io/module_parameter/read_input_item_output.cpp:1031-1048`
- Test: `source/source_io/test_serial/read_input_item_test.cpp:1555-1575`

- [ ] **Step 1: Change the serial input test to require order 8 and reject order 10**

Replace the current order-8 rejection block with:

```cpp
it->second.str_values = {"8"};
it->second.read_value(it->second, param);
EXPECT_EQ(param.input.sternheimer_fd_order, 8);

param.input.sternheimer_fd_order = 10;
testing::internal::CaptureStdout();
EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
output = testing::internal::GetCapturedStdout();
EXPECT_THAT(output, testing::HasSubstr("must be 2, 4, 6, or 8"));
```

- [ ] **Step 2: Stage a server-side red test**

Generate a patch from the worktree, apply it to a clean copy of the current source on `df_iopcas_ghj`, and build only `MODULE_IO_read_item_serial`:

```bash
git diff --binary > .codex_tmp/fd8_input_red.patch
cmake --build "$build" --parallel 40 --target MODULE_IO_read_item_serial
ctest --test-dir "$build" --output-on-failure -R '^MODULE_IO_read_item_serial$'
```

Expected: `MODULE_IO_read_item_serial` fails because order 8 is still rejected.

- [ ] **Step 3: Extend the accepted value set without changing the default**

Use this validation condition and message:

```cpp
if (para.input.sternheimer_fd_order != 2 && para.input.sternheimer_fd_order != 4
    && para.input.sternheimer_fd_order != 6 && para.input.sternheimer_fd_order != 8)
{
    ModuleBase::WARNING_QUIT("ReadInput", item.label + " must be 2, 4, 6, or 8.");
}
```

Update the description and comment to say `2, 4, 6, or 8`; keep `default_value = "2"` and `sternheimer_fd_order = 2` unchanged.

- [ ] **Step 4: Run the focused input test remotely**

```bash
cmake --build "$build" --parallel 40 --target MODULE_IO_read_item_serial
ctest --test-dir "$build" --output-on-failure -R '^MODULE_IO_read_item_serial$'
```

Expected: one selected test passes.

### Task 2: Add Failing Eighth-Order Hamiltonian Tests

**Files:**
- Test: `source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp:305-353`

- [ ] **Step 1: Add the eighth-order plane-wave test**

Add a test using the same 32-point, mode-3 periodic wave as the sixth-order test:

```cpp
TEST(SternheimerFDHamiltonian, EighthOrderPeriodicLaplacianFurtherReducesPlaneWaveError)
{
    constexpr int nx = 32;
    constexpr int mode = 3;
    const double length = 2.0 * std::acos(-1.0);
    const double spacing = length / nx;
    Hamiltonian::Grid grid{nx, 1, 1, spacing, 1.0, 1.0, true};
    const std::vector<double> potential(grid.size(), 0.0);
    Hamiltonian sixth_order(grid, potential, 1.0, nullptr, 6);
    Hamiltonian eighth_order(grid, potential, 1.0, nullptr, 8);

    Hamiltonian::Vector plane_wave(grid.size());
    for (int ix = 0; ix != nx; ++ix)
    {
        const double phase = static_cast<double>(mode) * spacing * ix;
        plane_wave[static_cast<std::size_t>(ix)] = Complex(std::cos(phase), std::sin(phase));
    }

    Hamiltonian::Vector sixth_action;
    Hamiltonian::Vector eighth_action;
    sixth_order.apply(plane_wave, sixth_action);
    eighth_order.apply(plane_wave, eighth_action);

    Complex sixth_rayleigh(0.0, 0.0);
    Complex eighth_rayleigh(0.0, 0.0);
    for (int ix = 0; ix != nx; ++ix)
    {
        sixth_rayleigh += std::conj(plane_wave[ix]) * sixth_action[ix];
        eighth_rayleigh += std::conj(plane_wave[ix]) * eighth_action[ix];
    }
    sixth_rayleigh /= static_cast<double>(nx);
    eighth_rayleigh /= static_cast<double>(nx);
    const double exact = static_cast<double>(mode * mode);

    EXPECT_EQ(eighth_order.finite_difference_order(), 8);
    EXPECT_LT(std::abs(eighth_rayleigh.real() - exact),
              0.1 * std::abs(sixth_rayleigh.real() - exact));
    EXPECT_NEAR(eighth_rayleigh.imag(), 0.0, 1.0e-12);
}
```

- [ ] **Step 2: Parameterize Hermiticity for periodic and nonperiodic order-8 grids**

Add a helper that constructs `Grid{6, 5, 4, 0.4, 0.5, 0.6, periodic}`, calls `dense_matrix(1000)`, and checks every `H[row][col] - conj(H[col][row])` component below `1e-12`. Call it from two tests with `periodic=true` and `periodic=false`.

- [ ] **Step 3: Extend the existing OpenMP equality test to order 8**

Construct the Hamiltonian in `ApplyUsesRequestedOpenMPThreadsAndMatchesSerial` with final argument `8`. Preserve its existing serial-versus-threaded `1e-12` element-wise check and thread-count assertion.

- [ ] **Step 4: Make the unsupported-order test reject 10**

```cpp
EXPECT_THROW(Hamiltonian(grid, std::vector<double>(grid.size(), 0.0), 1.0, nullptr, 10),
             std::invalid_argument);
```

- [ ] **Step 5: Run the Hamiltonian test remotely and confirm the red state**

```bash
cmake --build "$build" --parallel 40 --target MODULE_RI_sternheimer_fd_hamiltonian_test
ctest --test-dir "$build" --output-on-failure -R '^MODULE_RI_sternheimer_fd_hamiltonian_test$'
```

Expected: construction with order 8 fails before the new accuracy and Hermiticity assertions can pass.

### Task 3: Implement the Compact Shift Tables and Specialized Kernel

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h:5-74`
- Modify: `source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp:23-274`

- [ ] **Step 1: Add compact shift-table storage**

Add `<array>` and these private declarations:

```cpp
static constexpr int max_stencil_radius_ = 4;
using AxisShiftTable = std::array<std::vector<int>, max_stencil_radius_>;

void initialize_shift_tables();
template <int Radius>
void apply_local(const Vector& psi, Vector& hpsi, int* threads_used) const;

AxisShiftTable x_plus_;
AxisShiftTable x_minus_;
AxisShiftTable y_plus_;
AxisShiftTable y_minus_;
AxisShiftTable z_plus_;
AxisShiftTable z_minus_;
```

Remove the private `shifted_index()` declaration after all calls have moved to the tables. The six table groups contain exactly `2 * radius * (nx + ny + nz)` integers for the selected order.

- [ ] **Step 2: Validate order 8 and initialize tables in the constructor**

Use the accepted set `2/4/6/8`, update the exception text, and call `initialize_shift_tables()` after scalar/grid validation. For every offset `1..finite_difference_order_/2`, fill each one-dimensional table with a wrapped coordinate for periodic grids and `-1` outside nonperiodic grids:

```cpp
const auto shift = [this](const int coordinate, const int delta, const int extent) {
    const int shifted = coordinate + delta;
    if (grid_.periodic)
    {
        return (shifted % extent + extent) % extent;
    }
    return shifted < 0 || shifted >= extent ? -1 : shifted;
};
```

- [ ] **Step 3: Define fixed stencil coefficients**

Use a radius-specialized helper in the implementation file that returns center plus positive-offset coefficients:

```cpp
template <int Radius>
struct FDStencil;

template <>
struct FDStencil<1> { static std::array<double, 2> coefficients() { return {{-2.0, 1.0}}; } };
template <>
struct FDStencil<2> { static std::array<double, 3> coefficients() { return {{-5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0}}; } };
template <>
struct FDStencil<3> { static std::array<double, 4> coefficients() { return {{-49.0 / 18.0, 3.0 / 2.0, -3.0 / 20.0, 1.0 / 90.0}}; } };
template <>
struct FDStencil<4> { static std::array<double, 5> coefficients() { return {{-205.0 / 72.0, 8.0 / 5.0, -1.0 / 5.0, 8.0 / 315.0, -1.0 / 560.0}}; } };
```

- [ ] **Step 4: Replace the branched loop with `apply_local<Radius>`**

Inside the existing fused x/y/z loop, initialize the center term with `coefficients[0]`, then loop `offset=1..Radius`, reading precomputed coordinates and adding x, y, and z neighbor pairs. Preserve the current addition order: center; x+/x-; y+/y-; z+/z- for each successive offset. Skip a contribution only when its coordinate is `-1`.

Use a single switch before the grid loop:

```cpp
switch (finite_difference_order_)
{
    case 2: apply_local<1>(psi, hpsi, threads_used); break;
    case 4: apply_local<2>(psi, hpsi, threads_used); break;
    case 6: apply_local<3>(psi, hpsi, threads_used); break;
    case 8: apply_local<4>(psi, hpsi, threads_used); break;
    default: throw std::logic_error("Unsupported Sternheimer finite-difference order.");
}
```

Keep `nonlocal_projector_->add_to(psi, hpsi)` in `apply()` after this switch, so the local kernel remains the only optimized region.

- [ ] **Step 5: Run formatting and static diff checks locally without compiling**

```bash
clang-format -i source/source_io/module_parameter/input_parameter.h \
  source/source_io/module_parameter/read_input_item_output.cpp \
  source/source_io/test_serial/read_input_item_test.cpp \
  source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h \
  source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp \
  source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp
git diff --check
```

Expected: no formatting or whitespace errors. Do not compile ABACUS locally.

### Task 4: Run the Server-Only Unit and Regression Gate

**Files:**
- Create: `.codex_tmp/build_fd8_${short_commit}.slurm` (untracked)

- [ ] **Step 1: Commit the implementation with required attribution**

```bash
git add source/source_io/module_parameter/input_parameter.h \
  source/source_io/module_parameter/read_input_item_output.cpp \
  source/source_io/test_serial/read_input_item_test.cpp \
  source/source_lcao/module_ri/sternheimer_fd_hamiltonian.h \
  source/source_lcao/module_ri/sternheimer_fd_hamiltonian.cpp \
  source/source_lcao/module_ri/test/sternheimer_fd_hamiltonian_test.cpp
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m "feat(sternheimer): add efficient eighth-order grid Laplacian"
```

- [ ] **Step 2: Transfer an exact clean source bundle and record provenance**

Create a bundle containing the implementation commit, reconstruct it under a commit-named directory on `df_iopcas_ghj`, and record values generated by:

```bash
commit=$(git rev-parse HEAD)
short_commit=$(git rev-parse --short=10 HEAD)
source_dir=/data/home/df_iopcas_ghj/app/abacus/sternheimer-fd8-${short_commit}-src
printf 'source_commit=%s\nsource_status_clean=yes\nsource_dir=%s\n' \
  "$commit" "$source_dir" > SOURCE_PROVENANCE.txt
```

Reject the build if any line differs from the checked-out source.

- [ ] **Step 3: Build on one full normal node**

Use `#SBATCH --partition=p1`, one task, 40 cores, `--mem=190000`, and `--time=04:00:00`. Configure Release with `BUILD_TESTING`, MPI, LCAO, LibRI, LibComm, GreenX minimax, and the same dependency paths as the sixth-order build. Build:

```bash
cmake --build "$build" --parallel 40 --target \
  MODULE_RI_sternheimer_fd_hamiltonian_test \
  MODULE_IO_read_item_serial \
  MODULE_RI_sternheimer_fd_solver_test \
  MODULE_RI_sternheimer_delta_test \
  MODULE_RI_sternheimer_grid_diagnostics_test \
  abacus_3p
```

- [ ] **Step 4: Run all affected tests**

```bash
ctest --test-dir "$build" --output-on-failure \
  -R 'MODULE_RI_sternheimer_fd_hamiltonian_test|MODULE_IO_read_item_serial|MODULE_RI_sternheimer_fd_solver_test|MODULE_RI_sternheimer_delta_test|MODULE_RI_sternheimer_grid_diagnostics_test'
```

Expected: all five selected tests pass. Copy `abacus_3p` to a commit-named artifact directory and record its SHA256 before any physics run.

### Task 5: Run the H2 Numerical Gate

**Files:**
- Create: `.codex_tmp/run_fd8_h2_smoke_${short_commit}.slurm` (untracked)

- [ ] **Step 1: Prepare matched order-2/6/8 cases**

Copy the frozen one-frequency H2 inputs used for the sixth-order gate. Change only:

```text
sternheimer_grid_diagnostics         1
sternheimer_fd_order                 2  # then 6 and 8 in separate directories
```

Use one normal node with 40 OpenMP threads, `MKL_NUM_THREADS=1`, `OMP_PROC_BIND=spread`, and the new binary hash.

- [ ] **Step 2: Run ABACUS once for each order and validate outputs**

For each case require:

```bash
grep -q '^status success$' OUT.*/STERNHEIMER_CHI0.dat
grep -q '^all_converged yes$' OUT.*/STERNHEIMER_CHI0.dat
awk '$1 == "sternheimer_component_reconstruction_error_max" && $2 < 1e-10 {ok=1} END {exit !ok}' OUT.*/STERNHEIMER_CHI0.dat
```

The order-2 `v1_sternheimer_chi0_iq_1_ifreq_1_rank0.dat` must match the frozen baseline below relative `1e-11`. Record wall time and maximum resident set size for all three orders; require finite, nonidentical order-6 and order-8 responses.

- [ ] **Step 3: Stop on any H2 failure**

Do not submit N or N2 unless the binary hash, status, convergence, component reconstruction, and frozen order-2 checks all pass.

### Task 6: Run N Before N2

**Files:**
- Create: `.codex_tmp/run_n_fd8_grid_diag_40_50_${short_commit}.slurm` (untracked)

- [ ] **Step 1: Submit only N 40/50 Ry at order 8**

Use the existing N definition unchanged: 12 A box, `nspin=2`, `nelec=5`, `nbands=22`, TZDP-8au, PCA `1e-4`, 12 fixed GreenX frequencies, `ks_bands`, full-grid diagnostics, fixed full-Ewald Coulomb, 18 MPI ranks on 18 nodes, and 40 OpenMP threads per rank. Change only:

```text
sternheimer_fd_order                 8
ecutwfc                              40  # array member 0
ecutwfc                              50  # array member 1
```

- [ ] **Step 2: Run ABACUS and LibRPA end to end**

Require exact binary hashes, grids `96^3` and `108^3`, 12 chi0 files, six component matrices per spin/frequency definition already enforced by the sixth-order script, component reconstruction below `1e-10`, `libRPA finished successfully`, and a final `Total Sternheimer EcRPA` line.

- [ ] **Step 3: Evaluate the N stop condition**

Compute separate 40-to-50 Ry changes for total, SOS, Pulay, and Q-space hybrid components in kcal/mol. Stop before N2 if the order-8 Q-space change is not smaller than the sixth-order value or if total response exceeds `0.1 kcal/mol`.

### Task 7: Run the Conditional N2 Production Gate

**Files:**
- Create: `.codex_tmp/run_n2_fd8_grid_diag_40_50_${short_commit}.slurm` (untracked)

- [ ] **Step 1: Submit N2 only after N passes**

Use the fixed N2 definition: 16 A box, `nspin=1`, `nelec=10`, `nbands=44`, TZDP-8au, PCA `1e-4`, 12 fixed GreenX frequencies, `ks_bands`, fixed full-Ewald Coulomb, 18 ranks on 18 normal nodes, and 40 threads per rank. Run the `125^3` 40 Ry and `144^3` 50 Ry grids as two array members with `sternheimer_fd_order 8`.

- [ ] **Step 2: Validate each N2 case independently**

Require the same hash, status, convergence, file-count, reconstruction, grid, and LibRPA completion gates as for N. Record ABACUS and LibRPA wall time, peak memory, node count, and node-hours separately.

- [ ] **Step 3: Apply the production decision**

Use these exact rules:

```text
Keep order 6:
  abs(N2_total_fd8 - N2_total_fd6) < 0.02 kcal/mol
  and the N2 Q-space improvement is not operationally useful.

Recommend order 8:
  N total 40->50 < 0.1 kcal/mol
  and N2 total 40->50 < 0.1 kcal/mol
  and N2 Q hybrid 40->50 < 0.1 kcal/mol
  and matched-cutoff wall overhead relative to order 6 <= 20%.

Escalate to a spectral reference:
  order 8 does not reduce the residual Q-space sensitivity.
```

Do not infer convergence from `2E_c(N)-E_c(N2)` cancellation; gate N and N2 separately.

### Task 8: Document Results and Close the Branch

**Files:**
- Modify: `/Users/ghj/Downloads/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.tex`
- Verify: `/Users/ghj/Downloads/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project/main.pdf`

- [ ] **Step 1: Add a provenance and cost table**

Record the source commit, binary SHA256, build/test job, H2 job, N/N2 job IDs, grid dimensions, order, wall time, node count, and node-hours. Label unfinished or failed jobs explicitly and preserve their directories.

- [ ] **Step 2: Add component convergence plots**

Plot N and N2 separately against actual grid points per axis. Include total, SOS, Pulay, and Q-space hybrid 40-to-50 changes for orders 2/4/6/8. Do not use `ecutwfc` alone as the horizontal coordinate when two cutoffs map to the same FFT grid.

- [ ] **Step 3: State the bounded conclusion**

Write whether order 6 or 8 is recommended and quote both accuracy and node-hour evidence. If the Q-space gate fails, state that finite-difference order is not the sole remaining error and identify the spectral kinetic operator as the next diagnostic, not as a completed implementation.

- [ ] **Step 4: Compile and inspect the TeX document locally**

```bash
cd /Users/ghj/Downloads/同步空间/AITP_project/sternheimer_abacus/sternheimer_siab_project
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Expected: `main.pdf` is created with no fatal LaTeX errors; inspect the updated table and plots for clipping or overflow.

- [ ] **Step 5: Run final branch verification**

```bash
git diff --check
git status --short
git log -3 --format='%h %an <%ae> | %cn <%ce> | %s'
```

Expected: only intentionally untracked `.codex_tmp/` artifacts remain; committed source changes use Codex as author and AroundPeking as committer.

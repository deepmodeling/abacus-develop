# Independent Sternheimer Response Grid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decouple the Delta-Sternheimer response grid from the converged LCAO PBE grid, validate the reciprocal-space restriction on a small periodic Si system, and release a 2x2 graphene plus water feasibility pilot only after numerical gates pass.

**Architecture:** Add a compatibility-default response cutoff, construct replicated serial ABACUS fine and response FFT bases only when the requested cutoff generates a smaller grid, and spectrally restrict the gathered PBE local potential by matching integer reciprocal vectors. Reuse unchanged LCAO coefficients and rebuild sampled states, nonlocal projectors, auxiliary channels, the FD Hamiltonian, and the Delta subspace directly on the response grid. Keep same-grid execution on the existing path and validate reduced-grid physics separately.

**Tech Stack:** C++17, ABACUS `PW_Basis`, GoogleTest, MPI/OpenMP, reader-v1 response data, Slurm on `df_iopcas_ghj` or `df_dcu`, Python/NumPy response validators, LaTeX.

---

## File Map

- Modify `source/source_io/module_parameter/input_parameter.h`: store `sternheimer_response_ecutwfc` with compatibility default zero.
- Modify `source/source_io/module_parameter/read_input_item_output.cpp`: register and validate the response cutoff.
- Modify `source/source_io/test_serial/read_input_item_test.cpp`: cover the parser contract.
- Create `source/source_lcao/module_ri/sternheimer_response_grid.h`: declare the response-grid owner, cutoff decision, and spectral restriction interface.
- Create `source/source_lcao/module_ri/sternheimer_response_grid.cpp`: construct serial full-G bases and restrict fine real fields without aliasing.
- Create `source/source_lcao/module_ri/test/sternheimer_response_grid_test.cpp`: test identity, constants, retained modes, rejected modes, and mixed fields.
- Modify `source/source_lcao/module_ri/CMakeLists.txt`: compile the new production source.
- Modify `source/source_lcao/module_ri/test/CMakeLists.txt`: compile and register the focused grid test.
- Modify `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp`: select the response grid, transfer the potential, use it throughout the periodic Delta path, and report provenance.
- Create `.codex_tmp/validate_si_response_grid.py`: compare reader-v1 metadata, finite values, residuals, and relative Frobenius differences; keep untracked.
- Create `.codex_tmp/run_si_response_grid_${short_commit}.slurm`: execute the cost-ordered Si gate on `df_iopcas_ghj`; keep untracked.
- Create `.codex_tmp/run_graphene2x2_response_grid_${short_commit}.slurm`: execute only the post-Si feasibility gate; keep untracked.
- Modify `/Users/ghj/Downloads/同步空间/AITP_project/h2o_bn_rpa_adsorption/main.tex`: record accepted code, grids, comparisons, timings, and the next physical gate.

### Task 1: Add the Input Contract

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h:591-603`
- Modify: `source/source_io/module_parameter/read_input_item_output.cpp:1131-1150`
- Test: `source/source_io/test_serial/read_input_item_test.cpp:1643-1665`

- [ ] **Step 1: Write the failing parser test**

Insert after the `sternheimer_fd_order` block:

```cpp
{ // sternheimer_response_ecutwfc
    auto it = find_label("sternheimer_response_ecutwfc", readinput.input_lists);
    ASSERT_NE(it, readinput.input_lists.end());
    EXPECT_DOUBLE_EQ(param.input.sternheimer_response_ecutwfc, 0.0);

    it->second.str_values = {"30.0"};
    it->second.read_value(it->second, param);
    EXPECT_DOUBLE_EQ(param.input.sternheimer_response_ecutwfc, 30.0);

    param.input.sternheimer_response_ecutwfc = -1.0;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(it->second.check_value(it->second, param), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("must be zero or positive"));
}
```

- [ ] **Step 2: Run the focused test and verify the red state**

Run on the configured ABACUS test build:

```bash
cmake --build build --parallel 16 --target MODULE_IO_read_item_serial
ctest --test-dir build --output-on-failure -R '^MODULE_IO_read_item_serial$'
```

Expected: compilation fails because the parameter and label do not exist.

- [ ] **Step 3: Add the parameter and parser item**

Add to `Input_para`:

```cpp
double sternheimer_response_ecutwfc = 0.0; ///< response-grid cutoff in Ry; 0 reuses the PBE grid
```

Register an `Input_Item` with type `Real`, unit `Ry`, default `0.0`, and this check:

```cpp
if (para.input.sternheimer_response_ecutwfc < 0.0)
{
    ModuleBase::WARNING_QUIT("ReadInput", item.label + " must be zero or positive.");
}
```

Defer the cross-field `response <= ecutwfc` check to the response-grid selector, after `ecutwfc` reset logic has completed.

- [ ] **Step 4: Run the parser test to green**

Run the command from Step 2. Expected: the selected test passes.

- [ ] **Step 5: Commit the input contract**

```bash
git add source/source_io/module_parameter/input_parameter.h \
  source/source_io/module_parameter/read_input_item_output.cpp \
  source/source_io/test_serial/read_input_item_test.cpp
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'feat(sternheimer): add response grid cutoff input'
```

### Task 2: Specify Spectral Restriction with Failing Tests

**Files:**
- Create: `source/source_lcao/module_ri/test/sternheimer_response_grid_test.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Add a reusable test-basis builder**

Create serial full-complex bases with explicit dimensions:

```cpp
std::unique_ptr<ModulePW::PW_Basis> make_basis(const int nx, const int ny, const int nz)
{
    auto basis = std::make_unique<ModulePW::PW_Basis>("cpu", "double");
#ifdef __MPI
    basis->initmpi(1, 0, MPI_COMM_SELF);
#endif
    basis->initgrids(1.0, ModuleBase::Matrix3(8.0, 0.0, 0.0,
                                             0.0, 8.0, 0.0,
                                             0.0, 0.0, 8.0), nx, ny, nz);
    basis->initparameters(false, 1.0e12);
    basis->fft_bundle.initfftmode(0);
    basis->setuptransform();
    basis->collect_local_pw();
    return basis;
}
```

- [ ] **Step 2: Add identity and constant tests**

Require `restrict_sternheimer_real_field(fine, fine, field)` to return the input exactly through the identity branch. Require an `8^3 -> 4^3` restriction of the constant `3.25` to have maximum absolute error below `1e-12`.

- [ ] **Step 3: Add retained and rejected plane-wave tests**

Generate real fields as `cos(2*pi*(mx*x/nx + my*y/ny + mz*z/nz))`. Require mode `(1,1,0)` to survive `8^3 -> 4^3` with relative L2 error below `1e-11`; require mode `(3,0,0)` to have output norm below `1e-11` of its fine-grid norm.

- [ ] **Step 4: Add the no-alias mixed-field test**

Restrict `low + 0.3*high`, subtract the independently restricted `low`, and require the remaining output norm below `1e-11` of the restricted low-mode norm.

- [ ] **Step 5: Add cutoff-decision tests**

Require:

```cpp
EXPECT_FALSE(sternheimer_uses_independent_response_grid(0.0, 80.0));
EXPECT_TRUE(sternheimer_uses_independent_response_grid(30.0, 80.0));
EXPECT_THROW(sternheimer_uses_independent_response_grid(-1.0, 80.0), std::invalid_argument);
EXPECT_THROW(sternheimer_uses_independent_response_grid(90.0, 80.0), std::invalid_argument);
```

- [ ] **Step 6: Register and run the red test**

Add an executable linked with `planewave`, `parameter`, `base`, `device`, `container`, and `${math_libs}`. Run:

```bash
cmake --build build --parallel 16 --target MODULE_RI_sternheimer_response_grid_test
ctest --test-dir build --output-on-failure -R '^MODULE_RI_sternheimer_response_grid_test$'
```

Expected: compilation fails because the response-grid API does not exist.

### Task 3: Implement the Response Grid and Fourier Projection

**Files:**
- Create: `source/source_lcao/module_ri/sternheimer_response_grid.h`
- Create: `source/source_lcao/module_ri/sternheimer_response_grid.cpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`

- [ ] **Step 1: Define the owner and interface**

Use a non-copyable owner so the FFT basis lifetime covers every sampled response object:

```cpp
struct SternheimerResponseGrid
{
    std::unique_ptr<ModulePW::PW_Basis> serial_fine_basis;
    std::unique_ptr<ModulePW::PW_Basis> serial_response_basis;
    const ModulePW::PW_Basis* basis = nullptr;
    bool independent = false;
    double requested_ecutwfc = 0.0;
};

bool sternheimer_uses_independent_response_grid(double response_ecutwfc, double pbe_ecutwfc);
SternheimerResponseGrid make_sternheimer_response_grid(const ModulePW::PW_Basis& pbe_basis,
                                                       double pbe_ecutwfc,
                                                       double response_ecutwfc,
                                                       int fft_mode);
std::vector<double> restrict_sternheimer_real_field(const ModulePW::PW_Basis& fine_basis,
                                                    const ModulePW::PW_Basis& coarse_basis,
                                                    const std::vector<double>& fine_values);
```

- [ ] **Step 2: Implement cutoff and identity decisions**

Reject a negative response cutoff, a nonpositive PBE cutoff, and response cutoffs above PBE by more than `1e-12 * max(1, pbe)`. Return the PBE basis for cutoff zero. Build the candidate response basis from `4.0 * response_ecutwfc`; if all three dimensions equal the PBE grid, discard it and return the PBE path.

- [ ] **Step 3: Build serial complex-G fine and response bases**

Initialize both bases with `MPI_COMM_SELF` under MPI.  Build the fine basis with
the exact `pbe_basis.nx/ny/nz` dimensions and `initparameters(false,
4.0 * pbe_ecutwfc)`.  Build the response basis with `initgrids(lat0, latvec,
4.0 * response_ecutwfc)` and `initparameters(false,
4.0 * response_ecutwfc)`.  For both, set the FFT mode, call
`setuptransform()`, then `collect_local_pw()`.  Require `nrxx == nxyz` and one
local rank.  Do not enable the OFDFT `full_pw` grid-shape option because it can
change otherwise identical PBE dimensions.

- [ ] **Step 4: Implement reciprocal-index matching**

Transform the full fine real vector with `real2recip`. Build a map keyed by rounded integer `(gx,gy,gz)` from the fine `gdirect` values. For every coarse G vector, copy the matching coefficient or leave zero. Reject any `gdirect` component farther than `1e-10` from an integer. Inverse transform into a full real vector on the coarse basis.

- [ ] **Step 5: Run all response-grid tests**

Run the command from Task 2 Step 6. Expected: all focused tests pass.

- [ ] **Step 6: Commit the grid component**

```bash
git add source/source_lcao/module_ri/sternheimer_response_grid.* \
  source/source_lcao/module_ri/CMakeLists.txt \
  source/source_lcao/module_ri/test/sternheimer_response_grid_test.cpp \
  source/source_lcao/module_ri/test/CMakeLists.txt
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'feat(sternheimer): add independent response grid'
```

### Task 4: Integrate the Grid into Periodic Delta-Sternheimer

**Files:**
- Modify: `source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp:2400-3250,4290-4345`

- [ ] **Step 1: Write a same-grid integration assertion**

Before changing production flow, add a focused helper test that constructs a requested cutoff yielding the PBE dimensions and requires `independent == false` and `basis == &pbe_basis`.

- [ ] **Step 2: Select the response basis once**

At periodic Delta setup, construct `SternheimerResponseGrid response_grid`. Use `*response_grid.basis` for `make_sternheimer_fd_full_grid`, all LCAO/ABFS sampling, nonlocal projector construction, and the Hamiltonian.

- [ ] **Step 3: Transfer the local potential only when needed**

Gather the complete PBE potential once. If `response_grid.independent`, call `restrict_sternheimer_real_field`; otherwise preserve the existing potential path. Store a response-grid-sized vector used by every k-pair Hamiltonian.

- [ ] **Step 4: Keep unchanged physics inputs unchanged**

Do not alter LCAO eigenvalues/coefficients, occupations, k/q pair construction, auxiliary radial functions, frequency weights, FD order, solver/preconditioner, full 2D Ewald settings, or reader-v1 assembly.

- [ ] **Step 5: Add exact provenance**

Write PBE cutoff, requested/effective response cutoff, source, dimensions, and point counts to status output. Add grid dimensions to `full_grid_ready` details.

- [ ] **Step 6: Run focused and neighboring regression tests**

```bash
cmake --build build --parallel 16 --target \
  MODULE_IO_read_item_serial \
  MODULE_RI_sternheimer_response_grid_test \
  MODULE_RI_sternheimer_abacus_fd_adapter_test \
  MODULE_RI_sternheimer_delta_test \
  MODULE_RI_sternheimer_kq_test
ctest --test-dir build --output-on-failure \
  -R 'MODULE_IO_read_item_serial|MODULE_RI_sternheimer_(response_grid|abacus_fd_adapter|delta|kq)_test'
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit integration**

```bash
git add source/source_lcao/module_ri/sternheimer_abacus_st_smoke.cpp
GIT_AUTHOR_NAME=Codex GIT_AUTHOR_EMAIL=codex@openai.com \
GIT_COMMITTER_NAME=AroundPeking GIT_COMMITTER_EMAIL=gonghuanjing@iphy.ac.cn \
git commit -m 'feat(sternheimer): use independent periodic response grid'
```

### Task 5: Build and Run the Small-System Numerical Gate

**Files:**
- Create: `.codex_tmp/validate_si_response_grid.py` (untracked)
- Create: `.codex_tmp/run_si_response_grid_${short_commit}.slurm` (untracked)

- [ ] **Step 1: Build on `df_iopcas_ghj`**

Use a clean staged source archive from the committed branch. Require `ENABLE_MPI=ON`, `ENABLE_LCAO=ON`, `ENABLE_LIBRI=ON`, `ENABLE_LIBCOMM=ON`, `DEBUG_INFO=ON`, and an `abacus_3p` executable. Record source commit and executable SHA256.

- [ ] **Step 2: Run the same-grid compatibility case**

Run Si q001 at fixed one frequency once with `sternheimer_response_ecutwfc=0` and once with an explicit cutoff that generates the PBE dimensions. Require identical metadata and relative Frobenius response difference `<=1e-10`.

- [ ] **Step 3: Run a cost-ordered reduced-grid scan**

Keep the PBE cutoff fixed. Try response cutoffs from the finest candidate downward and stop at the first failed numerical gate. Run reduced channels first, then the complete Si channel set for the selected cutoff.

- [ ] **Step 4: Validate the response artifacts**

The validator must reject missing files, mismatched dimensions/q/frequency/weights, NaN/Inf values, any `converged=no`, or relative Frobenius difference above `1e-6`. Report maximum equation residual, wall time, node-hours, and MaxRSS.

- [ ] **Step 5: Record the decision**

If no reduced grid passes, stop and retain the PBE grid. If one passes, record its actual dimensions and promote it only to the 2x2 feasibility pilot.

### Task 6: Commit, Push, and Gate the 2x2 Pilot

**Files:**
- Create: `.codex_tmp/run_graphene2x2_response_grid_${short_commit}.slurm` (untracked)
- Modify: `/Users/ghj/Downloads/同步空间/AITP_project/h2o_bn_rpa_adsorption/main.tex`

- [ ] **Step 1: Run final source verification**

```bash
git diff --check
git status --short
git log -3 --format='%h %an <%ae> | %cn <%ce> | %s'
```

Require every new commit to show Codex as author and AroundPeking as committer.

- [ ] **Step 2: Push to `master_ghj` only after code and same-grid gates pass**

Fetch the remote branch, confirm it has not advanced incompatibly, then push the verified commit range. Do not force-push.

- [ ] **Step 3: Submit one duplicate-guarded 2x2 pilot only after Si passes**

Hash the executable and all inputs; search active and historical jobs for the same fingerprint. Submit only when no equivalent active or completed validated run exists. Use one q, one fixed frequency, and the full physical state/channel definitions; do not call it an adsorption-energy calculation.

- [ ] **Step 4: Apply feasibility gates**

Require `channels_ready`, `delta_subspace_ready`, two finite converged equation batches, no OOM, and a conservative runtime estimate. Do not release the remaining q points automatically.

- [ ] **Step 5: Update and compile the research note**

Record the exact commit, binary hash, PBE/response grids, Si matrix differences, timing/memory, pilot job ID, and gate status. Compile the PDF and verify the output before reporting delivery.

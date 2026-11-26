# PPCG validation

The CPU timing study below uses the explicitly identified historical revision. Native CUDA implementation checks are recorded separately at the end of this report.

## CPU configuration

Validated on 2026-11-24 with GNU 13.3 Release, OpenMPI, OpenMP, OpenBLAS, and an Intel Xeon Gold 6242. LCAO and CUDA were disabled. Production sources are from `dd385365908b088333030b4623d72e3533b521a6`; the accompanying test-data/documentation additions do not change the executable.

Each application ran with two MPI ranks, core binding, and `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Residual tracing was unset during comparisons. Executable SHA-256: `c649f8a557e48527d0712f79c6dec8fe9260ffe0807dfbe356e6520bf4457dac`.

The comparison uses fixtures `817_PW_PPCG`, `818_PW_PPCG_Si`, `819_PW_PPCG_Al`, `003_PW_UPF100_USPP_Fe`, `008_PW_UPF201_USPP_NaCl`, and `009_PW_UPF201_USPP`. Geometry, pseudopotentials, cutoffs, k meshes, smearing, and mixing settings are held fixed within each case. All solvers use `pw_seed=1`, `scf_thr_type=1`, `scf_thr=1e-10`, `scf_nmax=200`, and `pw_diag_nmax=200`. The initial `pw_diag_thr` is omitted, exercising its default of `0.01`; PPCG uses the default `ppcg_rr_step=1` and block size 4.

## Default diagonalization threshold

PPCG interprets the threshold as an absolute generalized residual norm, whereas the comparison solvers use different stopping rules. Equal parameter values therefore do not imply equal error or work. SCF automatically tightens the diagonalization threshold after the initial electronic iteration.

The following differences compare default-threshold PPCG against CG with initial `pw_diag_thr=1e-10`, at the same SCF settings. CG and PPCG converged in every case. Values are maximum absolute differences in total energy, a force component, and a stress component over the three measured PPCG runs. These are finite-SCF comparisons, not certified exact errors.

| System | Energy (eV) | Force (eV/Å) | Stress (kbar) |
|---|---:|---:|---:|
| GaAs (NCPP) | 1.073e-10 | 5.301e-04 | 4.067e-03 |
| Si (NCPP) | 9.464e-10 | 0.000e+00 | 5.026e-03 |
| Al (NCPP) | 1.805e-10 | 0.000e+00 | 4.786e-02 |
| Fe (USPP) | 1.219e-09 | 0.000e+00 | 2.871e-02 |
| NaCl (USPP) | 1.473e-08 | 1.132e-04 | 3.385e-02 |
| HCl (mixed NCPP/USPP) | 7.080e-10 | 2.715e-06 | 2.262e-02 |

To isolate the initial PPCG tolerance, the same calculations were repeated changing only `pw_diag_thr` from its default to `1e-10`:

| System | Energy change (eV) | Force change (eV/Å) | Stress change (kbar) |
|---|---:|---:|---:|
| GaAs (NCPP) | 2.547e-11 | 4.102e-07 | 8.075e-06 |
| Si (NCPP) | 5.969e-13 | 0.000e+00 | 1.143e-04 |
| Al (NCPP) | 8.868e-12 | 0.000e+00 | 1.739e-06 |
| Fe (USPP) | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| NaCl (USPP) | 1.228e-11 | 1.677e-05 | 2.650e-04 |
| HCl (mixed NCPP/USPP) | 5.156e-10 | 1.533e-07 | 2.290e-02 |

All 24 PPCG runs in the controlled default-threshold batch, including warm-ups, and all six tighter-initial-threshold runs converged without PPCG nonconvergence warnings. These results support the conservative residual interpretation on the tested systems. Force and stress accuracy also depends on SCF convergence; tightening only the initial eigensolver threshold does not remove that error. System-specific convergence checks remain necessary for more demanding property tolerances.

## Efficiency

Runs were sequential: one warm-up followed by three measured runs for every case/solver, rotating the order of CG, Davidson, Davidson-subspace, and PPCG. The complete batch contains 96 launches (72 measured). Each cell gives median wall time / median `HSolverPW::solve_psik` time, in seconds. The wall timer surrounds the MPI subprocess.

| System | CG | Davidson | Davidson-subspace | PPCG |
|---|---:|---:|---:|---:|
| GaAs (NCPP) | 1.969 / 0.690 | 1.719 / 0.510 | 1.719 / 0.500 | 5.331 / 4.140 |
| Si (NCPP) | 1.519 / 0.470 | 1.319 / 0.310 | 1.318 / 0.290 | 2.572 / 1.570 |
| Al (NCPP) | 1.569 / 0.680 | 1.218 / 0.370 | 1.118 / 0.290 | 2.171 / 1.370 |
| Fe (USPP) | 1.518 / 0.630 | failed | mismatch¹ | 1.719 / 0.830 |
| NaCl (USPP) | 1.368 / 0.060 | failed | 3.173 / 0.300 | 1.469 / 0.200 |
| HCl (mixed NCPP/USPP) | 1.218 / 0.270 | 1.218 / 0.280 | 1.168 / 0.240 | 1.870 / 0.940 |

¹ Davidson-subspace reported SCF convergence for Fe but its energy differed from the tightened CG result by about 0.881 eV. That result is excluded from valid efficiency/accuracy comparisons. Ordinary Davidson exited on Fe and NaCl with a linearly dependent initial-basis diagnostic (all three measured attempts for each case). The failed runs are not included in timings. No changes to those solvers are included here.

PPCG wall-time ranges for GaAs, Si, Al, Fe, NaCl, and HCl were respectively 5.329–5.681, 2.572–2.672, 2.170–2.171, 1.718–1.770, 1.418–1.469, 1.868–1.920 seconds. PPCG is slower than CG on these small systems. No general speedup, large-system scaling, or GPU performance claim is made.

## Independent USPP regression references

Autotest cases `821_PW_PPCG_USPP_Fe`, `822_PW_PPCG_USPP_NaCl`, and `823_PW_PPCG_USPP_HCl` retain the structures/pseudopotentials/k meshes of cases 003, 008, and 009. Their SCF threshold is tightened to `1e-16`, with `scf_nmax=200` and `pw_diag_nmax=200`; the initial PPCG diagonalization threshold and RR interval retain their defaults. Their reference files were generated independently using BPCG with these same numerical settings. Existing fixtures and their references are unchanged. Autotest tolerances are not relaxed.

| System | PPCG–BPCG energy (eV) | Max force difference (eV/Å) | Max stress difference (kbar) |
|---|---:|---:|---:|
| Fe (USPP) | 2.274e-12 | 0.000e+00 | 1.164e-04 |
| NaCl (USPP) | 0.000e+00 | 3.839e-07 | 1.921e-04 |
| HCl (mixed NCPP/USPP) | 2.842e-13 | 3.300e-09 | 9.088e-05 |

At this deliberately tight SCF setting, Fe and HCl retain some iteration-budget warnings: residuals can stagnate near the numerical floor of `1e-14` (maximum warning residual in the two-rank calibration `2.11e-14`). SCF continues under the documented status policy. The warnings are not suppressed or relabeled as certified eigenpair convergence; the independent observable checks pass. A separate `scf_thr=1e-14`, `pw_diag_thr=1e-12` CG probe exited with `psi_norm <= 0` on four systems and was not used as a reference. Changing the PW SCF metric to type 2 caused HCl SCF stagnation across solvers; the reported comparison keeps the standard PW type-1 metric.

## Reproduction and checks

Build and unit checks:

```bash
cmake -S . -B build-review-20261124 -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DENABLE_LCAO=OFF -DBUILD_TESTING=ON
cmake --build build-review-20261124 -j8 --target abacus_pw_para \
  MODULE_HSOLVER_ppcg MODULE_HSOLVER_ppcg_float MODULE_HSOLVER_pw \
  MODULE_HSOLVER_ppcg_parallel MODULE_IO_ppcg_input
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ctest --test-dir build-review-20261124 -V \
  -R '^(MODULE_HSOLVER_(pw|ppcg|ppcg_float|ppcg_parallel|ppcg_parallel_test)|MODULE_IO_ppcg_input)$'
```

MPI CTest passed 6/6: 54 double-precision, four float, eight PW, and five INPUT cases, plus MPI runs with 2/3/6 ranks and two pools. A separate `ENABLE_MPI=OFF` build passed its three PPCG/float/INPUT CTest entries. INPUT checks in the CPU-only build cover NSCF/GPU rejection and RR interval validation. CLI `--check-input`, help and version checks passed for the supported inputs; NSCF was rejected during both checking and actual startup.

For each timing calculation, make a fresh copy of a listed fixture, apply the common overrides above, set only `ks_solver` differently within each group, resolve `pseudo_dir`, and run:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  mpirun --bind-to core --map-by core -np 2 /absolute/path/to/abacus_pw_para
```

For native regression checks, select cases 817/818/819/821/822/823 in a temporary case list and run from `tests/01_PW`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ABACUS_MPIEXEC_PREFLAGS='--bind-to core --map-by core' \
  bash ../integrate/Autotest.sh -a /absolute/path/to/abacus_pw_para \
  -n 2 -o 1 -j 1 -f /absolute/path/to/selected_cases.txt
```

The six native cases passed 36/36 checks with two ranks; the three new USPP cases also passed 15/15 checks with four ranks. A separate traced GaAs run preserved all 16 solves across eight SCF iterations and two k points, with initial and final records for each solve. A Si run with `ppcg_rr_step=4` converged and produced block-sweep records, while interval 0 was rejected.

Local raw logs, INPUT hashes, timing samples, failed probes, and JSON summaries are archived under `/mnt/data/Cheerly/ppcg-review-20261124/`. The tables above identify the usable data and do not mix failed probes into successful timing or accuracy samples.

## Native CUDA validation

The native CUDA implementation in `99950ee5a` was validated on 2026-11-24 with an NVIDIA GeForce RTX 3090 (24 GB), GNU 13.3, CUDA 13.5, OpenMPI and OpenBLAS. The Release build enabled MPI/OpenMP/testing and disabled LCAO. Runtime checks used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

Wavefunctions, residuals, search directions and H/S images remain in device memory. CUDA kernels compute residuals, preconditioning, column products and scaling; device BLAS performs Gram products and rotations. Only projected matrices, their coefficients and per-band convergence data cross to the host for LAPACK and the existing MPI reduction wrappers. The preconditioner is uploaded once per diagonalization. There is no full-wavefunction CPU callback bridge. ROCm is not implemented.

The CUDA CTest target passed all six cases: standard/generalized complex double, infrequent RR with a locked prefix, standard/generalized complex float, and iteration-budget exhaustion. The generalized fixture has a known spectrum, complex dense Hamiltonian and nonuniform overlap. Checks independently evaluate eigenvalues, residuals, S-orthogonality and leading-dimension padding. Every H/S callback verifies that both operand pointers are device allocations. Compute Sanitizer reported **0 errors and 0 leaked bytes** across all six cases.

### Application regression

The six `tests/11_PW_GPU/scf_ppcg_*` Autotest fixtures use the CPU fixtures' geometry, pseudopotentials, references and acceptance tolerances. Only device selection and explicit double precision differ. All six SCF calculations converged and **36/36 checks passed** with two MPI ranks. The current GPU PW setup assigns one k-point pool per rank; this is k-point parallelism, not a distributed plane-wave-row GPU test.

| Fixture | SCF iterations | Final density deviation | PPCG iteration-limit warnings | Largest warning residual |
|---|---:|---:|---:|---:|
| GaAs | 8 | 6.06555e-7 | 0 | — |
| Si | 7 | 1.44325e-7 | 0 | — |
| Al | 14 | 7.91835e-11 | 0 | — |
| Fe (USPP) | 26 | 6.09993e-17 | 9 | 1.63674e-14 |
| NaCl (USPP) | 35 | 4.70540e-17 | 0 | — |
| HCl (mixed NCPP/USPP) | 20 | 3.08266e-17 | 24 | 2.21592e-14 |

Warning counts above include stderr from both MPI ranks, rather than only the world-root running log.

In the CUDA CI run of `11e39d8f3` on an NVIDIA L20 with CUDA 12.9, all six fixtures also passed their reference-observable checks. Fe emitted 12 iteration-limit warnings with a maximum residual of `4.81339e-9`; HCl emitted 24 with a maximum of `2.15364e-14`. These individual diagonalization calls did not satisfy their requested thresholds. Observable agreement does not certify those eigenpairs, and the larger Fe warnings are not classified as roundoff-level convergence. The local RTX 3090 table and these CI observations are reported separately.

Those historical strict USPP runs retained eigensolver warnings; passing observable checks did not certify those eigenpairs. The residual threshold was not established to be unattainable. An initial batch exceeded the local harness's 300-second aggregate limit. A fresh complete batch passed with a 900-second harness limit; numerical thresholds, iteration budgets and references were unchanged.

A separate GaAs calculation with `ppcg_rr_step=4`, two ranks and two k-point pools converged without PPCG warnings. The world-root trace preserved all eight of its solves across eight SCF iterations, including initial/final records and block sweeps. Other pools intentionally do not write that shared trace. GPU CLI help/version and valid `--check-input` checks passed; RR interval 0 and NSCF were rejected, with NSCF also rejected during actual startup. A four-rank probe was rejected by the existing k-point setup because this fixture has only two k points; the successful probe uses two ranks.

After the device-workspace refactor, the CPU focused MPI CTest suite passed **6/6** and the six CPU application regressions passed **36/36** again. These are correctness regressions; the earlier CPU timing table was not remeasured.

### Strict-residual refinement (2026-11-26)

When a successful Ritz solve stalls within four stopping tolerances, PPCG now
cycles through retained pairs after five non-improving iterations. This also
reduces errors in already converged neighbors that can constrain the remaining
pair. The correction is projected against the retained states and the
candidate is S-normalized. It is accepted only if a fresh H/S application
confirms a smaller residual, without a resolved energy increase or reordering.
A rejected candidate is restored. Each attempt consumes an iteration from the
existing budget; the convergence criterion and final fresh-residual check are
unchanged. The `residual_refinement` trace stage records these attempts.

The six original CPU fixtures and six original RTX 3090 CUDA fixtures passed
**72/72 reference checks** with **zero PPCG nonconvergence warnings across both
MPI ranks**, including Fe and HCl on both backends. Their inputs, references and tolerances were
unchanged; Fe retained `scf_thr=1e-16` and `pw_diag_nmax=200`.

Both Fe calculations converged in 26 SCF steps, with total energy approximately
`-2739.33531173184 eV` and final density deviation approximately `6.10e-17`. Its energy,
force and stress checks passed against the existing independent reference.
These local results do not by themselves reproduce or certify the older L20
residual rebound described above; the updated branch requires fresh CI evidence.

The CPU suite passed 61 tests in MPI and non-MPI builds. Seven new tests cover
standard/generalized correction geometry, refinement of a converged neighbor,
preservation of converged pairs, rejection of nonfinite directions, and exact
float/double cancellation remainders. The CUDA target passed eight cases,
including the two cancellation checks. Both residual kernels now use fused
multiply-add: separate product rounding could otherwise turn a nonzero residual
into zero. The MPI regression driver also passed its multi-rank configurations. The correction uses the existing device workspaces,
vector kernels and H/S interface; no extended-precision solver or host
wavefunction bridge was added. Local logs and executable hashes are retained in
`/mnt/data/Cheerly/ppcg-fe-stability-20261126/`.

### Transfer audit and reproduction

Nsight Systems traced `DiagoPPCGGPU.GeneralizedComplexDouble` (24 active rows, leading dimension 31, five bands). Its full wavefunction occupies 2,480 bytes. Exactly one upload initializes the fixture and one download checks its final result. Between those transfers, the solver made 680 host-to-device transfers (128,480 bytes total) and 896 device-to-host transfers (181,352 bytes total), with a maximum of **800 bytes per transfer**. No complete wavefunction was transferred to or from the host during the solve. Device-to-device copies are excluded from these host-transfer counts.

The GPU was shared with another workload, so these checks do not establish speedup. Small projected solves, synchronization and allocations still have overhead; no GPU performance advantage is claimed.

```bash
cmake -S . -B build-ppcg-gpu-20261124 -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DENABLE_LCAO=OFF \
  -DBUILD_TESTING=ON -DUSE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
  -DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include/cccl \
  -DCMAKE_CUDA_FLAGS=-I/usr/local/cuda/include/cccl
cmake --build build-ppcg-gpu-20261124 -j8 \
  --target abacus_pw_gpu MODULE_HSOLVER_ppcg_gpu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ctest --test-dir build-ppcg-gpu-20261124 -V -R '^MODULE_HSOLVER_ppcg_gpu$'
```

The extra CCCL include flags apply to this local CUDA 13.5 installation; they are not a repository-wide toolchain change. Run Compute Sanitizer with `--tool memcheck --leak-check full --error-exitcode 1` on the GPU test executable. Profile it with `nsys profile --trace=cuda --sample=none --cpuctxsw=none`, selecting `--gtest_filter=DiagoPPCGGPU.GeneralizedComplexDouble`, and export the report to SQLite to inspect `CUPTI_ACTIVITY_KIND_MEMCPY`.

For application checks, create a temporary case list containing the six `scf_ppcg_*` fixtures and run from `tests/11_PW_GPU`:

```bash
ABACUS_MPIEXEC_PREFLAGS='--bind-to core --map-by core' \
  bash ../integrate/Autotest.sh -a /absolute/path/to/abacus_pw_gpu \
  -n 2 -o 1 -j 1 -f /absolute/path/to/selected_gpu_cases.txt
```

The CUDA CI workflow now explicitly runs the GPU unit target, and its existing `11_PW_GPU` run includes all six new application fixtures. Local logs, the transfer audit, executable hashes and the failed-probe diagnostics are archived in `/mnt/data/Cheerly/ppcg-review-20261124/gpu-port/`.

## Formatting follow-up

On 2026-11-25, a formatting-only cleanup applied the repository's configured tools to new files and modified regions. All 16 new C++/CUDA source/header files passed `clang-format 23.1.3 --dry-run --Werror`. The MPI shell helper uses the exact shfmt output from CI and passes `bash -n`. CPU and CUDA applications were rebuilt; focused CPU MPI CTest passed 6/6, non-MPI CTest passed 3/3, and the GPU target passed all six cases again. Numerical settings, references and iteration budgets were unchanged.

The repository-wide pre-commit job is advisory (`continue-on-error: true`). Its raw clang-format, cmake-format and shfmt hooks failed even though GitHub displayed a successful job. The log included formatting changes in 2,475 files outside this PR and unsupported JSON-language configuration errors. A green job result therefore does not certify repository-wide formatting; the focused checks above describe this PR's validation. Existing files were formatted only in the modified regions to avoid an unrelated repository-wide rewrite.

A subsequent raw CI formatter audit of `71e6ab31a` changed no newly added PR file. Repository-wide formatter failures remained in existing files; the non-blocking workflow status is still not a full formatting pass.

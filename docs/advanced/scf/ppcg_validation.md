# CPU PPCG validation

## Configuration

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

MPI CTest passed 6/6: 54 double-precision, four float, eight PW, and five INPUT cases, plus MPI runs with 2/3/6 ranks and two pools. A separate `ENABLE_MPI=OFF` build passed its three PPCG/float/INPUT CTest entries. INPUT checks cover NSCF rejection, CPU-only support, and RR interval validation. CLI `--check-input`, help and version checks passed for the supported inputs; NSCF was rejected during both checking and actual startup.

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

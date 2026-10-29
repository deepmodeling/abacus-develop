# Boundary regression verification

Recorded on 2026-10-29 for PR #8229.

## Executable and configuration

- Fixed source patch: `37dea1ae4ba682d476e15e67d2df54715770687c`; compiled before committing
  that patch on top of `7b3a835bc70a5f1c7db9ac1fea0abf149cda6f18`.
  `abacus_std_para --version` reports `v3.11.0-beta11`.
- Fixed executable SHA-256: `2dde0853955017f62dca9ddb39611367f894028b701678dba0ad4f02a9680bca`.
- Negative control: unmodified develop
  `e5ad53046d4f4c4502bd04657b7fbeba60d97ad4`,
  executable SHA-256 `aaff3e368eb0d6b5eb8630af53c97f82598978e3b26c77a5bbb985cf256dc1f9`.
- CPU double-precision LCAO; Intel icpx 2022.1, Intel MPI 2021.6,
  MKL 2022.1, ELPA 2021.11.002, Libxc 6.2.2 and LibRI enabled.
  MPI, OpenMP, testing and float FFTW were enabled; CUDA was disabled.
- Compilation used at most 20 jobs. Runtime defaults were
  `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_PIN=0`.
  Explicit OpenMP tests used four threads, at most four MPI ranks, and
  `MKL_NUM_THREADS=1`; builds and runtime tests did not overlap.

The existing configured build was updated with:

```bash
module load abacus
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_PIN=0
cmake --build build-fix --target MODULE_MD_base_mpi MODULE_MD_func \
  MODULE_MD_verlet MODULE_MD_dump_mpi MODULE_MD_run \
  MODULE_CELL_NEIGHBOR_mdcell_migrate_mpi abacus_std_para --parallel 20
ctest --test-dir build-fix -V \
  -R '^MODULE_MD_(base_mpi(_np4|_np[14]_omp4)?|func|verlet|dump_mpi(_np4)?|run)$|^MODULE_CELL_NEIGHBOR_(mdcell_migrate_mpi|domain_decomposition_np[14])$'
```

All 12 CTests passed. The new tests failed with the old boundary implementation:
the shared Psi cursor changed, and hybrid H/S disagreed with fresh folding.
The fixed tests cover Gamma and three active non-Gamma k points, both crossing
directions, nonorthogonal lattice vectors, distributed and gathered H/S, and
1/4 MPI ranks with 1/4 OpenMP threads. The hybrid oracle folds real-space
atom pairs again after shifting their lattice indices, using the saved field
endpoint; a distinct next endpoint ensures this is a time-sensitive oracle.

## Integration commands and results

Commands were run on fresh copies of tracked case files, with shared PP_ORB
and integrate resources. Here `fixed` is the absolute executable path and
`repo` is the checkout root.

```bash
cd "$repo/tests/05_rtTDDFT"
bash ../integrate/Autotest.sh -a "$fixed" -n 4 -o 1 -j 1 -f CASES_CPU.txt
# Also run cases 06 and 10 for each (ranks, threads) in (1,1), (1,4), (4,1), (4,4):
bash ../integrate/Autotest.sh -a "$fixed" -n "$ranks" -o "$threads" -j 1 \
  -f CASES_CPU.txt -r '^(06_NO_dir_TDDFT|10_NO_MD_boundary_shift)$'
```

All 28 RT cases passed (113 property checks), as did the additional 06/10
rank/thread combinations. Existing references were unchanged.

The previous CI failures were rerun with the same harness using these lists:

```bash
cd "$repo/tests/01_PW"
printf '%s\n' 093_PW_MSST 094_PW_MSST2 md_nvt 102_PW_MD_1O 103_PW_MD_2O > /tmp/pr8229-pw.list
bash ../integrate/Autotest.sh -a "$fixed" -n 4 -o 1 -j 1 -f /tmp/pr8229-pw.list
cd "$repo/tests/03_NAO_multik"
printf '%s\n' md_out_wf md_chg_extra > /tmp/pr8229-nao.list
bash ../integrate/Autotest.sh -a "$fixed" -n 2 -o 1 -j 1 -f /tmp/pr8229-nao.list
```

All seven cases passed. CUDA runtime was not tested locally.

## Case 06 reference and negative control

The replacement reference was introduced in `d373187`, using the boundary
implementation in `24fad25`. It was reproduced from the fixed source above,
with the committed INPUT/STRU/KPT, by running this in a scratch copy:

```bash
bash ../integrate/Autotest.sh -a "$fixed" -n 4 -o 1 -j 1 \
  -f CASES_CPU.txt -r '^06_NO_dir_TDDFT$' -g
# In a separate scratch copy, with the committed reference and unmodified develop:
bash ../integrate/Autotest.sh -a "$develop" -n 4 -o 1 -j 1 \
  -f CASES_CPU.txt -r '^06_NO_dir_TDDFT$'
```

The standard property collector generates energy/per-atom energy, force sum,
stress sum and timing. No custom validator or tolerance change is used.
Timing is hardware dependent and is not a numerical pass criterion.

| Run | Final energy (eV) | Force sum | Stress sum |
| --- | ---: | ---: | ---: |
| Committed reference | -30.89742338434547 | 1.490274 | 1.942380 |
| Fixed regeneration, 4 ranks / 1 thread | -30.89742338437829 | 1.490274 | 1.942380 |
| Unmodified develop, same crossing input | -13.12053220665643 | 31.992600 | 50.138461 |
| Fixed noncrossing translation control | -30.89742338437838 | 1.490274 | 1.942380 |

Develop exited normally but the harness failed, with a 17.7768911777 eV
energy error. In the independent noncrossing control, both fractional x
coordinates were translated by `20/45` modulo one (an integer shift on
the 45-point FFT grid). Its energy agrees with the fixed crossing run to
less than 1e-9 eV; force and stress sums agree at the reported precision.

The fresh-fold hybrid unit test verifies the boundary transformation only.
A separate Gamma-point, nonzero-file-A diagnostic still showed an
origin-translation discrepancy with fully frozen atoms, equally before and
after this patch. It is not used as a passing reference for this fix.

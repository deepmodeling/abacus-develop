# 17_DS_DFTU — DeltaSpin & DFT+U Integration Test Suite

This directory contains integration test cases for **DeltaSpin (spin-constrained DFT)** and **DFT+U** functionality in ABACUS,
covering LCAO and PW basis sets, collinear/noncollinear spin, DFT+U, DeltaSpin, and their combinations.

## Test List (47 cases)

### I. LCAO Spin (01-02)

| # | Test Case | Description |
|---|------|------|
| 01 | LCAO_SPIN_S2_Z | Verify basic SCF convergence of collinear spin with LCAO basis, serves as baseline for LCAO magnetic calculations |
| 02 | LCAO_SPIN_S4_XYZ | Verify basic SCF convergence of noncollinear spin with LCAO basis, covers LCAO noncollinear calculation path |

### II. LCAO DFT+U (03-05)

| # | Test Case | Description |
|---|------|------|
| 03 | LCAO_DFTU_S2_Z | Verify coupling of DFT+U (U=5.0eV, l=2) with collinear spin in LCAO basis, ensures correct DFT+U occupation matrix calculation in LCAO path |
| 04 | LCAO_DFTU_S4_XY | Verify coupling of DFT+U with noncollinear spin (XY magnetization) in LCAO basis, covers nspin=4 occupation matrix calculation in LCAO path |
| 05 | LCAO_DFTU_S4_XYZ | Verify coupling of DFT+U with noncollinear spin (XYZ magnetization) in LCAO basis, covers the most complete occupation matrix scenario in LCAO path |

### III. PW Spin (06-07)

| # | Test Case | Description |
|---|------|------|
| 06 | PW_SPIN_S2_Z | Verify basic SCF convergence of collinear spin with PW basis, serves as baseline for PW magnetic calculations |
| 07 | PW_SPIN_S4_XYZ | Verify basic SCF convergence of noncollinear spin with PW basis, covers PW noncollinear calculation path |

### IV. PW DFT+U (08-09, 11)

| # | Test Case | Description |
|---|------|------|
| 08 | PW_DFTU_S2_Z | Verify coupling of DFT+U (U=5.0eV, l=2) with collinear spin in PW basis, ensures correct DFT+U effective potential calculation in PW path |
| 09 | PW_DFTU_S4_XY | Verify coupling of DFT+U with noncollinear spin (XY magnetization) in PW basis, covers onsite projection matrix for nspin=4 in PW path |
| 11 | PW_DFTU_S2_FeO | Verify correctness of DFT+U on FeO system with PW basis, ensures DFT+U correction for Fe-3d orbitals is effective |

### V. PW DeltaSpin (12, 14-16)

| # | Test Case | Description |
|---|------|------|
| 12 | PW_DS_S2_Z | Verify coupling of DeltaSpin with collinear spin in PW basis, ensures correct DeltaSpin iterative optimization of magnetization to target values |
| 14 | PW_DS_S4_XYZ | Verify iterative optimization of noncollinear DeltaSpin under XYZ three-direction magnetization constraint, covers the most complete spin constraint scenario |
| 15 | PW_DS_S4_Z | Verify behavior of noncollinear DeltaSpin when constraining only Z-direction magnetization, ensures uniaxial constraint does not introduce unphysical XY components in noncolin=1 framework |
| 16 | PW_DS_S4_XY | Verify iterative optimization of noncollinear DeltaSpin under XY magnetization constraint with a different crystal structure, verifies generalization of noncollinear DeltaSpin XY constraint under different lattices |

### VI. PW DFT+U + DeltaSpin (18-19, 21)

| # | Test Case | Description |
|---|------|------|
| 18 | PW_DFTU_DS_S2_Z | Verify coupling of DFT+U with DeltaSpin combined (collinear spin) in PW basis, ensures U correction and magnetization constraint do not conflict |
| 19 | PW_DFTU_DS_S4_XY | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XY magnetization constraint, covers joint iteration of both methods in nspin=4 path |
| 21 | PW_DFTU_DS_S4_Z | Verify behavior of noncollinear DFT+U+DeltaSpin combined when constraining only Z-direction magnetization, ensures correct superposition of uniaxial constraint with DFT+U effective potential |

### VII. LCAO DeltaSpin (24, 26-28)

| # | Test Case | Description |
|---|------|------|
| 24 | LCAO_DS_S2_Z | Verify coupling of DeltaSpin with collinear spin in LCAO basis, ensures correct spin constraint optimization in LCAO density matrix path |
| 26 | LCAO_DS_S4_XYZ | Verify iterative optimization of noncollinear DeltaSpin under XYZ three-direction magnetization constraint in LCAO basis, covers the most complete constraint scenario in LCAO path |
| 27 | LCAO_DS_S4_Z | Verify behavior of noncollinear DeltaSpin when constraining only Z-direction magnetization in LCAO basis, ensures correctness of uniaxial constraint in noncolin=1 framework |
| 28 | LCAO_DS_S4_XY | Verify iterative optimization of noncollinear DeltaSpin under XY magnetization constraint in LCAO basis with a different crystal structure, verifies generalization of LCAO noncollinear DeltaSpin XY constraint under different lattices |

### VIII. LCAO DFT+U + DeltaSpin (30-33)

| # | Test Case | Description |
|---|------|------|
| 30 | LCAO_DFTU_DS_S2_Z | Verify coupling of DFT+U with DeltaSpin combined (collinear spin) in LCAO basis, ensures U correction and magnetization constraint do not conflict in density matrix path |
| 31 | LCAO_DFTU_DS_S4_XY | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XY magnetization constraint in LCAO basis, covers joint constraint in LCAO density matrix path |
| 32 | LCAO_DFTU_DS_S4_XYZ | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XYZ three-direction magnetization constraint in LCAO basis, covers the most complete joint scenario in LCAO path |
| 33 | LCAO_DFTU_DS_S4_Z | Verify behavior of noncollinear DFT+U+DeltaSpin combined when constraining only Z-direction magnetization in LCAO basis, ensures correct superposition of uniaxial constraint with DFT+U density matrix |

### IX. PW DeltaSpin Special Parameters (36-41)

| # | Test Case | Description |
|---|------|------|
| 36 | PW_DS_S2_ReadLam_Z | Verify correctness of `nsc=1` mode (read lambda file directly without iterative optimization), ensures DeltaSpin correctly computes magnetization in non-self-consistent lambda mode |
| 37 | PW_DS_S4_ReadLam_XY | Verify `nsc=1` mode for noncollinear DeltaSpin, covers non-self-consistent lambda path under XY magnetization constraint |
| 38 | PW_DS_S2_Thr1e10_Z | Verify stability of DeltaSpin under strict convergence threshold (sc_scf_thr=1e-10), ensures iterative optimization converges to high-precision solution |
| 39 | PW_DS_S4_Thr1e10_XY | Verify stability of noncollinear DeltaSpin under strict convergence threshold (sc_scf_thr=1e-10), covers XY magnetization constraint scenario |
| 40 | PW_DS_S2_Thr10_Z | Verify behavior of DeltaSpin under loose convergence threshold (sc_scf_thr=10), tests algorithm robustness and out_alllog log output under low precision requirements |
| 41 | PW_DS_S4_Thr10_XY | Verify behavior of noncollinear DeltaSpin under loose convergence threshold (sc_scf_thr=10), covers low precision scenario with XY magnetization constraint |

### X. PW DFT+U + DeltaSpin Special Parameters (42-45)

| # | Test Case | Description |
|---|------|------|
| 42 | PW_DFTU_DS_S2_Thr1e10_Z | Verify iterative stability of DFT+U with DeltaSpin combined under strict convergence threshold (sc_scf_thr=1e-10), ensures convergence when both methods are coupled |
| 43 | PW_DFTU_DS_S4_Thr1e10_XY | Verify coupling stability of noncollinear DFT+U+DeltaSpin under strict convergence threshold (sc_scf_thr=1e-10), covers XY magnetization constraint |
| 44 | PW_DFTU_DS_S2_Thr10_Z | Verify behavior of DFT+U with DeltaSpin combined under loose convergence threshold (sc_scf_thr=10), tests coupled algorithm robustness under low precision requirements |
| 45 | PW_DFTU_DS_S4_Thr10_XY | Verify behavior of noncollinear DFT+U+DeltaSpin under loose convergence threshold (sc_scf_thr=10), covers low precision scenario with XY magnetization constraint |

### XI. FeO Atom Ordering (50-51)

| # | Test Case | Description |
|---|------|------|
| 50 | FeO_O_first_Fe_second | Verify correctness of DFT+U in FeO system with O atom type first and Fe second, ensures atom type ordering does not affect DFT+U onsite projection |
| 51 | FeO_Fe_first_O_second | Verify correctness of DFT+U in FeO system with Fe atom type first and O second, compare with 50 to ensure eff_pot_pw_index indexing is independent of atom type ordering |

### XII. NSCF Mode (55, 60-64)

**Note:** These test cases have been converted to use a **SCF+NSCF workflow**.
The pre-converged charge density files have been removed. To run these tests:

```bash
# Use the workflow script (runs SCF first, then NSCF)
cd tests/17_DS_DFTU/55_PW_DS_NSCF_S4_XY
bash ../run_scf_nscf.sh <abacus_path> 4
```

| # | Test Case | Description |
|---|------|------|
| 55 | PW_DS_NSCF_S4_XY | Verify DeltaSpin functionality in non-self-consistent (nscf) calculation mode, ensures lambda constraint is applied correctly without charge update |
| 60 | PW_DFTU_DS_NSCF_Band_XY | Verify DFT+U+DeltaSpin in NSCF band structure calculation, tests band output with spin constraints on high-symmetry k-point path |
| 61 | LCAO_DS_NSCF_S4_XY | Verify LCAO DeltaSpin functionality in nscf calculation mode |
| 62 | LCAO_DFTU_NSCF_Band_XY | Verify LCAO DFT+U (without DeltaSpin) in NSCF band structure calculation; note: runs as `calculation = scf` with `scf_nmax = 1` using pre-converged charge density |
| 63 | LCAO_DFTU_DS_NSCF_Band_XY | Verify LCAO DFT+U+DeltaSpin in NSCF band structure calculation, tests band output with spin constraints |
| 64 | PW_DFTU_NSCF_Band_XY | Verify DFT+U (without DeltaSpin) in NSCF band structure calculation, tests band output with Hubbard U correction |

**SCF+NSCF Workflow:**
1. `scf/INPUT` — SCF input (calculation=scf, init_chg=atomic, out_chg=1)
2. `scf/STRU`, `scf/KPT` — SCF structure and k-points
3. `run_scf_nscf.sh` — Script that runs SCF, copies charge density, then runs NSCF
4. CI tests are **disabled** for these cases (see CASES_CPU.txt)

### XIII. sc_direction_only Constraint (56-59)

| # | Test Case | Description |
|---|------|------|
| 56 | PW_DS_S4_DirectionOnly_XY | Verify `sc_direction_only=1` mode: only magnetization direction is constrained while magnitude is free to relax, projects lambda perpendicular to target direction |
| 57 | PW_DFTU_DS_S4_DirectionOnly_XY | Verify `sc_direction_only=1` combined with DFT+U, tests direction-only constraint superposition with Hubbard U correction |
| 58 | LCAO_DS_S4_DirectionOnly_XY | Verify `sc_direction_only=1` in LCAO basis, ensures direction-only constraint works correctly in LCAO density matrix path |
| 59 | LCAO_DFTU_DS_S4_DirectionOnly_XY | Verify `sc_direction_only=1` combined with DFT+U in LCAO basis, tests full direction-only constraint in LCAO path |

## Running Tests

```bash
# Run all tests
cd tests/17_DS_DFTU
bash ../integrate/Autotest.sh -a <abacus_path> -n 4

# Run a single test
cd 08_PW_DFTU_S2_Z
bash ../../integrate/run_debug.sh ""
```

## CI-Disabled Tests

The following test cases are disabled in `CASES_CPU.txt` (commented out with `#`) and excluded from CI testing due to **convergence and numerical stability issues**. They can be manually unskipped for local testing by removing the `#` prefix.

| # | Test Case | Reason |
|---|------|--------|
| 02 | LCAO_SPIN_S4_XYZ | Convergence / numerical stability |
| 04 | LCAO_DFTU_S4_XY | Convergence / numerical stability |
| 05 | LCAO_DFTU_S4_XYZ | Convergence / numerical stability |
| 24 | LCAO_DS_S2_Z | Convergence / numerical stability |
| 26 | LCAO_DS_S4_XYZ | Convergence / numerical stability |
| 27 | LCAO_DS_S4_Z | Convergence / numerical stability |
| 28 | LCAO_DS_S4_XY | Convergence / numerical stability |
| 30 | LCAO_DFTU_DS_S2_Z | Convergence / numerical stability |
| 31 | LCAO_DFTU_DS_S4_XY | Convergence / numerical stability |
| 32 | LCAO_DFTU_DS_S4_XYZ | Convergence / numerical stability |
| 33 | LCAO_DFTU_DS_S4_Z | Convergence / numerical stability |
| 44 | PW_DFTU_DS_S2_Thr10_Z | Convergence / numerical stability |
| 58 | LCAO_DS_S4_DirectionOnly_XY | Convergence / numerical stability |
| 59 | LCAO_DFTU_DS_S4_DirectionOnly_XY | Convergence / numerical stability |
| 62 | LCAO_DFTU_NSCF_Band_XY | Convergence / numerical stability; genelpa eigenvalue inconsistency across thread counts (scalapack_gvx consistent); **also disabled for SCF+NSCF workflow conversion** |
| 63 | LCAO_DFTU_DS_NSCF_Band_XY | Convergence / numerical stability; **also disabled for SCF+NSCF workflow conversion** |
| 55 | PW_DS_NSCF_S4_XY | **Disabled for SCF+NSCF workflow conversion** — run manually with `run_scf_nscf.sh` |
| 60 | PW_DFTU_DS_NSCF_Band_XY | **Disabled for SCF+NSCF workflow conversion** — run manually with `run_scf_nscf.sh` |
| 61 | LCAO_DS_NSCF_S4_XY | **Disabled for SCF+NSCF workflow conversion** — run manually with `run_scf_nscf.sh` |
| 64 | PW_DFTU_NSCF_Band_XY | **Disabled for SCF+NSCF workflow conversion** — run manually with `run_scf_nscf.sh` |

## Test Condition Notes

- 09 (PW DFT+U + noncollinear): Only supports **2-process MPI** execution, `result.ref` reference files provided
- The following test cases set `kpar=2` in INPUT and require at least **2 MPI processes** to run: 11, 12, 14, 15, 16, 18, 19, 21, 37, 39, 41, 43, 45
- 62 (LCAO_DFTU_NSCF_Band_XY): Single-thread and multi-thread results are inconsistent; investigation shows HR, HK, and SK are consistent across threads, but eigenvalues from genelpa differ; switching to scalapack_gvx produces consistent results across thread counts. Note: this test is named "NSCF" but actually runs with `calculation = scf` (`scf_nmax = 1`), using pre-shipped charge density and onsite.dm files as initial guess
- All NSCF tests (55, 60, 61, 62, 63, 64) have been **converted to SCF+NSCF workflow**:
  - Pre-converged `autotest-CHARGE-DENSITY.restart` and `onsite.dm` files have been removed
  - Each test directory contains a `scf/` subdirectory with SCF input files
  - Run with: `bash ../run_scf_nscf.sh <abacus_path> [mpi_np]`
  - These tests are **disabled in CI** (commented out in CASES_CPU.txt)
- All LCAO basis tests use `ks_solver = genelpa`. The genelpa eigenvalue inconsistency across thread counts observed in test 62 may potentially affect other LCAO tests as well
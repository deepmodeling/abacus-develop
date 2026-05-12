# 17_DS_DFTU — DeltaSpin & DFT+U Integration Test Suite

This directory contains all integration test cases for **DeltaSpin (spin-constrained DFT)** and **DFT+U** functionality in ABACUS,
covering LCAO and PW basis sets, collinear/noncollinear spin, DFT+U, DeltaSpin, and their combinations.

## Test List (60 cases)

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

### IV. PW DFT+U (08-11)

| # | Test Case | Description |
|---|------|------|
| 08 | PW_DFTU_S2_Z | Verify coupling of DFT+U (U=5.0eV, l=2) with collinear spin in PW basis, ensures correct DFT+U effective potential calculation in PW path |
| 09 | PW_DFTU_S4_XY | Verify coupling of DFT+U with noncollinear spin (XY magnetization) in PW basis, covers onsite projection matrix for nspin=4 in PW path |
| 10 | PW_DFTU_S4_XY | Same parameters as 09 but different crystal structure, verifies generalization of PW DFT+U noncollinear under different lattices |
| 11 | PW_DFTU_S2_FeO | Verify correctness of DFT+U on FeO system with PW basis, ensures DFT+U correction for Fe-3d orbitals is effective |

### V. PW DeltaSpin (12-17)

| # | Test Case | Description |
|---|------|------|
| 12 | PW_DS_S2_Z | Verify coupling of DeltaSpin with collinear spin in PW basis, ensures correct DeltaSpin iterative optimization of magnetization to target values |
| 13 | PW_DS_S4_XY | Verify iterative optimization of noncollinear DeltaSpin under XY magnetization constraint, covers lambda update for nspin=4 path |
| 14 | PW_DS_S4_XYZ | Verify iterative optimization of noncollinear DeltaSpin under XYZ three-direction magnetization constraint, covers the most complete spin constraint scenario |
| 15 | PW_DS_S4_Z | Verify behavior of noncollinear DeltaSpin when constraining only Z-direction magnetization, ensures uniaxial constraint does not introduce unphysical XY components in noncolin=1 framework |
| 16 | PW_DS_S4_XY | Same parameters as 13 but different crystal structure, verifies generalization of noncollinear DeltaSpin XY constraint under different lattices |
| 17 | PW_DS_S4_XYZ | Same parameters as 14 but different crystal structure, verifies generalization of noncollinear DeltaSpin XYZ constraint under different lattices |

### VI. PW DFT+U + DeltaSpin (18-23)

| # | Test Case | Description |
|---|------|------|
| 18 | PW_DFTU_DS_S2_Z | Verify coupling of DFT+U with DeltaSpin combined (collinear spin) in PW basis, ensures U correction and magnetization constraint do not conflict |
| 19 | PW_DFTU_DS_S4_XY | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XY magnetization constraint, covers joint iteration of both methods in nspin=4 path |
| 20 | PW_DFTU_DS_S4_XYZ | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XYZ three-direction magnetization constraint, covers the most complete joint constraint scenario |
| 21 | PW_DFTU_DS_S4_Z | Verify behavior of noncollinear DFT+U+DeltaSpin combined when constraining only Z-direction magnetization, ensures correct superposition of uniaxial constraint with DFT+U effective potential |
| 22 | PW_DFTU_DS_S4_XY | Same parameters as 19 but different crystal structure, verifies generalization of noncollinear DFT+U+DeltaSpin combined under different lattices |
| 23 | PW_DFTU_DS_S4_XYZ | Same parameters as 20 but different crystal structure, verifies generalization of noncollinear DFT+U+DeltaSpin combined XYZ constraint under different lattices |

### VII. LCAO DeltaSpin (24-29)

| # | Test Case | Description |
|---|------|------|
| 24 | LCAO_DS_S2_Z | Verify coupling of DeltaSpin with collinear spin in LCAO basis, ensures correct spin constraint optimization in LCAO density matrix path |
| 25 | LCAO_DS_S4_XY | Verify iterative optimization of noncollinear DeltaSpin under XY magnetization constraint in LCAO basis, covers magnetization projection for nspin=4 in LCAO path |
| 26 | LCAO_DS_S4_XYZ | Verify iterative optimization of noncollinear DeltaSpin under XYZ three-direction magnetization constraint in LCAO basis, covers the most complete constraint scenario in LCAO path |
| 27 | LCAO_DS_S4_Z | Verify behavior of noncollinear DeltaSpin when constraining only Z-direction magnetization in LCAO basis, ensures correctness of uniaxial constraint in noncolin=1 framework |
| 28 | LCAO_DS_S4_XY | Same parameters as 25 but different crystal structure, verifies generalization of LCAO noncollinear DeltaSpin XY constraint under different lattices |
| 29 | LCAO_DS_S4_XYZ | Same parameters as 26 but different crystal structure, verifies generalization of LCAO noncollinear DeltaSpin XYZ constraint under different lattices |

### VIII. LCAO DFT+U + DeltaSpin (30-35)

| # | Test Case | Description |
|---|------|------|
| 30 | LCAO_DFTU_DS_S2_Z | Verify coupling of DFT+U with DeltaSpin combined (collinear spin) in LCAO basis, ensures U correction and magnetization constraint do not conflict in density matrix path |
| 31 | LCAO_DFTU_DS_S4_XY | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XY magnetization constraint in LCAO basis, covers joint constraint in LCAO density matrix path |
| 32 | LCAO_DFTU_DS_S4_XYZ | Verify coupling of noncollinear DFT+U+DeltaSpin combined under XYZ three-direction magnetization constraint in LCAO basis, covers the most complete joint scenario in LCAO path |
| 33 | LCAO_DFTU_DS_S4_Z | Verify behavior of noncollinear DFT+U+DeltaSpin combined when constraining only Z-direction magnetization in LCAO basis, ensures correct superposition of uniaxial constraint with DFT+U density matrix |
| 34 | LCAO_DFTU_DS_S4_XY | Same parameters as 31 but different crystal structure, verifies generalization of LCAO DFT+U+DeltaSpin combined under different lattices |
| 35 | LCAO_DFTU_DS_S4_XYZ | Same parameters as 32 but different crystal structure, verifies generalization of LCAO DFT+U+DeltaSpin combined XYZ constraint under different lattices |

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

### XI. BFGS Lambda Strategy (46-49)

| # | Test Case | Description |
|---|------|------|
| 46 | PW_DS_S2_Thr1e10_Z_bfgs | Verify convergence behavior of DeltaSpin using BFGS strategy (sc_lambda_strategy=bfgs), tests BFGS optimizer correctness in spin-constrained SCF |
| 47 | PW_DS_S4_Thr1e10_XY_bfgs | Verify convergence behavior of noncollinear DeltaSpin using BFGS strategy, covers correctness of BFGS optimizer under XY magnetization constraint |
| 48 | PW_DFTU_DS_S2_Thr10_Z_bfgs | Verify convergence behavior of DFT+U with DeltaSpin combined using BFGS strategy, tests BFGS correctness in DFT+U+DS coupled scenario |
| 49 | PW_DFTU_DS_S4_Thr10_XY_bfgs | Verify convergence behavior of noncollinear DFT+U+DeltaSpin combined using BFGS strategy, covers correctness of BFGS optimizer under XY magnetization constraint |

### XII. FeO Atom Ordering (50-51)

| # | Test Case | Description |
|---|------|------|
| 50 | FeO_O_first_Fe_second | Verify correctness of DFT+U in FeO system with O atom type first and Fe second, ensures atom type ordering does not affect DFT+U onsite projection |
| 51 | FeO_Fe_first_O_second | Verify correctness of DFT+U in FeO system with Fe atom type first and O second, compare with 50 to ensure eff_pot_pw_index indexing is independent of atom type ordering |

### XIII. SOC + DFT+U (52)

| # | Test Case | Description |
|---|------|------|
| 52 | PW_DFTU_SO | Verify compatibility when DFT+U and spin-orbit coupling (SOC) are enabled simultaneously, ensures DFT+U onsite projection correctly couples with SOC spin mixing |

### XIV. Magnetic Moment & Lambda Verification (53-54)

| # | Test Case | Description |
|---|------|------|
| 53 | PW_DS_S4_XY_MagMomCheck | Verify atomic magnetic moments converge to target values and lambda values are within expected range for PW DeltaSpin |
| 54 | PW_DFTU_DS_S4_XY_MagMomCheck | Verify atomic magnetic moments and lambda values for PW DFT+U+DeltaSpin combined, ensures both corrections coexist correctly |

### XV. NSCF Mode (55, 60)

| # | Test Case | Description |
|---|------|------|
| 55 | PW_DS_NSCF_S4_XY | Verify DeltaSpin functionality in non-self-consistent (nscf) calculation mode, ensures lambda constraint is applied correctly without charge update |
| 60 | PW_DFTU_NSCF_Band_XY | Verify DFT+U+DeltaSpin in NSCF band structure calculation, tests band output with spin constraints on high-symmetry k-point path |

### XVI. sc_direction_only Constraint (56-59)

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

## Known Issues

- 19-23: PW DFT+U + DeltaSpin + noncollinear → both port and zdy-tmp crash (upstream bug)

## Test Condition Notes

- 09/10 (PW DFT+U + noncollinear): Only supports **2-process MPI** execution, `result.ref` reference files provided

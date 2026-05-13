# DeltaSpin ~1000 eV Energy Deviation Debug Log

## Basic Information

| Item | Value |
|---|---|
| **Debug Date** | 2026-05-08 |
| **Current Branch** | `feat/dftu-pw-port-v2` |
| **Current Version** | ABACUS v3.11.0-beta.1 |
| **Baseline Branch** | tmp (zdy) |
| **Baseline Version** | ABACUS v3.7.1 |
| **Current Build** | `/root/abacus-dftu-pw-port/build/abacus_pw_para` |
| **Baseline Build** | `/root/abacus-zdy-tmp/build_zdy/abacus_pw` |
| **Remote Repository** | `dyzheng/abacus-develop.git` (zdy) |

---

## Test Case Under Investigation

| Item | Value |
|---|---|
| **Test Path** | `/root/abacus-dftu-pw-port/tests/17_DS_DFTU/13_PW_DS_S4_XY/` |
| **Test Name** | PW_DS_S4_XY |
| **System** | 2 Fe atoms, PW basis, nspin=4, kpar=4, 4-process MPI |
| **Constraint Type** | XY-plane spin constraint (target: atom1=[2,0,0], atom2=[-2,0,0]) |
| **tmp Branch Reference Output** | `/root/abacus-dftu-pw-port/tests/17_DS_DFTU/13_PW_DS_S4_XY/log2-1` |
| **tmp Branch Test Output** | `/root/test_tmp_cmp/13_PW_DS_S4_XY/OUT.autotest/running_scf.log` |

---

## Related Source Files

### DeltaSpin Module (compared, identical)

| File Path | Status |
|---|---|
| `source/source_lcao/module_deltaspin/lambda_loop.cpp` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/lambda_loop_helper.cpp` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/cal_mw_from_lambda.cpp` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/cal_mw.cpp` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/basic_funcs.cpp` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/spin_constrain.h` | ✓ diff with tmp branch is empty |
| `source/source_lcao/module_deltaspin/update_psi_charge.cpp` | ✓ diff with tmp branch is empty |

### PW Layer Related Code (compared, identical)

| File Path | Status |
|---|---|
| `source/source_pw/module_pwdft/deltaspin_pw.cpp` | ✓ diff with tmp branch is empty |
| `source/source_pw/module_pwdft/deltaspin_pw.h` | ✓ diff with tmp branch is empty |
| `source/source_pw/module_pwdft/op_pw_proj.cpp` | ✓ diff with tmp branch is empty |
| `source/source_pw/module_pwdft/onsite_proj.cpp` | ✓ diff with tmp branch is empty |
| `source/source_pw/module_pwdft/dftu_pw.cpp` | ✓ diff with tmp branch is empty |

### Parallel-Related Code (fixed)

| File Path | Changes |
|---|---|
| `source/source_lcao/module_dftu/dftu_pw.cpp` | `PARAM.globalv.nproc_in_pool` → `GlobalV::NPROC_IN_POOL` (3 places) |
| `source/source_lcao/module_deltaspin/cal_mw.cpp` | `PARAM.globalv.nproc_in_pool` → `GlobalV::NPROC_IN_POOL` (1 place) |
| `source/source_lcao/module_deltaspin/cal_mw_from_lambda.cpp` | `PARAM.globalv.nproc_in_pool` → `GlobalV::NPROC_IN_POOL` (1 place) |
| `source/source_lcao/module_deltaspin/lambda_loop_helper.cpp` | Fixed `print` undeclared compilation error |

### Diagnostic Code (added)

| File Path | Changes |
|---|---|
| `source/source_lcao/module_deltaspin/lambda_loop.cpp` | Added `ds_diag` namespace with expected values from tmp branch log2-1, outputs `[DS-DIAG]` diagnostics at 4 key checkpoints |

---

## Debug Process

### 1. nproc_in_pool Fix (completed ✓)

**Problem**: `PARAM.globalv.nproc_in_pool` defaults to 1, `divide_pools()` only updates `GlobalV::NPROC_IN_POOL`, causing `reduce_double_allpool` to double-count locale values when kpar>1.

**Fix**: Replaced 5 occurrences of `PARAM.globalv.nproc_in_pool` with `GlobalV::NPROC_IN_POOL`.

**Verification**: test 09 (PW_DFTU_S4_XY, pure DFT+U without DeltaSpin) results match:
- tmp: ETOT = -6364.254 eV
- current: ETOT = -6364.266 eV
- Difference: ~0.01 eV, acceptable

### 2. DeltaSpin ~1000 eV Deviation Debug (in progress ✗)

**Symptom**: All S4+DS test cases (test 13-23, 41, 43, 45) have ETOT deviating from tmp branch by ~1000 eV.

### Diagnostic Data (test 13 single-process run)

#### SCF Iteration Comparison

| Iteration | Current TMAGX | tmp TMAGX | Current ETOT (eV) | tmp ETOT (eV) | Difference (eV) |
|---|---|---|---|---|---|
| DS1 | -9.91e-07 | -1.59e-06 | -6358.66605 | -6358.64595 | -0.020 |
| DS2 | 1.49e-04 | -1.37e-03 | -6362.79540 | -6362.80760 | +0.012 |
| DS3 | -5.00e-04 | -2.17e-03 | -6368.63682 | -6368.53835 | -0.098 |
| DS4 | 4.07e-03 | -8.84e-03 | -6367.85915 | -6367.60040 | -0.259 |
| DS5 | -2.24e-01 | -5.00e-01 | -6370.49481 | -6370.42715 | -0.068 |
| DS6 | 8.75e-02 | -4.02e-02 | -6370.57871 | -6370.54764 | -0.031 |
| DS7 | -1.42e-01 | 1.25e-01 | -6370.60386 | -6370.60460 | +0.001 |
| DS8 | 3.71e-02 | 1.12e-02 | -6370.62546 | -6370.61781 | -0.008 |
| DS9* | -4.90e-02 | 2.67e-02 | -6372.32067 | -6370.62092 | **-1.700** |
| DS10 | 2.13e-02 | -5.17e-03 | -5348.77865 | -6361.08653 | **~1000** |
| DS11 | 6.35e-02 | -8.95e-03 | -5315.65026 | -6369.02313 | **~1000** |

> DS9: Iteration after SCF restart, also the entry point to the lambda optimization loop for the first time

#### Lambda Loop Diagnostic Output `[DS-DIAG]`

| Checkpoint | Current Value | tmp Expected Value | Relative Error | Status |
|---|---|---|---|---|
| **Outer=9, Inner=-1 (initial spin)** | | | | |
| spin[0].x | -3.13049 | -3.26060 | 3.99% | **FAIL** |
| spin[0].y | 0.000221 | -0.005538 | 104.0% | **FAIL** |
| spin[0].z | -0.002899 | -0.006936 | 58.2% | **FAIL** |
| spin[1].x | 3.36201 | 3.24301 | 3.67% | **FAIL** |
| spin[1].y | -0.001333 | 0.006276 | 121.2% | **FAIL** |
| spin[1].z | 0.001583 | 0.004310 | 63.3% | **FAIL** |
| lambda (all) | 0 | 0 | 0% | PASS |
| **Outer=9, Inner=0 (first step)** | | | | |
| RMS | 5.24753 | 5.25182 | 0.08% | PASS |
| alpha_opt | 0.337484 | 0.570277 | 40.8% | **FAIL** |
| **Outer=9, Inner=1** | | | | |
| RMS | 1.01044 | 0.818055 | 23.5% | **FAIL** |
| **Outer=9, Inner=2** | | | | |
| RMS | 0.316889 | 0.235267 | 34.7% | **FAIL** |
| **Outer=9, Inner=3** | | | | |
| RMS | 0.0129247 | 0.0274929 | 53.0% | **FAIL** |
| **Outer=9, Inner=4 (converged)** | | | | |
| spin[0].x | 2.00026 | 2.00037 | 0.0055% | PASS |
| spin[0].y | -8.69e-06 | 3.13e-04 | 102.8% | **FAIL** |
| spin[0].z | -1.47e-04 | 2.75e-04 | 153.3% | **FAIL** |
| spin[1].x | -2.00069 | -1.99996 | 0.037% | PASS |
| lambda[0].x (eV) | -3.5935 | -3.70354 | 2.97% | **FAIL** |
| lambda[0].y (eV) | 1.89e-04 | -4.69e-03 | 104.0% | **FAIL** |
| lambda[0].z (eV) | -2.38e-03 | -5.82e-03 | 59.0% | **FAIL** |
| RMS | 0.0005335 | 0.0005114 | 4.33% | FAIL |
| **Outer=9, Inner=999 (Current RMS after SCF)** | | | | |
| Current RMS | 0.003254 | 0.124618 | 97.4% | **FAIL** |

---

## Analysis Conclusions

### Ruled Out

1. ✅ DeltaSpin core code (lambda_loop, cal_mw, basic_funcs, etc.) is **identical** to tmp branch
2. ✅ DFT+U code (dftu_pw.cpp) fixed and verified with test 09
3. ✅ Parallel reduction logic (reduce_double_allpool) fixed

### Key Findings

1. **SCF phase DS1-DS8**: ETOT is mostly consistent (difference < 0.1 eV), but TMAGX/Y/Z shows small deviations from DS1, accumulating to 1.7 eV by DS9
2. **Lambda loop initial conditions differ**: Spin values after SCF restart differ by 3-4%, causing the entire lambda optimization path to diverge
3. **Lambda loop itself converges**: RMS drops to ~0.0005 in both versions, but converged lambda/spin values differ (y/z components differ significantly)
4. **DS10 energy jump ~1000 eV**: After lambda loop returns, `update_psi_charge` updates wavefunctions and charge density, resulting in completely wrong energy
5. **Current RMS after SCF**: current=0.0033 vs tmp=0.125, differs by 97.4%

### Investigation Directions

| Priority | Direction | Description |
|---|---|---|
| P0 | `update_psi_charge_pw_cpu` and underlying functions | The base code that `update_psi_charge` depends on (Hamilt, ElecState, Charge, etc.) may have changed between v3.7.1→v3.11.0 |
| P0 | SCF restart logic and charge mixing parameters | DS9 energy differs by 1.7 eV after SCF restart; charge mixing implementation or parameters may differ |
| P1 | HSolver (diagonalization solver) | PW diagonalization implementation differences may yield different wavefunctions for the same Hamiltonian |
| P1 | Charge/density reduction logic | Whether `psiToRho`, `reduce_diff_pools`, etc. behave consistently under multi-process |
| P2 | DFT+U and DeltaSpin interaction | Whether the superposition of DFT+U and DeltaSpin in `OnsiteProj::cal_ps_delta_spin` is correct |

---

## Git Commit History

```
ab1863479 Remove deprecated sc_scf_nmin parameter
61f044c8e Update nspin=4 test references after nproc_in_pool fix
15459c692 Update test references for kpar>1 fix
4d567a4cd Fix DeltaSpin reduce_double_allpool using GlobalV::NPROC_IN_POOL
73163aa3e Fix DFT+U locale double-counting with kpar>1
```

## Diagnostic Code Notes

Diagnostic code added to the `ds_diag` namespace in `lambda_loop.cpp`, including:
- Expected values (spin, lambda, RMS, alpha) for outer steps 9-12 extracted from tmp branch log2-1
- 4 checkpoints calling `ds_diag::check_step()`:
  1. Initial spin/lambda (inner=-1)
  2. After each step RMS calculation
  3. After alpha_opt calculation
  4. Current RMS after SCF convergence

Search for `[DS-DIAG]` tag in runtime output to quickly locate diagnostic information.

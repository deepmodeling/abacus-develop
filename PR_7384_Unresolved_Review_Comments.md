# PR #7384 Unresolved Review Comments

> **PR:** [Feature: DeltaSpin for LCAO and PW base and DFTU for PW, both collinear and noncollinear spin](https://github.com/deepmodeling/abacus-develop/pull/7384)
> **Reviewer:** @dzzz2001 (May 28, 2026)
> **Total Comments:** 30 (15 Blockers, 15 Should-fix)
> **All comments are currently unresolved (no replies from author)**

---

## Blockers (15)

### 1. CUDA onsite_op.cu — DFT+U kernel missing npol==1 branch
- **File:** `source/source_pw/module_pwdft/kernels/cuda/onsite_op.cu` (line 46)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629318
- **Issue:** The CPU `onsite_ps_op` DFT+U overload gained an `if(npol == 2) {…} else {…}` branch, but the second `__global__ void onsite_op(…)` variant (DFT+U kernel taking `vu_iat`) was **not** patched. For npol==1 it reads past `becp`, writes `ps[psind+1]` into the next atom's row, and indexes off-diagonal `vu` entries that don't exist. Same issue in the ROCm mirror (`kernels/rocm/onsite_op.hip.cu`).

### 2. CUDA force_op.cu — Device kernels still hardcode spin factor 2
- **File:** `source/source_pw/module_pwdft/kernels/cuda/force_op.cu` (line 457)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629333
- **Issue:** Only the function signature gained `const int& npol`; the device kernels (`cal_force_npw_op` / `cal_force_nl_op`) still hardcode `ib * 2`, `(ib2+1) * nkb`, and `ipol * nbands * 2 * nkb`. For npol==1 these will read OOB on `dbecp`/`becp`. The ROCm mirror (`kernels/rocm/force_op.hip.cu`) doesn't appear to have the new `npol` parameter in its signature at all, causing a link failure on ROCm builds.

### 3. CUDA stress_op.cu — Device kernels still hardcode spin factor 2
- **File:** `source/source_pw/module_pwdft/kernels/cuda/stress_op.cu` (line 1054)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629337
- **Issue:** Same as force_op.cu: only the signature gained `npol`; the kernels at lines 317, 946, 1008 still hardcode `ib * 2` and `(ib2+1) * nkb`. For npol==1 this will OOB-read `becp`/`dbecp` and double-count. Same for `kernels/rocm/stress_op.hip.cu`.

### 4. cusolver.h — Non-existent cuSOLVER symbols in CUDA < 11 fallback
- **File:** `source/source_base/module_container/base/third_party/cusolver.h` (line 57)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629339
- **Issue:** The `#if CUDA_VERSION < 11000` fallback calls `cusolverDnStrtri_bufferSize` / `cusolverDnStrtri` (and `D`/`C`/`Z` variants). These symbols **do not exist** in any cuSOLVER release. Any build with CUDA < 11 will fail to link.
- **Options:** (1) Drop the fallback and require CUDA ≥ 11; (2) Route to `cublas{S,D,C,Z}trtri_batched` with `batchCount=1`.

### 5. force_op.cpp — npol==1 branch missing spin_sign
- **File:** `source/source_pw/module_pwdft/kernels/force_op.cpp` (line 411)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629342
- **Issue:** The npol==1 branch (`local_force[ipol] -= fac * lambda[iat*3+2] * dbb;`) omits the spin sign that is correctly applied in the Hamiltonian at `op_pw_proj.cpp:1448-1456`. For nspin==2, the force from spin-down k-points (`isk==1`) is added with the wrong sign, so force ≠ `-dE/dR`.

### 6. stress_op.cpp — npol==1 branch missing spin_sign
- **File:** `source/source_pw/module_pwdft/kernels/stress_op.cpp` (line 360)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629352
- **Issue:** Same as force_op.cpp: `coefficients0 = lambda[iat*3+2]` is missing the `spin_sign`. For nspin==2 the stress contribution from spin-down k-points will have the wrong sign; stress tensor inconsistent with `dE/dε`.

### 7. dspin_lcao.h — set_current_spin shadows non-virtual base method
- **File:** `source/source_lcao/module_operator_lcao/dspin_lcao.h` (line 73)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629358
- **Issue:** `DeltaSpin<...>::set_current_spin(int)` shadows (does not override) `OperatorLCAO<TK,TR>::set_current_spin`, which is not `virtual`. The dispatch chain calls the **base** method via `dynamic_cast`, so the derived reset `this->sc_hr_done = false` never fires. The cross-spin accumulation bug is still live.
- **Fix:** Make the base method `virtual` or reset `sc_hr_done` from a different hook (e.g. inside `contributeHR`).

### 8. dspin_lcao.cpp — zgemm called with wrong matrix layout
- **File:** `source/source_lcao/module_operator_lcao/dspin_lcao.cpp` (line 620)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629361
- **Issue:** `D_I` is filled row-major as `D_I[lm * nbands_global + jb_global]`, but passed to column-major BLAS `zgemm` with `lda=r`. The column-major view corresponds to the *transpose* of the stored data. Result: computes `D · D^H` instead of `D^H · D`.
- **Fix:** Fill `D_I` column-major (`D_I[lm + jb_global * r] = ...`) or call zgemm with `lda = nbands_global` and correct transpose flags.

### 9. lambda_loop.cpp — direction_only mode permanently mutates target_mag_
- **File:** `source/source_lcao/module_deltaspin/lambda_loop.cpp` (line 154)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629366
- **Issue:** `this->target_mag_[ia] += parallel * dir;` permanently mutates the user-supplied target magnetization on every inner step and is **never restored**. After the first inner step, `target_mag_` no longer reflects user input, and corruption persists across SCF iterations.
- **Fix:** Keep a fresh copy of the original target, or compute the perpendicular component without modifying `target_mag_`.

### 10. dftu_lcao.cpp — nspin=4 read-from-file branch has wrong access pattern
- **File:** `source/source_lcao/module_operator_lcao/dftu_lcao.cpp` (line 196)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629373
- **Issue:** The access pattern treats `locale[iat][L][0][0]` as a `(2*tlp1) × (2*tlp1)` matrix indexed as a 2×2 block of Pauli channels. But every writer stores the four Pauli channels as **four stacked `tlp1*tlp1` blocks** at offsets `0, tlp1², 2*tlp1², 3*tlp1²`. Restarts with `init_chg=file` + nspin=4 will silently produce a corrupted DFT+U potential.

### 11. dftu.cpp — Early continue skips per-atom resize calls
- **File:** `source/source_lcao/module_dftu/dftu.cpp` (line 94)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629380
- **Issue:** The early `continue` for atoms without correlated orbitals runs **before** `locale[iat].resize(...)`, `iatlnmipol2iwt[iat]...resize(...)`, and `eff_pot_pw_index[iat] = pot_index`. These per-atom structures are left default-constructed (size 0), which will segfault if any future code accesses them without the `orbital_corr[it] == -1` guard.
- **Fix:** Hoist the guard so only Hubbard-specific work is skipped, but resize calls still run.

### 12. lambda_strategy_integration.cpp — Incomplete file will not compile
- **File:** `source/source_lcao/module_deltaspin/lambda_strategy_integration.cpp` (line 6)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629387
- **Issue:** The file's own header comment states "Status: INCOMPLETE … will not compile as-is." It references `strategy_type_`, `strategy_`, `LambdaStrategyType`, `set_strategy_type`, `set_strategy_params` — none declared in `spin_constrain.h`. Companion `lambda_update_strategies.{h,cpp}` (~700 lines) is labelled "NOT compiled into the library."
- **Fix:** Either finish + register with CMake + unit-test, or drop from the PR.

### 13. read_atoms_helper.cpp — Stray diagnostic std::cout left in production code
- **File:** `source/source_cell/read_atoms_helper.cpp` (line 456)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629388
- **Issue:** Stray diagnostic `std::cout << "[DS-DIAG] STRU parse: lambda[…]…"` left in the production STRU parse path. Fires on every STRU read, on every MPI rank, clutters every log.
- **Fix:** Remove or gate behind a verbosity flag.

### 14. charge_mixing.cpp — mix_uom has empty timer and no MPI sync
- **File:** `source/source_estate/module_charge/charge_mixing.cpp` (line 276)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629396
- **Issues:**
  1. **Timer is empty:** `timer::start`/`timer::end` called with no work between them. Actual work is outside the timer span.
  2. **No MPI synchronization:** The DFT+U occupation matrix `uom_in` is mixed locally on each rank with no `MPI_Bcast`/`Allreduce`. After a few Pulay steps, occupation matrices on different ranks will silently diverge.

### 15. Test files — Tests don't exercise production code
- **File:** `source/source_lcao/module_dftu/test/dftu_pw_test.cpp` (line 32)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629400
- **Issue:** None of the new DFT+U/DeltaSpin test files (`deltaspin_core_test.cpp`, `deltaspin_pw_test.cpp`, `dftu_core_test.cpp`, `dftu_operator_test.cpp`, `dftu_pw_test.cpp`, ≈4,200 lines combined) include or call **any** production header from the modules they purport to test. Each file reimplements the algorithm as a static helper, then asserts the reimplementation matches hand-computed numbers. The production code (`Plus_U` / `SpinConstrain` / `OnsiteProjector` APIs) has **zero** direct unit-test coverage.
- **Fix:** Rewrite tests to `#include` the relevant headers and exercise the real APIs.

---

## Should-fix (15)

### 1. cal_mw_from_lambda.cpp — DMR may be stale in LCAO path
- **File:** `source/source_lcao/module_deltaspin/cal_mw_from_lambda.cpp` (line 479)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629413
- **Issue:** The previous code called `cal_dm_psi` + `cal_DMR()` to refresh the DM/DMR after wavefunction update. The new `cal_mi_lcao` reads `get_DMR_pointer(1)` directly. If `HSolverLCAO::solve(..., /*skip_charge=*/true)` does not rebuild DM/DMR internally, the moment is computed from stale DMR.

### 2. spin_constrain.cpp — target_mag_in[iat].z change needs matching UnitCell fix
- **File:** `source/source_lcao/module_deltaspin/spin_constrain.cpp` (line 553)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629419
- **Issue:** Changing `target_mag_in[iat].x` → `.z` is correct only if the producer side (`UnitCell::get_target_mag()`) was also fixed to write the nspin==2 target in `.z`. Without that, every existing nspin==2 STRU file silently switches meaning.

### 3. cal_mw_from_lambda.cpp — GPU memory leak in update_psi_charge_pw_gpu
- **File:** `source/source_lcao/module_deltaspin/cal_mw_from_lambda.cpp` (line 68)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629427
- **Issue:** `h_tmp` and `s_tmp` are allocated via `resize_memory_op<...DEVICE_GPU>` but never freed. Each invocation leaks `2 * nbands^2` complex doubles on the device.
- **Fix:** Add matching `delete_memory_op<...DEVICE_GPU>` calls before return.

### 4. dftu_pw.cpp — Hubbard energy may be doubled for nspin==1
- **File:** `source/source_lcao/module_dftu/dftu_pw.cpp` (line 250)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629432
- **Issue:** For nspin==1, `weight_eu = 1.0` and the formula is `energy_u += u_value * weight_eu * locale * locale`. But `wg` already carries the factor-of-2 spin degeneracy. The PW formula may give `U * n^2 * 1.0` vs the LCAO reference of `U * n^2 / 2` — i.e. doubled.
- **Fix:** Verify against a known reference (single-atom Hubbard vs LCAO).

### 5. dftu_pw.cpp — cal_VU_pot_pw is a stub
- **File:** `source/source_lcao/module_dftu/dftu_pw.cpp` (line 352)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629439
- **Issue:** `cal_VU_pot_pw` is now a `(void)spin;` stub but the public declaration in `dftu.h` still advertises "calculate the local DFT+U effective potential matrix for PW base." Any caller expecting a side effect now silently does nothing.
- **Fix:** Either remove the public declaration or have it do the in-place update.

### 6. dftu_occup.cpp — copy_locale / mix_locale inconsistency
- **File:** `source/source_lcao/module_dftu/dftu_occup.cpp` (line 9)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629443
- **Issue:** The new `copy_locale` only writes `locale_save[iat][target_l][0][...]`. But `mix_locale` still iterates **all** `l` and **all** `n` and mixes `locale[iat][l][n][0/1]`. For non-target `(l, n)`, `locale_save` is stale.
- **Fix:** Either restrict `mix_locale` to the same target `(l, n=0)` slot, or keep `copy_locale` saving everything.

### 7. read_input_item_other.cpp — sc_mag_switch validator is commented out
- **File:** `source/source_io/module_parameter/read_input_item_other.cpp` (line 196)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629447
- **Issue:** The `sc_mag_switch` validator's `check_value` has all its nspin / calculation / nupdown enforcement commented out. A user can set `sc_mag_switch=true` with `nspin=1` and silently get nonsense.
- **Fix:** At minimum, restore the `nspin >= 2` requirement on `sc_mag_switch`.

### 8. read_input_item_other.cpp — New keywords lack validation, tests, and docs
- **File:** `source/source_io/module_parameter/read_input_item_other.cpp` (line 201)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629453
- **Issue:** Five new keywords (`sc_lambda_strategy`, `sc_direction_only`, `sc_scan_lambda_start`, `sc_scan_lambda_end`, `sc_scan_steps`) have no `check_value` validation, no test entries, and no user-facing documentation. Also, the previously-existing `sc_scf_nmin` keyword was deleted but `Input::sc_scf_nmin` (`input_parameter.h`) is still declared as an orphan.
- **Fix:** Add validation, tests, docs; remove orphan `sc_scf_nmin`.

### 9. module_deltaspin/test/CMakeLists.txt — Disabled test not re-enabled
- **File:** `source/source_lcao/module_deltaspin/test/CMakeLists.txt` (line 47)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629458
- **Issue:** `deltaspin_lambda_loop_helper_test` is disabled with a TODO "needs API update to match current SpinConstrain<TK> signature", but the PR **does** modify `lambda_loop_helper_test.cpp` itself.
- **Fix:** Either re-enable the test (and verify it passes), or drop the source edits.

### 10. module_dftu/test/CMakeLists.txt — 37 duplicate remove_definitions lines
- **File:** `source/source_lcao/module_dftu/test/CMakeLists.txt` (line 1)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629460
- **Issue:** 37 consecutive identical `remove_definitions(-D__CUDA)` lines — a clear copy-paste accident.
- **Fix:** Collapse to one line.

### 11. cusolver.h — No RAII guard for cudaMalloc
- **File:** `source/source_base/module_container/base/third_party/cusolver.h` (line 67)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629464
- **Issue:** New functions allocate `d_work` / `d_info` via `cudaMalloc` and free at the end. No RAII guard, so a `CHECK_CUSOLVER` failure between alloc and free leaks both.
- **Fix:** Add a scope guard, or check status code explicitly and free before throwing.

### 12. onsite_proj.cpp — becp cache is never read
- **File:** `source/source_pw/module_pwdft/onsite_proj.cpp` (line 412)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629471
- **Issue:** The new cache (`becp_ready_` / `ik_becp_`) is set/invalidated but never **read** anywhere. `force_onsite` / `stress_onsite` still recompute `becp` unconditionally.
- **Fix:** Either wire the consumer to skip recomputation, or remove the dead state.

### 13. force_op.h — Unnecessary #include of parameter.h in header
- **File:** `source/source_pw/module_pwdft/kernels/force_op.h` (line 3)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629476
- **Issue:** `#include "source_io/module_parameter/parameter.h"` in a kernel-op header causes compile-time bloat. Parameter values can be passed as function arguments.
- **Fix:** Move include into `.cpp` files only.

### 14. stress_op.h — Unnecessary #include of parameter.h in header
- **File:** `source/source_pw/module_pwdft/kernels/stress_op.h` (line 3)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629480
- **Issue:** Same as force_op.h: `#include "source_io/module_parameter/parameter.h"` pulls in the full parameter parser.
- **Fix:** Move include into `.cpp` files only.

### 15. elecstate_lcao.h — Unnecessary includes in header
- **File:** `source/source_estate/elecstate_lcao.h` (line 6)
- **URL:** https://github.com/deepmodeling/abacus-develop/pull/7384#discussion_r3316629485
- **Issue:** `parallel_orbitals.h` and `klist.h` are included in the header but only used in the `.cpp` file.
- **Fix:** Move them into the implementation file.

# DeltaSpin ~1000 eV Bug Investigation Plan

## Data Flow Overview

```
SCF Loop (esolver_ks_pw.cpp:hamilt2rho_single)
│
├── iter 1-8: HSolverPW::solve → psiToRho → normal (ETOT difference < 0.1 eV)
│
├── iter 9 (DS9): drho < sc_scf_thr → triggers DeltaSpin
│   │
│   └── run_deltaspin_lambda_loop(iter=8)
│       │
│       ├── [CP1] cal_mw_from_lambda(i_step=-1)  ← initial spin already has 3-4% deviation
│       │   ├── updateHk → cal_hs_subspace → save sub_h/sub_s/becp
│       │   ├── diag_responce → update ekb
│       │   ├── calculate_weights → update wg
│       │   └── accumulate_Mi_from_becp → Mi_
│       │
│       ├── [CP2] Lambda optimization loop (i_step=0..nsc)
│       │   ├── calculate_delta_hcc → h_tmp += delta_H(lambda)
│       │   ├── diag_responce → rotate becp, update ekb
│       │   ├── calculate_weights → update wg
│       │   ├── accumulate_Mi_from_becp → new spin
│       │   ├── cal_alpha_opt → alpha_opt (40% deviation!)
│       │   └── check_rms_stop → convergence check
│       │
│       ├── [CP3] update_psi_charge(dnu_last_step)  ← key to ~1000 eV energy jump!
│       │   └── update_psi_charge_pw_cpu:
│       │       ├── restore H/S from sub_h_save/sub_s_save
│       │       ├── calculate_delta_hcc(h_tmp, becp_save, delta_lambda=dnu_last_step)
│       │       ├── diag_subspace_psi → rotate psi, update ekb
│       │       ├── delete sub_h_save/sub_s_save/becp_save
│       │       └── psiToRho(*psi_t) → update charge density
│       │
│       └── [CP4] cal_mi_pw() → verify spin (Current RMS: 0.003 vs tmp 0.125)
│           ├── tabulate_atomic → overlap_proj_psi → get_h_becp
│           └── accumulate_Mi_from_becp
│
├── iter 10 (DS10): skip_solve=true (mag_converged)
│   └── skip HSolverPW::solve, use charge updated by CP3 directly
│       → Symmetry_rho → iter_finish → ETOT = -5348 eV (deviation ~1000 eV!)
│
└── iter 11+: continue with deviation
```

---

## Investigation Checklist

### Phase 1: Locate Initial Deviation Source Before SCF

**Objective**: Confirm the source of small spin deviations during DS1-DS8

- [ ] **1.1** Compare `ekb` at end of DS8 (eigenvalues of first 5 bands for all k-points)
  - Location: After `hsolver_pw_obj.solve` returns in `hamilt2rho_single`
  - Output: `[CHK-1.1] ik=%d band=%d ekb=%.10f`
  - Expected: Difference from tmp branch < 1e-6 Ry

- [ ] **1.2** Compare `wg` at end of DS8 (occupation numbers)
  - Location: Same as above
  - Output: `[CHK-1.2] ik=%d band=%d wg=%.10f`
  - Expected: Consistent with tmp branch

- [ ] **1.3** Compare `drho` value when DeltaSpin is triggered
  - Location: Entry of `run_deltaspin_lambda_loop`
  - Output: `[CHK-1.3] iter=%d drho=%.10e sc_scf_thr=%.10e`
  - Expected: Both versions trigger at the same iter

---

### Phase 2: Lambda Loop Initial Conditions

**Objective**: Confirm input and output of `cal_mw_from_lambda(i_step=-1)`

- [ ] **2.1** Check Hamiltonian subspace matrix `sub_h_save` after `updateHk`
  - Location: `cal_mw_from_lambda`, `initial_hs=1` branch, after `cal_hs_subspace`
  - Output: `[CHK-2.1] ik=%d h_diag[0..4]=%.10f` (diagonal elements)
  - Expected: Difference from tmp branch < 1e-8 Ry

- [ ] **2.2** Check `becp_save` (projection coefficients)
  - Location: After `memcpy(becp_k, onsite_p->get_becp(), ...)`
  - Output: `[CHK-2.2] ik=%d becp[0..3]=(%.8f,%.8f)` (first few complex numbers)
  - Expected: Consistent with tmp branch

- [ ] **2.3** Check `ekb` and rotated `becp_tmp` after `diag_responce`
  - Location: After `diag_responce` call
  - Output: `[CHK-2.3] ik=%d ekb[0..4]=%.10f`
  - Expected: Consistent with tmp branch

- [ ] **2.4** Check `wg` and `ef` after `calculate_weights`
  - Location: After `calculate_weights` call
  - Output: `[CHK-2.4] ik=%d wg[0..4]=%.10f ef=%.10f`
  - Expected: Consistent with tmp branch

- [ ] **2.5** Check input and output of `accumulate_Mi_from_becp`
  - Location: After loop ends, before/after `reduce_double_allpool`
  - Output: `[CHK-2.5] Mi_before_reduce[0]=(%.8f,%.8f,%.8f)` and `Mi_after_reduce`
  - Expected: For single process, before=after; difference from tmp branch < 1%

---

### Phase 3: Inside Lambda Optimization Loop

**Objective**: Locate the root cause of 40% `alpha_opt` deviation

- [ ] **3.1** Check input of first `calculate_delta_hcc`
  - Location: Inside `cal_mw_from_lambda(i_step=0, delta_lambda)`
  - Output: `[CHK-3.1] delta_lambda[0]=(%.8f,%.8f,%.8f) delta_lambda[1]=(%.8f,%.8f,%.8f)`
  - Expected: Consistent with tmp branch (since lambda starts at 0, delta_lambda is determined by search*alpha_trial)

- [ ] **3.2** Check `alpha_trial` and `search` vector
  - Location: In `run_lambda_loop`, after `check_restriction`
  - Output: `[CHK-3.2] i_step=%d alpha_trial=%.10f search[0]=(%.8f,%.8f,%.8f)`
  - Expected: Consistent with tmp branch

- [ ] **3.3** Check intermediate quantities `sum_k` and `sum_k2` in `cal_alpha_opt`
  - Location: Inside `cal_alpha_opt` (already has `[ALPHA-OPT]` print, need to enable `print=true`)
  - Output: `[CHK-3.3] sum_k=%.10e sum_k2=%.10e alpha_trial=%.10e`
  - Expected: Consistent with tmp branch

- [ ] **3.4** Check `spin_plus` (spin after alpha_trial step)
  - Location: After `cal_mw_from_lambda` returns, before `cal_alpha_opt`
  - Output: `[CHK-3.4] spin_plus[0]=(%.8f,%.8f,%.8f) spin_plus[1]=(%.8f,%.8f,%.8f)`
  - Expected: Difference from tmp branch < 1%

---

### Phase 4: update_psi_charge (Key to Energy Jump)

**Objective**: Confirm correctness of each step in `update_psi_charge_pw_cpu`

- [ ] **4.1** Check incoming `delta_lambda` (i.e., `dnu_last_step`)
  - Location: Entry of `update_psi_charge_pw_cpu`
  - Output: `[CHK-4.1] dnu_last_step[0]=(%.8f,%.8f,%.8f) dnu_last_step[1]=(%.8f,%.8f,%.8f)`
  - Expected: Consistent with tmp branch

- [ ] **4.2** Check diagonal elements of restored `sub_h_save`
  - Location: After `memcpy(h_tmp, h_k, ...)`
  - Output: `[CHK-4.2] ik=%d h_tmp_diag[0..4]=%.10f`
  - Expected: Consistent with CP2.1 values (sub_h_save saved at i_step=-1, unchanged thereafter)

- [ ] **4.3** Check `h_tmp` diagonal element change after `calculate_delta_hcc`
  - Location: After `calculate_delta_hcc` call
  - Output: `[CHK-4.3] ik=%d h_tmp_diag_after[0..4]=%.10f delta=%.10e`
  - Expected: Delta magnitude reasonable (consistent with lambda magnitude, ~0.1 Ry)

- [ ] **4.4** Check `ekb` after `diag_subspace_psi`
  - Location: After `diag_subspace_psi` call
  - Output: `[CHK-4.4] ik=%d ekb_after[0..4]=%.10f`
  - Expected: Consistent with tmp branch

- [ ] **4.5** Check total charge density after `psiToRho`
  - Location: After `psiToRho` call
  - Output: `[CHK-4.5] nelec_from_rho=%.10f` (integral of rho)
  - Expected: Equals total system electron count (16 for 2 Fe)

- [ ] **4.6** Compare behavior difference between `diag_subspace_psi` vs `diag_responce`
  - Key question: `diag_responce` rotates becp but not psi; `diag_subspace_psi` rotates psi
  - Check: Are H/S matrices used by both functions identical?
  - Output: `[CHK-4.6] Are h_tmp and s_tmp the same in both functions`

---

### Phase 5: psiToRho and Energy Calculation

**Objective**: Confirm correctness of energy calculation after charge density update

- [ ] **5.1** Check if `wg` is correct inside `psiToRho`
  - Location: Entry of `psiToRho`
  - Output: `[CHK-5.1] sum(wg)=%.10f` (should equal total electron count)
  - Expected: sum(wg) = 16.0

- [ ] **5.2** Check positivity and magnitude of rho after `psiToRho`
  - Location: After `psiToRho` returns
  - Output: `[CHK-5.2] rho_min=%.10e rho_max=%.10e rho_sum=%.10f`
  - Expected: rho_min >= 0, rho_sum ≈ 16

- [ ] **5.3** Check energy components in `iter_finish` for DS10
  - Location: After `cal_delta_eband` in `iter_finish`
  - Output: `[CHK-5.3] eband=%.6f deband=%.6f etxc=%.6f ehart=%.6f`
  - Expected: Each component has reasonable magnitude

---

### Phase 6: cal_mi_pw Verification (Abnormal Current RMS)

**Objective**: Understand why current=0.003 vs tmp=0.125

- [ ] **6.1** Check projection basis after `tabulate_atomic` in `cal_mi_pw`
  - Location: Inside `cal_mi_pw`, before `overlap_proj_psi`
  - Output: `[CHK-6.1] ik=%d nkb=%d`
  - Expected: nkb consistent with lambda loop

- [ ] **6.2** Compare becp from `cal_mi_pw` and lambda loop `accumulate_Mi_from_becp`
  - Key question: Lambda loop uses becp rotated by `diag_responce`; `cal_mi_pw` recomputes via `overlap_proj_psi`
  - If `diag_subspace_psi` correctly rotated psi, both should be consistent
  - Output: `[CHK-6.2] Mi_from_lambda_loop=(%.6f,%.6f,%.6f) Mi_from_cal_mi_pw=(%.6f,%.6f,%.6f)`

- [ ] **6.3** Check if `update_psi_charge_pw_cpu` correctly updated psi
  - Key: Did `diag_subspace_psi` actually rotate the wavefunction in `psi_t`?
  - Verify: Print psi norm before and after `diag_subspace_psi`
  - Output: `[CHK-6.3] ik=%d psi_norm_before=%.10f psi_norm_after=%.10f`

---

### Phase 7: Version Difference Investigation

**Objective**: Confirm behavioral changes in key functions between v3.7.1→v3.11.0

- [ ] **7.1** Compare `diag_subspace_psi` implementation
  - File: `source/source_hsolver/diago_iter_assist.cpp`
  - Check: Any API changes, parameter meaning changes

- [ ] **7.2** Compare `psiToRho` implementation
  - File: `source/source_estate/elecstate_pw.cpp`
  - Check: Any changes in `wg` usage

- [ ] **7.3** Compare `calculate_weights` implementation
  - File: `source/source_estate/elecstate_tools.cpp`
  - Check: Any changes in Fermi level search algorithm

- [ ] **7.4** Compare `cal_hs_subspace` implementation
  - File: `source/source_hsolver/diago_iter_assist.cpp`
  - Check: Subspace H/S matrix construction method

---

## Execution Priority

1. **Phase 4 (CP3)** — Direct cause of energy jump, most likely to find the bug
2. **Phase 6** — Abnormal Current RMS suggests psi update may have issues
3. **Phase 2** — Initial condition deviation propagates to all subsequent steps
4. **Phase 7** — If Phase 4/6 reveals no logic errors, then underlying function behavior has changed
5. **Phase 5** — If charge density itself is problematic
6. **Phase 3** — alpha_opt deviation may be a result of initial condition deviation
7. **Phase 1** — Lowest priority, small deviations may be normal numerical noise

---

## Quick Verification Plan

Before adding full diagnostics, perform these quick checks first:

### Quick Check A: Is sub_h_save correctly preserved?

Add at entry of `update_psi_charge_pw_cpu`:
```cpp
printf("[QUICK-A] sub_h_save[0]=(%.10f,%.10f) sub_s_save[0]=(%.10f,%.10f)\n",
       sub_h_save[0].real(), sub_h_save[0].imag(),
       sub_s_save[0].real(), sub_s_save[0].imag());
```

### Quick Check B: Is dnu_last_step reasonable?

Add where `check_rms_stop` returns true:
```cpp
for(int ia=0; ia<nat; ia++)
    printf("[QUICK-B] dnu_last_step[%d]=(%.8f,%.8f,%.8f)\n",
           ia, dnu_last_step[ia].x, dnu_last_step[ia].y, dnu_last_step[ia].z);
```

### Quick Check C: Electron number conservation after psiToRho

Add after `psiToRho` in `update_psi_charge_pw_cpu`:
```cpp
double nelec = 0;
for(int is=0; is<PARAM.inp.nspin; is++)
    for(int ir=0; ir<pelec->charge->nrxx; ir++)
        nelec += pelec->charge->rho[is][ir];
nelec *= ucell.omega / pelec->charge->nrxx;
printf("[QUICK-C] nelec_from_rho=%.10f expected=%.1f\n", nelec, pelec->nelec);
```

### Quick Check D: ekb before and after diag_subspace_psi

```cpp
printf("[QUICK-D] ik=%d ekb_before: %.10f %.10f %.10f\n", ik,
       pelec->ekb(ik,0), pelec->ekb(ik,1), pelec->ekb(ik,2));
// ... diag_subspace_psi ...
printf("[QUICK-D] ik=%d ekb_after:  %.10f %.10f %.10f\n", ik,
       pelec->ekb(ik,0), pelec->ekb(ik,1), pelec->ekb(ik,2));
```

---

## Key Hypotheses and Verification

| # | Hypothesis | Verification Method | If True |
|---|------|----------|----------|
| H1 | `diag_subspace_psi` did not correctly rotate psi | Phase 6.3: Check psi norm | psi not updated → psiToRho uses old psi |
| H2 | `dnu_last_step` has wrong value when passed to `update_psi_charge` | Quick-B | Lambda loop converged but wrong parameter passed |
| H3 | `sub_h_save` was accidentally modified during loop | Quick-A + Phase 4.2 | H matrix no longer has initial values |
| H4 | `wg` not conserved after `calculate_weights` | Phase 5.1 | Electron count wrong → charge density incorrect |
| H5 | `psiToRho` implementation behavior changed in v3.11 | Phase 7.2 | Need to adapt to new interface |
| H6 | `cal_mi_pw` and lambda loop use different psi | Phase 6.2 | psi state inconsistent after update_psi_charge |

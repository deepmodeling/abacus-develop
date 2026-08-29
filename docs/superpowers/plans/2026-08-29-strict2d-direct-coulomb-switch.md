# Strict-2D direct Coulomb switch implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan.

**Goal:** Productize the WS2 direct mixed-Fourier strict-2D Coulomb route as a fail-closed, default-off ABACUS reader-v1 option.

**Architecture:** Preserve the current Ewald producer path for exchange, Cs, shrink, and analytic head normalization.  When the explicit direct method is selected, construct positive Gram matrices from auxiliary Fourier amplitudes and overwrite only the reader-v1 full-Coulomb payload, with a separate provenance sidecar.  Validate input combinations before calculation and keep the old route unchanged by default.

**Tech Stack:** C++17, ABACUS input registry, LibRI auxiliary bases, MPI, reader-v1 binary output, GoogleTest/CTest, Slurm on df.

---

### Task 1: Add failing input-contract tests

**Files:**
- Modify: `source/source_io/test_serial/read_input_item_test.cpp`
- Test: `source/source_io/test_serial/read_input_item_test.cpp`

1. Add tests for the default method and parameters.
2. Add parsing tests for `direct_mixed_fourier`, cutoff 110 Ry, kz order 64, and Gamma order 8.
3. Add death tests for an unknown method, missing/negative direct parameters, odd Gamma order, non-2D Ewald, reader version 0, and `rpa=false`.
4. Build the serial input test remotely and record the expected RED failure because the new fields/items do not exist.

### Task 2: Add failing numerical and metadata tests

**Files:**
- Create: `source/source_lcao/module_ri/test/direct_2d_coulomb_test.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

1. Test positive quadrature weights and deterministic reciprocal enumeration.
2. Test the synthetic Gram accumulation for Hermiticity and nonnegative eigenvalues without clipping.
3. Test that increasing cutoff/order nests or converges the generated quadrature.
4. Test metadata formatting and exact method/parameter fields.
5. Build the target remotely and record the expected RED failure because the production header/source do not exist.

### Task 3: Implement the input interface

**Files:**
- Modify: `source/source_io/module_parameter/input_parameter.h`
- Modify: `source/source_io/module_parameter/read_input_item_output.cpp`
- Modify: `source/source_io/test_serial/read_input_item_test.cpp`

1. Add the four fields with legacy-preserving defaults.
2. Register the four input items and descriptions.
3. Add fail-closed cross-field validation for the direct method.
4. Run the focused input test remotely until GREEN.
5. Commit the input contract.

### Task 4: Implement the direct mixed-Fourier builder

**Files:**
- Create: `source/source_lcao/module_ri/direct_2d_coulomb.h`
- Create: `source/source_lcao/module_ri/direct_2d_coulomb.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`
- Test: `source/source_lcao/module_ri/test/direct_2d_coulomb_test.cpp`

1. Implement reciprocal enumeration, positive kz quadrature, even-order Gamma-plane quadrature, and Gram accumulation.
2. Expose numerical diagnostics without eigenvalue clipping.
3. Implement deterministic provenance formatting.
4. Run the focused numerical test remotely until GREEN.
5. Commit the numerical core.

### Task 5: Integrate reader-v1 output

**Files:**
- Modify: `source/source_lcao/module_ri/RPA_LRI.hpp`
- Modify: `source/source_lcao/module_ri/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

1. Keep the existing Ewald producer calculation and Cs/shrink flow intact.
2. When selected, write direct matrices using the exact full-grid q map and active auxiliary basis to `v1_coulomb_full_iq_*`.
3. Always retain the established strict-2D head sidecar.
4. Write `librpa_2d_coulomb_method.dat` on rank zero.
5. Add tests for q ordering, file sizes, method selection, and default-route compatibility.
6. Run the focused and module RI tests remotely until GREEN.

### Task 6: Remote build and regression

**Files:**
- Verify: remote staged source and build manifests

1. Synchronize committed source to a fresh df staging tree without touching the canonical build.
2. Configure with MPI, LibRI, LibComm, debug-info, and the validated df dependency profile.
3. Run focused CTest targets, then the relevant reader-v1/2D-Ewald regression set.
4. Verify `abacus_3p`, CMake cache flags, executable SHA256, and linked libraries.
5. Run a matched WS2 producer probe with the committed input interface and compare its matrices, exchange, and near-Fermi EXX band with the accepted prototype result.
6. Commit the integrated implementation only after fresh verification.

### Task 7: Controlled MoS2/WS2 bilayer production

**Files:**
- Create remotely: a fresh direct-method run root with `FORMAL_JOB_ID`, input hashes, and parent provenance
- Update locally: `.codex_tmp/mos2_ws2_gw_pbe_20260822/run-report.md`

1. Audit all existing job IDs and reject duplicate fingerprints.
2. Freeze the accepted structure, PP, NAO, stable external ABFS, symmetry, shrink, k mesh, and producer charge definition.
3. Change only the reader-v1 full-Coulomb method and record 110/64/8 as the initial WS2-derived test point.
4. Run `sbatch --test-only`, then submit exactly one producer without dependencies.
5. Require scheduler, application, reader-v1, exchange/Cs, Hermiticity, eigenvalue, and convergence gates.
6. Reuse the accepted producer-derived 121-point band only if its hashes and Gamma-state agreement pass; otherwise submit exactly one new matched band.
7. Generate and validate one full-grid PyATB payload, then one strict-2D LibRPA GW job sequentially.
8. Report low-energy GW bands, layer character, and reviewer-facing comparison only after numerical and physical acceptance.

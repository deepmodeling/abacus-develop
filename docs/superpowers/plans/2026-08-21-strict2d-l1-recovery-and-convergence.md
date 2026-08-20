# Strict-2D Analytic L1 Recovery and Convergence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the omitted analytic finite-q L=1 term in the ABACUS two-dimensional Ewald producer, prove that it reproduces the validated producer contract, and replace every affected MoS2, hBN, and phosphorene convergence result before reorganizing the research note.

**Architecture:** Keep the analytic two-dimensional Gaussian primitives in one header and route only the Coulomb L=1 channel through them; retain the existing GH32 path for L>=2. Use a focused unit test as the code gate, then one P4 10x14 and two MoS2 N20 producer/GW regressions as end-to-end gates before expanding any campaign. Treat every post-fix binary as a new immutable series.

**Tech Stack:** C++17, CMake/CTest, ABACUS reader-v1, LibRPA strict-2D head/wing, PyATB, Slurm, Python analysis, XeLaTeX.

---

### Task 1: Restore and test the analytic L=1 primitive

**Files:**
- Create: `source/source_lcao/module_ri/gaussian_abfs_2d_integrals.h`
- Create: `source/source_lcao/module_ri/test/gaussian_abfs_2d_integrals_test.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`
- Modify: `source/source_lcao/module_ri/gaussian_abfs.cpp`

- [ ] Add tests for the in-plane L=1 limit, the out-of-plane derivative identity, and finite small-q behavior.
- [ ] Configure/build the focused test and verify RED because the analytic header is absent.
- [ ] Add the analytic primitives and route `K_LM_2d(L=1,power=-2)` through `coulomb_l1`.
- [ ] Verify the focused test passes and the existing RI tests remain green.
- [ ] Compare the four modified source files against the validated 2026-08-04 source snapshot.

### Task 2: Commit, merge, and push master_ghj

**Files:**
- Modify: Git history only.

- [ ] Verify the feature worktree is clean except for the intended L=1 files and this plan.
- [ ] Commit with author `Codex <codex@openai.com>` and committer `AroundPeking <gonghuanjing@iphy.ac.cn>`.
- [ ] Rebase the feature branch onto the current `origin/master_ghj` and rerun tests.
- [ ] Fast-forward `master_ghj` in an isolated clean worktree and push it.
- [ ] Verify remote head, attribution, and clean status.

### Task 3: Build and run end-to-end producer gates on df_dcu

**Files:**
- Create: immutable build provenance and test receipts under `/work1/ghj/2d_gw_gth_multimaterial_validation_20260806/validation_20260821/`.
- Create: one producer/GW gate each for P4 10x14, MoS2 ONCV N20, and MoS2 GTH N20.

- [ ] Build ABACUS directly on the df_dcu login node from the pushed `master_ghj` commit and record executable SHA256.
- [ ] Run focused CTest gates and record source/binary provenance.
- [ ] Check the scheduler and ledgers to prevent duplicate submissions.
- [ ] Generate exactly one producer for each gate with reader-v1, symmetry, full Ewald, and the corrected L=1 path.
- [ ] Validate producer-local `stru_out`, LRI/shrink, symmetry, Coulomb coverage, and same-producer PyATB.
- [ ] Run one GW per gate with the pinned LibRPA binary and fixed Padé-12 for comparison where saved self-energy is reprocessed.
- [ ] Require P4 to reproduce the accepted 10x14 result within the stated numerical tolerance and MoS2 GTH to recover the CP2K-scale N20 gap.

### Task 4: Recompute only affected convergence data

**Files:**
- Create: post-fix ledgers, CSVs, plots, and validation records under the existing material campaign root.

- [ ] Audit P4 vacuum inputs and retain the existing series if its producers already contain the analytic L=1 implementation.
- [ ] Recompute only P4 mesh points produced by the incomplete ABACUS binary.
- [ ] Recompute MoS2 ONCV and GTH k/vacuum producer/GW points from the incomplete binary; use one continuation definition per curve.
- [ ] Recompute hBN k/vacuum producer/GW points from the incomplete binary and preserve the first-PBE-indexed conduction criterion.
- [ ] Generate smooth KS/GW overlays and on/off convergence plots with complete provenance.

### Task 5: Reorganize the research note

**Files:**
- Modify: `/Users/ghj/同步空间/AITP_project/2d_gw/derivations/main.tex`
- Modify/Create: material-specific figures and data under `/Users/ghj/同步空间/AITP_project/2d_gw/derivations/figures` and `data`.

- [ ] Create one complete MoS2 chapter containing setup, two pseudopotential chains, k convergence, vacuum convergence, final GW bands, CP2K comparison, and limitations.
- [ ] Create one complete hBN chapter containing setup, k convergence, vacuum convergence, final bands, and literature comparison.
- [ ] Create one complete phosphorene chapter containing setup, retained/recomputed vacuum evidence, k convergence, final bands, and literature comparison.
- [ ] Move development chronology, failed jobs, binary hashes, and detailed artifact contracts to appendices.
- [ ] Rebuild with XeLaTeX, require resolved references, render all pages, and visually inspect the three chapters and changed appendices.

### Task 6: Resume the four-TMDC CP2K comparison

**Files:**
- Create: converged GTH/SIAB campaign data and four overlaid GW band figures.

- [ ] Confirm every PBE gate and producer uses the post-fix ABACUS and pinned LibRPA profile.
- [ ] Establish k and vacuum convergence before selecting each final band calculation.
- [ ] Compare scalar-relativistic no-SOC bands against the CP2K companion data using the same plotting convention.
- [ ] Add four final overlays and a quantitative gap/Delta-v table to the MoS2/TMDC document section.

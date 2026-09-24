# SOC symmetry LibRPA export regression Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task with verification checkpoints.

**Goal:** Add a fast ABACUS regression test that protects the common `stru_out` spatial/magnetic symmetry export used by LibRPA, then validate it remotely before rebasing the branch.

**Architecture:** Move the existing symmetry-row serialization and near-integer validation into a small header-only `RpaLriDetail` helper. `RPA_LRI::out_struc` keeps all runtime gating and supplies ordinary and antiunitary symmetry arrays to that helper. A standalone module_ri test serializes deterministic operations and parses the output contract without PP, orbital, MPI, or full ABACUS runtime data.

**Tech Stack:** C++11-compatible ABACUS headers, CMake/CTest, remote GCC/MPI build on `df_iopcas_ghj`, Git.

---

### Task 1: Add the failing symmetry-format test

**Files:**
- Create: `source/source_lcao/module_ri/test/librpa_stru_symmetry_test.cpp`
- Modify: `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Add a test that calls the intended helper API.**

Create a standalone test with deterministic identity, inversion, and antiunitary translation operations. It should include `../librpa_stru_symmetry.h`, serialize to `std::ostringstream`, split the first line and three rows, and exit nonzero unless:

```cpp
require(lines[0] == "3 row", "ordinary and antiunitary operations share one row block");
require(lines[1].find("1   0   0") != std::string::npos, "identity operation is first");
require(lines[3].find("1.250000000000000e-01") != std::string::npos,
        "antiunitary translation is preserved");
require(output.find("spin_symmetry") == std::string::npos,
        "legacy spin-specific trailer must not be emitted");
```

For each operation, parse whitespace-separated fields and require exactly 12 fields; parse fields 1–9 as integers and fields 10–12 as finite doubles. Register `MODULE_RI_librpa_stru_symmetry_test` with the same standalone `add_executable`, `install`, and `add_test` pattern as `MODULE_RI_librpa_stru_units_test`.

- [ ] **Step 2: Run the new target remotely before adding production code.**

On `df_iopcas_ghj`, configure/build only this target against the current source and run it. Expected result: compilation fails because `librpa_stru_symmetry.h` and its helper are not yet present; this confirms the test exercises the new behavior rather than existing code.

### Task 2: Extract the production formatter

**Files:**
- Create: `source/source_lcao/module_ri/librpa_stru_symmetry.h`
- Modify: `source/source_lcao/module_ri/RPA_LRI.hpp:20-35,793,3114-3145`

- [ ] **Step 1: Define the focused helper.**

The header will include `source_base/matrix3.h`, `source_base/vector3.h`, `<cmath>`, `<iomanip>`, `<ostream>`, `<stdexcept>`, `<string>`, and `<vector>`. In namespace `RpaLriDetail`, define a plain-data `LibRpaSymmetryOperation` containing nine rotation values and three translation values, plus:

```cpp
inline int checked_near_int(const double value, const std::string& context);
inline LibRpaSymmetryOperation make_librpa_symmetry_operation(
    const ModuleBase::Matrix3& rotation,
    const ModuleBase::Vector3<double>& translation);
inline void write_librpa_symmetry_rows(
    std::ostream& output,
    const std::vector<LibRpaSymmetryOperation>& unitary,
    const std::vector<LibRpaSymmetryOperation>& antiunitary);
```

The production call converts ABACUS `Matrix3`/`Vector3` entries into the plain records, then the formatter writes the combined count, all ordinary rows first, all antiunitary spatial rows second, and the existing 9-integer/3-scientific-field formatting. The helper emits no spin trailer and remains independently compilable by the fast unit test.

- [ ] **Step 2: Reuse the helper from `RPA_LRI::out_struc`.**

Include the new header from `RPA_LRI.hpp`, remove the old local `write_op` lambda, and call:

```cpp
std::vector<RpaLriDetail::LibRpaSymmetryOperation> unitary(symm.nrotk);
for (int isym = 0; isym < symm.nrotk; ++isym)
    unitary[isym] = RpaLriDetail::make_librpa_symmetry_operation(symm.gmatrix[isym], symm.gtrans[isym]);
std::vector<RpaLriDetail::LibRpaSymmetryOperation> antiunitary(n_anti);
for (int isym = 0; isym < n_anti; ++isym)
    antiunitary[isym] = RpaLriDetail::make_librpa_symmetry_operation(symm.gmatrix_anti[isym], symm.gtrans_anti[isym]);
RpaLriDetail::write_librpa_symmetry_rows(ofs, unitary, antiunitary);
```

Remove the old `RpaLriDetail::checked_near_int` definition from `RPA_LRI.hpp` so there is one implementation and one error contract.

- [ ] **Step 3: Re-run the focused target remotely.**

Expected result: `MODULE_RI_librpa_stru_symmetry_test` builds and prints a passing message. Run the existing `MODULE_RI_librpa_stru_units_test` and `MODULE_RI_librpa_bz_sampling_test` in the same remote build to catch include or namespace regressions.

### Task 3: Validate the producer-level SOC gate

**Files:**
- No source changes; use the existing remote BN case under `/data/home/df_iopcas_ghj/app/abacus/soc_symmetry_compare-20260924`.

- [ ] **Step 1: Run the post-change tight SOC comparison.**

Run the already-defined `lspinorb=1` BN inputs for `symmetry=1` and `symmetry=-1`, with the same executable and convergence settings used for the pre/post merge evidence. Extract `!FINAL_ETOT_IS`, irreducible-k-point counts, and wall times.

- [ ] **Step 2: Apply the acceptance gates.**

Require both symmetry modes to finish, absolute energy difference no greater than `1e-6 eV`, and the symmetry-enabled run to use fewer irreducible k points. Record that reduced-k-point floating-point summation can produce a small non-bitwise energy difference; do not claim exact equality.

### Task 4: Commit the regression change

**Files:**
- `source/source_lcao/module_ri/librpa_stru_symmetry.h`
- `source/source_lcao/module_ri/RPA_LRI.hpp`
- `source/source_lcao/module_ri/test/librpa_stru_symmetry_test.cpp`
- `source/source_lcao/module_ri/test/CMakeLists.txt`

- [ ] **Step 1: Inspect the diff and run static checks.**

Check formatting, include guards, no `spin_symmetry` production output, and no unrelated changes. Do not compile locally.

- [ ] **Step 2: Commit with the required attribution.**

```bash
GIT_AUTHOR_NAME='Codex' GIT_AUTHOR_EMAIL='codex@openai.com' \
GIT_COMMITTER_NAME='AroundPeking' GIT_COMMITTER_EMAIL='gonghuanjing@iphy.ac.cn' \
git commit -m 'test(rpa): protect SOC symmetry stru_out export'
```

Verify the commit shows `Author: Codex <codex@openai.com>` and `Committer: AroundPeking <gonghuanjing@iphy.ac.cn>`.

### Task 5: Rebase after all gates pass

**Files:**
- Git history only; preserve the working tree and create a backup ref before rewriting.

- [ ] **Step 1: Fetch and record rebase inputs.**

Fetch `upstream/develop`, record its commit and the current branch tip, and create a backup ref. Confirm the worktree is clean.

- [ ] **Step 2: Rebase using the smallest history-preserving strategy.**

Attempt the planned upstream rebase only after the regression commit and remote gates pass. Preserve the full `master_ghj` feature history; resolve conflicts in favor of current upstream file naming and separated EXX/RPA interfaces while reapplying the PR3/PR4 and regression commits. If Git reports a conflict set that would require dropping the user’s 2D/Delta-ST history, stop and report the exact conflict boundary rather than silently discarding it.

- [ ] **Step 3: Re-run static and remote gates on the rebased tip.**

Verify the branch is clean, inspect the final diff against `upstream/develop`, rerun the focused module_ri tests and the short BN SOC comparison remotely, then push the rebased branch and update `origin/master_ghj` only after all gates pass.

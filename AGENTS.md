# ABACUS Agent Instructions

This file is the entry point for AI agents, automated review tools, and human
contributors who want the short operational version of the ABACUS development
rules. Read the complete governance document before making or reviewing changes:

- `docs/developers_guide/agent_governance.md`

## Required Baseline

- Follow the nine ABACUS coding rules summarized from the project governance:
  1. Do not increase cross-layer control through `GlobalV`, `GlobalC`, or
     `PARAM`; pass dependencies explicitly where practical. Migration-neutral
     moves must keep the PR-level global dependency budget non-increasing and
     explain the remaining global usage.
  2. Do not hide workflow switches in mutable member variables that can be
     changed from multiple places.
  3. Keep header dependencies minimal.
  4. Avoid adding `.hpp` implementation headers or propagating them through
     other headers unless there is a narrow reason.
  5. Do not add default arguments to existing interfaces; update call sites or
     design a clearer extension.
  6. Add focused tests for key features, bug fixes, INPUT behavior changes,
     heterogeneous kernels, and core-module refactors.
  7. Keep code compatible with the repository C++11 baseline.
  8. Declare one variable per line; do not use comma-separated declarations.
  9. Do not call MPI routines directly; use the internally-guarded wrappers
     (e.g., `Parallel_Reduce::reduce_*`, `Parallel_Common::bcast_*`) instead.
  10. Do not write new `#define private public` or `#define protected public`
      access hacks in test files. If a unit test needs to inspect internal
      state, either promote the member visibility explicitly or add a
      public test-only accessor.
  11. New unit test source files shall be named `test_<module_name>.cpp`,
      matching the source file they exercise. For example, the test for
      `rhog_io.cpp` shall be `test_rhog_io.cpp`. This naming keeps the
      file-to-test relationship discoverable and consistent across the
      repository. Historical tests are not required to be renamed.
- Use LF line endings for text files. Only `.bat` and `.cmd` files may use CRLF.
- Keep source file additions deterministic: update the relevant `CMakeLists.txt`
  or explain why the file is generated or included indirectly.
- INPUT parameter behavior changes must update `docs/parameters.yaml` and
  `docs/advanced/input_files/input-main.md`, or the PR must state why no update
  is required.
- Report the exact verification performed, and match the depth of verification
  to the risk of the change. `Verification Tiers` below states which evidence
  each kind of change owes. Never present an unrun command as if it had passed,
  and never silently claim a lower tier than the change belongs to.
- Prefer `std::vector` over raw `new`/`delete` for dynamic arrays; before
  converting class members, confirm no external code consumes them as raw
  pointers (e.g., `std::vector<bool>` has no `.data()`), and use
  `std::fill`/`std::copy` instead of `ZEROS`/`COPYARRAY` on vector buffers.

## Repository Map

- Core C++ implementation lives under `source/`; source additions must be wired
  through the relevant `CMakeLists.txt`.
- INPUT parsing and help metadata live under `source/source_io/`; user-facing
  INPUT docs live in `docs/parameters.yaml` and
  `docs/advanced/input_files/input-main.md`.
- Unit tests are colocated under module `test/` directories such as
  `source/source_md/test/`; integration and workflow tests are selected through
  CTest labels and patterns.
- Developer and user build/install references live in `docs/quick_start/`,
  `docs/advanced/`, `toolchain/`, `Dockerfile.gnu`, `Dockerfile.intel`, and
  `Dockerfile.cuda`.

## Build And Test Entry Points

- Prefer the repository CMake/CTest flow already used by CI. For focused local
  checks, use commands such as `ctest --test-dir build -V -R MODULE_MD` after a
  usable build exists.
- The executable is not named `abacus`. CMake names it
  `abacus_<level>_<para-or-device>` (`ABACUS_BIN_NAME` in the top-level
  `CMakeLists.txt`), where level is `pw`, `basic`, `ml`, `std`, or `max`, and
  `install` keeps that name. Locate it the way `.github/workflows/test.yml`
  does, and use the resolved path in recorded commands:

  ```bash
  ABACUS_BIN=$(find build -name "abacus_*" -type f -executable | head -1)
  ```

- For INPUT-related changes, verify both documentation and CLI behavior when an
  executable is available: `"$ABACUS_BIN" -h <parameter>` and
  `"$ABACUS_BIN" --check-input` from a valid case directory.
- For executable identity, record `"$ABACUS_BIN" --version` together with the
  binary name it resolved to.
- Reuse existing Docker and toolchain assets. Do not add a new container,
  compiler setup, or calculation-task skill unless the PR explicitly requires
  and justifies it.

## Verification Tiers

The verification a change owes depends on what the change can break, not on how
large the diff is. Pick the highest tier that applies, produce its evidence, and
record the exact commands in the PR.

- **Tier 1, static checks only.** Text no compiler reads: documentation,
  comments, PR templates, governance files, whitespace or line-ending
  normalization. Evidence: the governance check from `Local Commands` below,
  plus regenerated `docs/parameters.yaml` and
  `docs/advanced/input_files/input-main.md` when those artifacts are in scope.
  A local build is not expected, and configuring one to satisfy this tier is
  wasted effort.
- **Tier 2, must compile.** Changes that can break the build without changing
  any number: include and header-dependency edits, moving code between
  translation units, splitting or merging files, renames, mechanical call-site
  updates, interface signature changes, `new`/`delete` to `std::vector`
  conversions, and `CMakeLists.txt` linkage. Evidence: a build of the affected
  targets plus the unit tests of the touched modules, such as
  `ctest --test-dir build -R MODULE_MD`. Header-dependency trimming belongs
  here even when it only deletes lines, because a removed include can break an
  unrelated translation unit that relied on it transitively.
- **Tier 3, must reproduce behavior.** Changes that can move numbers or alter
  control flow: algorithms, operators, heterogeneous kernels, eigensolvers,
  mixing, occupations, symmetry, INPUT semantics or defaults. Evidence: Tier 2,
  plus at least one relevant case from `tests/`, plus the CLI checks above for
  INPUT and command-line changes. Say which reference values were compared, not
  only that the case ran.

For a multi-step Tier 2 or Tier 3 refactor, such as splitting a large `.cpp`
into several files, build and commit after each step rather than batching every
change before the first build. This keeps the blast radius small when a step
surfaces a missing include or instantiation error.

When the evidence a tier requires cannot be produced locally, say so
specifically: name the tier, name the command that could not be run, and give
the reason, such as no Linux toolchain in the checkout, no configured build
tree, a sandbox restriction, or a missing optional dependency. ABACUS builds on
Linux in practice, so a Windows or otherwise unsupported checkout is a
legitimate reason to stop at static checks and leave compile and runtime
evidence to CI. A disclosed gap is acceptable and reviewable. An undisclosed gap
is not, and neither is provisioning a new toolchain or container to reach a tier
the change never required; reuse the existing Docker and toolchain assets or
defer to CI.

## Local Runtime Testing

- Set `OMP_NUM_THREADS=1` for ABACUS runtime, integration, and MPI tests unless
  a test explicitly requires another value.
- Run MPI/runtime tests outside restricted sandboxes when process visibility,
  sockets, or MPI launch behavior matters.
- Treat OpenMPI `opal_ifinit: socket() failed errno=1` warnings from sandboxed
  MPI-linked builds or runs as expected sandbox artifacts; rerun outside the
  sandbox before treating them as ABACUS failures.
- Do not relax existing tests or references merely to make a failure pass.
  Update references only when the intended behavior changed and the PR explains
  why.

## Review And Exception Flow

- Mechanical blockers are enforced by hook and CI only for new files, changed
  files, or diff-added lines. Historical untouched code is not a default blocker.
- Warnings from CI or AI review require reviewer attention but do not block by
  themselves.
- Semantic questions such as module ownership, member-variable workflow state,
  test sufficiency, and exception approval require human review.
- Exceptions must be recorded in the PR with reason, scope, risk, and a follow-up
  cleanup plan.
- After a refactor, propose brief lessons worth recording in this file, then
  ask the developer whether to write them in; be cautious and skip unclear
  or unverified lessons.

## Refactoring Patterns

- Member -> free function: inventory `this->` reads; pass as params (const
  for config, ref for mutable state); move only when body is `this`-free;
  keep thin wrapper; compile each step when a usable build exists (Tier 2).
- Extract a base-class nested-vector member in three steps (hold + forward,
  switch writers, delete legacy) so no commit mixes old-storage writes with
  new-storage reads.

## Local Commands

```bash
python3 tools/03_code_analysis/agent_governance_check.py --staged
python3 tools/03_code_analysis/agent_governance_check.py --base upstream/develop --head HEAD --format text
pre-commit run abacus-agent-governance --all-files
# Score changed C++ files for quality debt (pass line is 60):
python3 tools/03_code_analysis/code_quality_score.py $(git diff --name-only upstream/develop...HEAD | grep -E '\.(cpp|h)$')
```

The repository text files have been normalized to LF once. Day-to-day line
ending enforcement should rely on staged/changed-file hooks and CI; rerun the
full mixed-line-ending hook only for intentional repository-wide normalization.

## PR Self-Check

- Confirm the PR body names the verification tier reached, states the exact
  commands run and whether they passed or failed, and explains why any check
  expected at that tier could not be run.
- Keep warning rationales concrete. For example, a header include warning can be
  acceptable when the header owns a value member that requires the complete type.
- Keep historical-debt notes separate from new deterministic errors introduced
  by the PR.

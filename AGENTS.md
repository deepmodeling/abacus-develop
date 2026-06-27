# ABACUS Agent Instructions

This file is the entry point for AI agents, automated review tools, and human
contributors who want the short operational version of the ABACUS development
rules. Read the complete governance document before making or reviewing changes:

- `docs/developers_guide/agent_governance.md`

## Required Baseline

- Follow the seven ABACUS coding rules summarized from the project governance:
  1. Do not introduce new cross-layer control through `GlobalV`, `GlobalC`, or
     `PARAM`; pass dependencies explicitly.
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
- Use LF line endings for text files. Only `.bat` and `.cmd` files may use CRLF.
- Keep source file additions deterministic: update the relevant `CMakeLists.txt`
  or explain why the file is generated or included indirectly.
- INPUT parameter behavior changes must update `docs/parameters.yaml` and
  `docs/advanced/input_files/input-main.md`, or the PR must state why no update
  is required.
- Report the exact verification performed. Do not claim completion without
  fresh test or check output.

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
- For INPUT-related changes, verify both documentation and CLI behavior when an
  executable is available: `./build/abacus -h <parameter>` and
  `./build/abacus --check-input` from a valid case directory.
- For executable identity, record `./build/abacus --version` or the equivalent
  installed `abacus --version` command used during verification.
- Reuse existing Docker and toolchain assets. Do not add a new container,
  compiler setup, or calculation-task skill unless the PR explicitly requires
  and justifies it.

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

## Local Commands

```bash
python3 tools/03_code_analysis/agent_governance_check.py --staged
python3 tools/03_code_analysis/agent_governance_check.py --base upstream/develop --head HEAD --format text
pre-commit run abacus-agent-governance --all-files
```

The repository text files have been normalized to LF once. Day-to-day line
ending enforcement should rely on staged/changed-file hooks and CI; rerun the
full mixed-line-ending hook only for intentional repository-wide normalization.

## PR Self-Check

- Confirm the PR body states exact commands run, whether they passed or failed,
  and why any expected check could not be run.
- Keep warning rationales concrete. For example, a header include warning can be
  acceptable when the header owns a value member that requires the complete type.
- Keep historical-debt notes separate from new deterministic errors introduced
  by the PR.

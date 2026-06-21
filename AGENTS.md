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
pre-commit run abacus-agent-governance --all-files
```

The repository text files have been normalized to LF once. Day-to-day line
ending enforcement should rely on staged/changed-file hooks and CI; rerun the
full mixed-line-ending hook only for intentional repository-wide normalization.

# ABACUS Agent Governance

This document is the shared development and review contract for human
contributors, general AI agents, CodeRabbit, and GitHub CI. `AGENTS.md` is the
short entry point; this file is the complete rule source for implementation and
review.

## Source Materials

The governance rules consolidate project guidance from @mohanchen's coding
rules, the developer guide, contribution guide, PR template, existing CI, and
historical Chinese development notes. Historical notes are used only after path
modernization and rule grading; they are not copied directly into automated
checks.

## Diff Scope

Automated checks default to the PR diff:

- New files.
- Diff-added lines.
- Newly introduced symbols or includes.
- Changed text files for line-ending checks.

Untouched historical code is not a default blocker. Review tools may mention
historical debt when it is relevant to the changed area, but they must separate
that from blocking findings on new changes.

## Core Coding Rules

- Do not introduce new cross-layer control through `GlobalV`, `GlobalC`, or
  `PARAM`. Prefer explicit parameters or narrow local interfaces.
- Do not store workflow switches in mutable member variables that can be changed
  implicitly from multiple places.
- Keep header dependencies minimal; avoid adding includes to headers unless the
  declaration truly requires them.
- Avoid adding `.hpp` implementation headers and avoid including `.hpp` files
  from other headers unless the PR explains why header-only implementation is
  needed.
- Do not add default arguments to existing interfaces. Update call sites
  explicitly or design a clearer overload/configuration object.
- Add short, focused tests for key functionality, bug fixes, INPUT behavior
  changes, heterogeneous kernels, and core-module refactors.
- Keep default and general-purpose C++ changes compatible with the repository
  C++11 baseline. Backend-specific or dependency-constrained paths may use the
  higher standard already selected by existing CMake configuration.
- Use LF line endings for text files. `.bat` and `.cmd` are the CRLF exception.

AI agents have additional workflow obligations:

- Inspect existing interfaces before using them.
- State uncertainty instead of inventing business rules or APIs.
- Report exact verification results and any checks that could not be run.

## Rule Grading Matrix

The first implementation phase separates deterministic mechanical checks from
review-only rules. "Phase-one mechanical" means the local hook or CI checker can
act on the PR diff without semantic judgment. "AI review" and "human
confirmation" items must be reviewed, but the checker does not hard-code those
decisions.

| Rule category | Typical rule | Phase-one status | Default executor | Severity | Default action | Detection scope | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Basic text format | LF line endings | phase-one mechanical | hook + CI | medium | block | full changed text file | `.bat` and `.cmd` keep CRLF |
| Language baseline | C++11 compatibility | build/toolchain | CI | high | block | build/static tooling | Actual compiler/toolchain result wins |
| New global dependency | Added `GlobalV`/`GlobalC`/`PARAM` as cross-layer control | phase-one mechanical + AI review | CI + AI review | high | block | added code lines | Historical untouched usage and documentation mentions are not blocked |
| New default parameter | Header declaration adds a default argument | phase-one mechanical + AI review | CI + AI review | high | block | header diff | High misuse risk |
| `.hpp` propagation | New `.hpp` or header includes `.hpp` | phase-one mechanical warning | CI + AI review | medium | warn | new files and added includes | Exception can be recorded in PR |
| Header dependency growth | Header diff adds includes | phase-one mechanical warning + AI review | CI + AI review | medium | warn | added header includes | Necessity is semantic and not mechanically decided |
| Member variable workflow switch | Key flow state hidden as mutable member state | AI review + human confirmation | AI + human review | high | human confirmation | semantic review | Static matching is unreliable |
| Module path and build linkage | New source path and `CMakeLists.txt` linkage | phase-one mechanical | CI | medium | block | new source files and build-script diff | Deterministic path/build check only |
| Module semantic ownership | Best module/submodule placement | AI review + human confirmation | AI + human review | medium | human confirmation | semantic review | Final call belongs to maintainers |
| Heterogeneous code linkage | CUDA/ROCM/kernel source and `CMakeLists.txt` linkage | phase-one mechanical | CI + AI review | medium | block | new heterogeneous files and linkage | Mechanical path/linkage only |
| Heterogeneous test evidence | CUDA/ROCM/kernel change has test evidence or reason | phase-one mechanical warning + AI review | CI + AI review | high | warn | changed paths and PR body | Sufficiency is human-reviewed |
| Test existence | Source change has test evidence or reason | phase-one mechanical warning + AI review | CI + AI review | high | warn | PR body and changed paths | Sufficiency is human-reviewed |
| Test sufficiency | Tests cover important behavior | AI review + human confirmation | AI + human review | medium | human confirmation | semantic review | Not mechanically blocked |
| INPUT behavior linkage | Parameter metadata/default/type/parser behavior updates YAML and docs | phase-one mechanical + AI review | CI + AI review | high | block | behavior-field diff plus docs/PR body | Comment-only parameter-file changes are not blocked |
| Documentation sync | Behavior/interface docs updated | phase-one mechanical warning + AI review | CI + AI review | medium | warn | changed paths and PR body | Major behavior changes escalate to reviewers |
| PR metadata completeness | Issue, tests, behavior, INPUT, core impact, exceptions | phase-one mechanical | CI or GitHub bot | medium | block | PR template fields | Not run by local hook |
| AI workflow | Interface lookup, uncertainty, verification report | AI review | AI review | high | warn | review transcript/output | Applies to AI agents |
| Exceptions | Reason, scope, risk, follow-up plan | human confirmation | human review + CI | high | human confirmation | PR exception section | CI checks presence, not approval |

Implementations must use this matrix as the default baseline. Upgrading warnings
to blockers or converting human-confirmation rules into mechanical blockers
requires an explicit governance change.

## Automation Responsibilities

Local hooks:

- Fix or block deterministic local issues such as mixed line endings.
- Run staged governance checks with `agent_governance_check.py --staged`.
- Do not require PR metadata because it is unavailable locally.

CI:

- Run diff-level governance checks with PR base/head SHAs.
- Check PR body completeness and INPUT documentation linkage.
- Publish Markdown summaries that humans and AI reviewers can consume.

AI review:

- Explain governance findings in actionable terms.
- Add semantic review for module ownership, header dependency growth, test
  sufficiency, documentation sync, and AI workflow discipline.
- Use the output format below for actionable findings.

Human review:

- Approve or reject exceptions.
- Confirm module boundaries and high-risk design decisions.
- Decide whether tests are sufficient for the scientific and numerical risk.

## AI PR Review Integration

ABACUS uses a layered review model:

- `Agent Governance` is the deterministic GitHub Actions check for low-noise
  diff rules. Repository maintainers may make this workflow a required check in
  branch protection.
- CodeRabbit is the minimal PR-triggered AI reviewer for semantic review hints.
  Its repository configuration lives in `.coderabbit.yaml` and uses this
  document plus `AGENTS.md` as review guidelines.
- CodeRabbit comments are advisory by default. They do not replace maintainer
  approval, exception approval, or numerical/test sufficiency review.

To activate CodeRabbit on real PRs, a repository or organization administrator
must install the CodeRabbit GitHub App for `deepmodeling/abacus-develop` and
grant it pull request review access. After installation, non-draft pull requests
and new commits should receive automatic review according to `.coderabbit.yaml`;
if automatic review does not start, maintainers may request it with
`@coderabbitai review`.

Copilot code review, Copilot coding-agent setup steps, Qodo, and PR-Agent are
not part of the phase-one baseline. They require separate organization settings,
secrets, or setup workflows and should be added only through a later governance
change.

## INPUT Parameter Changes

Changes to parameter metadata, default values, type, availability, description,
or parsing behavior must include both:

- `docs/parameters.yaml`
- `docs/advanced/input_files/input-main.md`

If the diff touches parameter internals but does not change user-visible INPUT
behavior, the PR must state why no documentation update is required.

## Exception Template

Use this template in the PR when a rule must be bypassed temporarily:

```markdown
### Governance Exception
- Rule:
- Reason:
- Scope:
- User or maintenance risk:
- Why the normal rule cannot be followed now:
- Follow-up cleanup plan:
- Requested approver:
```

## AI Review Finding Format

AI and bot review findings should use this shape:

```markdown
Rule: <rule name>
Severity: error | warning | info
Location: <file:line>
Reason: <specific trigger>
Suggested action: <concrete next step>
Exception: allowed | not allowed | human approval required
```

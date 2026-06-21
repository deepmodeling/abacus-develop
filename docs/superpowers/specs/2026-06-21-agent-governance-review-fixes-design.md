# Agent Governance Review Fixes Design

## Context

The governance checker added in `tools/03_code_analysis/agent_governance_check.py`
is intended to block deterministic issues only on changed files or diff-added
lines. Three review findings were verified locally:

- CRLF-only normalization diffs are treated as semantic added lines by
  `git diff -U0`, so historical default arguments and global dependencies can
  be reported as new changes.
- A pull request event with an empty `pull_request.body` is treated the same as
  no pull request event, so required PR metadata can be skipped.
- The GitHub Actions workflow passes the PR base tip directly as `--base`.
  After the target branch advances, that two-dot comparison no longer matches
  the merge-base-to-head PR diff that GitHub presents.

## Goals

- Keep semantic governance rules scoped to meaningful diff-added lines.
- Keep full changed-file LF enforcement intact, including detection of CRLF in
  the target file content.
- Treat an empty PR description in a PR event as incomplete metadata.
- Scope CI governance checks to the PR merge base and head commit.
- Add focused regression tests for each verified failure mode.

## Non-Goals

- Do not redesign the checker CLI or change the meaning of `--base` and
  `--head`.
- Do not convert warning-only governance rules into blockers.
- Do not change the PR template or governance policy text beyond what is needed
  to align implementation behavior with the existing policy.

## Recommended Approach

Use the minimal complete fix:

1. Make `added_lines()` use a CR-at-EOL-insensitive diff:
   - For staged mode, run `git diff --cached --ignore-cr-at-eol -U0`.
   - For base/head mode, run `git diff --ignore-cr-at-eol -U0 <base> <head>`.
2. Leave `changed_paths()` and `check_line_endings()` unchanged so changed-file
   discovery and byte-level LF enforcement still see files whose content changed
   only by line endings.
3. Make PR body reading distinguish these cases:
   - No `--event-path`: skip PR metadata completeness, because local/staged runs
     do not have PR metadata.
   - Invalid or non-PR event payload: skip PR metadata completeness.
   - PR event with `body: ""` or `body: null`: run metadata completeness and
     block as incomplete.
4. In `.github/workflows/agent_governance.yml`, compute
   `merge_base=$(git merge-base "$base_sha" "$head_sha")` after checkout and pass
   that value as `--base`.

This keeps responsibilities simple: the workflow chooses the correct comparison
range for PR CI, while the checker continues to compare two explicit commits.

## Component Changes

### Semantic Diff Collection

`added_lines()` is the only source of diff-added lines for semantic rules such
as no new globals, no new default parameters, `.hpp` propagation, header include
growth, and INPUT behavior checks. Adding `--ignore-cr-at-eol` there prevents
line-ending normalization from fabricating added semantic lines.

The separate changed-file path list remains based on regular `git diff
--name-status`. This is important because line-ending-only changes still need to
be included in the full-file LF check.

### PR Metadata Reading

Introduce a small internal distinction between "metadata unavailable" and
"metadata available but empty". A simple option is to let `read_pr_body()` return
`Optional[str]`, where `None` means no PR metadata is available and `""` means
the PR body is blank.

`check_pr_metadata()` should skip only when the body is `None`. For an empty
string, it should report all required sections as missing, using the existing
`PR metadata completeness` finding and `allow_exception=False`.

Other PR-body consumers should handle `None` as an empty string for warning
logic, because test evidence and documentation evidence cannot be inferred when
metadata is unavailable.

### Workflow Diff Base

The workflow already uses `actions/checkout@v4` with `fetch-depth: 0`, so the
merge base can be computed locally:

```bash
base_sha="${{ github.event.pull_request.base.sha }}"
head_sha="${{ github.event.pull_request.head.sha }}"
merge_base="$(git merge-base "$base_sha" "$head_sha")"
python3 tools/03_code_analysis/agent_governance_check.py \
  --base "$merge_base" \
  --head "$head_sha" \
  --event-path "$GITHUB_EVENT_PATH" \
  --format markdown
```

If `git merge-base` fails, the workflow should naturally fail before publishing
a misleading clean governance result. The existing summary fallback can still
report that the checker failed before producing a summary.

## Error Handling

- Git command failures inside the checker should continue to surface through
  `GitError` and exit code 2.
- Invalid JSON or non-PR event payloads should not block local or accidental
  non-PR invocations.
- Blank PR bodies in valid PR payloads should block with the existing metadata
  completeness rule.

## Tests

Add focused tests to `tools/03_code_analysis/test_agent_governance_check.py`:

- A CRLF-to-LF-only change in a header containing a historical default argument
  should not trigger `No new default parameters`.
- A valid PR event with `pull_request.body` set to an empty string should trigger
  `PR metadata completeness` and return non-zero.
- A merge-base-scoped comparison in a diverged history should not include
  base-branch-only changes, while the existing two-dot base-tip comparison would
  demonstrate the risk.

Existing tests for CRLF detection in changed text files must continue to pass.

## Verification

After implementation, run:

```bash
python3 -m unittest tools/03_code_analysis/test_agent_governance_check.py
python3 tools/03_code_analysis/agent_governance_check.py --staged
```

If dependencies are available, also run:

```bash
pre-commit run abacus-agent-governance --all-files
```

Report exact command output or failure details before claiming completion.

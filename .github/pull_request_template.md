### Reminder
- [ ] I have read `AGENTS.md` and `docs/developers_guide/agent_governance.md`.
- [ ] I have linked an issue or explained why this PR does not need one.
- [ ] I have added adequate unit tests and/or case tests, or explained why not.
- [ ] I have listed the exact verification commands run and their results.
- [ ] I have described user-visible behavior changes, including INPUT parameter changes.
- [ ] I have explained core-module impact for ESolver, HSolver, ElecState, Hamilt, Operator, Psi, or other `source/` changes.
- [ ] I have requested any needed governance exception below.

### Linked Issue
Fix #...

### Unit Tests and/or Case Tests for my changes
- A unit test is added for each new feature or bug fix.

### Exact Verification Performed
- Commands run:
- Result summary:
- Checks not run, with reason:

### What's changed?
- Example: My changes might affect the performance of the application under certain conditions, and I have tested the impact on various scenarios...

### Governance Checklist
- Global dependencies: no net increase in `GlobalV`, `GlobalC`, or `PARAM` code references, or exception requested below with reason, scope, risk, and cleanup plan.
- Default parameters: no new default arguments added to existing interfaces, or exception requested below.
- Headers: no unnecessary header dependencies or `.hpp` propagation, or rationale provided below.
- Line endings: text files use LF; only `.bat` and `.cmd` use CRLF.
- Build linkage: new source files are listed in the relevant `CMakeLists.txt`, or rationale provided below.
- Documentation: behavior/interface changes include documentation updates, or no documentation update is required because ...
- CodeRabbit: if automatic review has not started and the repository has CodeRabbit installed, request `@coderabbitai review`.

### INPUT Parameter Changes
- Parameters added/removed/changed:
- `docs/parameters.yaml` updated: yes/no/not applicable
- `docs/advanced/input_files/input-main.md` updated: yes/no/not applicable
- If not updated, explain why no INPUT documentation update is required:

### Core Module Impact
- Affected core modules:
- Risk summary:
- Compatibility or performance impact:

### Governance Exception
- Rule:
- Reason:
- Scope:
- User or maintenance risk:
- Why the normal rule cannot be followed now:
- Follow-up cleanup plan:
- Requested approver:

---
applyTo: "**"
---

# ABACUS Governance Review Instructions

Apply these instructions when reviewing or changing ABACUS code:

- Use `AGENTS.md` and `docs/developers_guide/agent_governance.md` as the
  authoritative project baseline.
- Keep review scope diff-oriented: new files, diff-added lines, new includes,
  newly introduced symbols, and changed text files for line-ending checks.
- Do not treat untouched historical debt as a default blocker. Mention it only
  when it affects the changed area, and label it as advisory.
- Flag PRs that increase `GlobalV`, `GlobalC`, or `PARAM` code references as
  blocker-level governance issues. Flag migration-neutral added usage as a
  warning that requires reason, scope, risk, and cleanup/follow-up rationale.
  Prefer explicit dependencies or narrow local interfaces.
- Flag new default arguments in existing header interfaces. Prefer explicit
  call-site updates, overloads, or a clearer configuration object.
- Review header include growth and `.hpp` propagation carefully. These are
  usually warnings unless the PR records a narrow reason.
- Require LF line endings for text files. `.bat` and `.cmd` files are the CRLF
  exceptions.
- For INPUT parameter behavior changes, require synchronized updates to
  `docs/parameters.yaml` and `docs/advanced/input_files/input-main.md`, or a
  clear no-update explanation in the PR.
- Check that new source files are linked through the relevant `CMakeLists.txt`
  unless the PR explains generated or indirect inclusion.
- Keep default C++ changes compatible with the repository C++11 baseline.
- Ask for focused tests or explicit test rationale for feature changes, bug
  fixes, INPUT behavior changes, heterogeneous kernels, and core-module
  refactors.
- Treat CI governance findings as deterministic evidence and semantic review
  findings as advisory until maintainers approve them.

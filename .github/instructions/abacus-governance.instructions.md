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
- Review header include growth and `.hpp` propagation carefully. Do not ask
  authors to replace explicit `#include` directives with forward declarations
  solely to avoid header include warnings.
- Require LF line endings for text files. `.bat` and `.cmd` files are the CRLF
  exceptions.
- Flag new `#pragma once` in header files as a blocker-level governance issue.
  Require traditional `#ifndef`/`#define`/`#endif` include guards instead.
- For INPUT parameter behavior changes, require corresponding C++ `Input_Item`
  metadata updates, or a clear no-update explanation in the PR. Do not require
  generated YAML or Markdown documentation files to be committed.
- Check that new source files are linked through the relevant `CMakeLists.txt`
  unless the PR explains generated or indirect inclusion.
- Apply `docs/developers_guide/cmake_guidelines.md` for CMake edits: require
  `cmake-format` on changed CMake files; flag newly added `${math_libs}`,
  `USE_LIBXC`, `ENABLE_RAPIDJSON`, `ENABLE_LIBCOMM`, `remove_definitions()`,
  directory-scoped compile/include/link or `find_package`/`find_library`
  commands under `source/`, and `install()` or `file(COPY ...)` in unit-test
  CMake files as blocking; and confirm new options or external dependencies
  update `_abacus_feature_definitions`, `abacus_add_target_compile_requirements()`,
  and `cmake/CollectBuildInfoVars.cmake`.
- Unit tests under `source/` belong in a flat `unittests/` directory, with
  sources named `test_<file under test>.cpp`. The checker blocks the layout
  and the `test_` prefix; confirm in review that the rest of the name matches
  the file under test.
- Keep default C++ changes compatible with the repository C++11 baseline.
- Ask for focused tests or explicit test rationale for feature changes, bug
  fixes, INPUT behavior changes, heterogeneous kernels, and core-module
  refactors.
- Treat CI governance findings as deterministic evidence and semantic review
  findings as advisory until maintainers approve them.

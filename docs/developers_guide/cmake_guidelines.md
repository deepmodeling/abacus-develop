# CMake Build Guidelines

This guide documents the architecture, conventions, and automated checks for
ABACUS's CMake build system. It complements
`docs/developers_guide/agent_governance.md` and `AGENTS.md`.

## Responsibility Boundaries

ABACUS separates build configuration from target construction across four
layers:

1. **Root `CMakeLists.txt`** — Project declaration, user-facing `option()` /
   `cmake_dependent_option()` declarations, C++ standard selection, and
   `find_package()` calls. Do not construct or link ABACUS targets here.
2. **`cmake/modules/Find<Pkg>.cmake` and `Setup<Pkg>.cmake`** — Package and
   backend discovery. Each module should export standard `IMPORTED` targets
   (e.g., `Pkg::Target` carrying `INTERFACE_INCLUDE_DIRECTORIES` and
   `INTERFACE_LINK_LIBRARIES`) instead of leaking raw library paths or global
   macros into `source/`.
3. **`source/CMakeLists.txt`** — The single place that assembles:
   - shared compile usage requirements (`abacus::compile_requirements`,
     configured via `cmake/CompilerConfiguration.cmake`);
   - global feature macro definitions (`_abacus_feature_definitions`, applied
     per target by `abacus_apply_build_options()`);
   - ordered external link closures (`abacus::feature_libs`,
     `abacus::linalg_libs`, `abacus::link_libs`);
   - the main executable and its ordered internal target linkage.
4. **Module `source/<module>/CMakeLists.txt` and `unittests/CMakeLists.txt`** —
   Own only their local targets (`OBJECT` libraries or `target_sources()`
   additions) and unit tests registered via `AddTest()`. They must not modify
   global compiler flags, global include paths, or discover external packages.

Remember that ABACUS still maintains two parallel build surfaces alongside the
main CMake build. When adding, moving, splitting, or deleting C++ source files,
keep all three in sync:

- the relevant C++ `CMakeLists.txt`;
- `source/Makefile.Objects`;
- `python/pyabacus/src/*/CMakeLists.txt` (when shared files are affected).

## Formatting

All CMake files (`CMakeLists.txt` and `*.cmake`) follow `cmake-format` with the
repository configuration in `cmake-format.yaml` (2-space indentation, markup
disabled, `AddTest` keyword metadata configured).

Run locally before committing:

```bash
pre-commit run cmake-format --files <changed-cmake-files>
```

The `Precommit` workflow enforces `cmake-format` as a blocking check on the PR
diff.

## Forbidden And Deprecated Patterns

The governance checker (`tools/03_code_analysis/agent_governance_check.py`)
**blocks** newly added lines that reintroduce legacy patterns:

- **Removed variables and macros** (main build tree: root `CMakeLists.txt`,
  `cmake/`, `source/`, `tests/`):
  - `${math_libs}` — removed; `AddTest()` links `abacus::linalg_libs`
    automatically and `source/CMakeLists.txt` owns the executable link closure.
  - `USE_LIBXC` — renamed to `__LIBXC`.
  - `ENABLE_RAPIDJSON` / `USE_RAPIDJSON` — removed; use `ENABLE_JSON`
    (`__JSON`, backed by `nlohmann_json`).
  - `ENABLE_LIBCOMM` — removed; LibComm is discovered with `ENABLE_LIBRI`.
  - `remove_definitions(...)` — ineffective against target-scoped compile
    definitions; if a legacy test target genuinely must drop a feature macro,
    use `abacus_disable_feature_definitions(...)` (discouraged for new tests).
- **Directory-scoped and package-discovery commands under `source/`**:
  - `add_definitions()`, `add_compile_definitions()`, `add_compile_options()`,
    `include_directories()`, `link_directories()`, `link_libraries()` — use
    `target_*()` commands, `cmake/CompilerConfiguration.cmake`, or
    `_abacus_feature_definitions` in `source/CMakeLists.txt`.
  - `find_package()`, `find_library()` — place dependency discovery in root
    `CMakeLists.txt` and `cmake/modules/Find<Pkg>.cmake` (or
    `cmake/Testing.cmake` for test-only packages).
- **Legacy test staging in unit-test CMake files under `source/`**:
  - `install(...)` and `file(COPY ...)` — use `abacus_copy_test_files(...)`
    from `cmake/Testing.cmake`, which prepares runtime scripts and data at
    configure time via `configure_file(... COPYONLY)` and tracks directory
    changes with `CONFIGURE_DEPENDS`.

## Standard Templates

### 1. Adding an Optional External Dependency

Follow this checklist so options, macros, linkage, and `--version` / build
summary stay consistent:

1. **Declare the option in root `CMakeLists.txt`**:
   - Use `ENABLE_<FEATURE>` for optional features (or `USE_<BACKEND>` for
     hardware/math backends).
   - If the option depends on another option (e.g. `ENABLE_MPI` or
     `ENABLE_LCAO`), use `cmake_dependent_option()` rather than manual `if()`
     overrides.
2. **Discover the package in root `CMakeLists.txt` + `cmake/modules/`**:
   - Prefer upstream `*Config.cmake` packages when available.
   - If a custom module is needed, add `cmake/modules/Find<Pkg>.cmake` and
     export an `IMPORTED` target `Pkg::Target`.
3. **Wire compile requirements and linkage in `source/CMakeLists.txt`**:
   - Add `$<$<BOOL:${ENABLE_<FEATURE>}>:__<FEATURE>>` to
     `_abacus_feature_definitions` (use the `__<FEATURE>` macro convention, not
     subdirectory-local `USE_<FEATURE>` definitions).
   - Add `$<$<BOOL:${ENABLE_<FEATURE>}>:Pkg::Target>` to `_abacus_feature_libs`
     (or `_abacus_linalg_libs` for numerical backends).
   - Add `Pkg::Target` to the `abacus_add_target_compile_requirements(...)` loop
     so object libraries and unit tests inherit its headers and compile flags.
4. **Update build summary**:
   - Record the option status and detected version/path in
     `cmake/CollectBuildInfoVars.cmake` (and `cmake/BuildInfo.cmake` if a new
     placeholder is printed).

### 2. Adding a Module or Submodule

Choose between two patterns:

- **Submodule that extends an existing library** (tightly coupled to parent, no
  separate link boundary needed):

```cmake
target_sources(parent_target PRIVATE foo.cpp bar.cpp)

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
```

- **Standalone `OBJECT` library** (distinct link boundary, optional feature, or
  isolated so lightweight unit tests can link it cleanly):

```cmake
add_library(my_module OBJECT foo.cpp bar.cpp)

if(ENABLE_COVERAGE)
  add_coverage(my_module)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
```

When adding a new `OBJECT` library:

- Link it in `source/CMakeLists.txt` in **consumer-before-provider** order
  (static linkers resolve symbols left-to-right).
- If callers span both PW and LCAO, keep the library unconditional rather than
  gating it inside `if(ENABLE_LCAO)`.

### 3. Adding Unit Tests

Per `AGENTS.md` rule 11, unit tests under `source/` live in a flat `unittests/`
directory with source files named `test_<file_under_test>.cpp`:

```cmake
abacus_copy_test_files(support)

AddTest(
  TARGET MODULE_FOO_bar
  LIBS my_module base parameter
  SOURCES test_bar.cpp)
```

Guidelines for unit-test `CMakeLists.txt` files:

- `AddTest()` automatically links `GTest::gtest_main`, `GTest::gmock_main`, and
  `abacus::linalg_libs` — do not pass `${math_libs}` or duplicate `GTest`
  targets.
- Prefer linking module `OBJECT` targets in `LIBS` rather than recompiling
  `../foo.cpp` or `../../bar.cpp` directly in `SOURCES`; direct cross-directory
  source lists break downstream tests whenever production files are split or
  moved.
- Avoid `abacus_disable_feature_definitions(...)` in new unit tests
  (**discouraged**). It only strips macros from targets defined in that test
  directory; linked production `OBJECT` libraries keep their original MPI/GPU
  compile definitions, so reusing them still requires compatible feature
  assumptions and initialization. Any test that genuinely needs those features
  disabled must link or compile objects built with matching settings. In legacy
  directories that still call `abacus_disable_feature_definitions(...)`, an
  individual `AddTest(...)` call may pass `KEEP_FEATURE_DEFINITIONS __MPI` to
  opt back in.
- For tests that require MPI, guard `add_subdirectory(unittests)` or the
  specific `AddTest()` target with `if(ENABLE_MPI)` so serial builds
  (`-DENABLE_MPI=OFF -DBUILD_TESTING=ON`) stay green.

## Governance Enforcement Summary

| Check | Executor | Severity | Action |
| --- | --- | --- | --- |
| `cmake-format` on changed CMake files | CI (`Precommit` workflow) + local hook | medium | block |
| New source / heterogeneous file CMake linkage | `agent_governance_check.py` (hook + CI) | medium | block |
| Deprecated CMake variables, macros, and `remove_definitions()` | `agent_governance_check.py` (hook + CI) | high | block |
| Directory-scoped commands and `find_package`/`find_library` under `source/` | `agent_governance_check.py` (hook + CI) | high | block |
| `install()` / `file(COPY ...)` in unit-test CMake files | `agent_governance_check.py` (hook + CI) | medium | block |
| New root options or `cmake/modules/` files without `CollectBuildInfoVars.cmake` sync | `agent_governance_check.py` (hook + CI) + AI review | medium | warn |
| Unit test `unittests/` directory layout and `test_` prefix | `agent_governance_check.py` (hook + CI) | high | block |

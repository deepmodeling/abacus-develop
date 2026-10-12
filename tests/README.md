# ABACUS tests

## Purpose

The purposes for this test directory are:

1. Cover most features of ABACUS.
2. Provide an autotest script to check if the version is correct.
   (Reference results are calculated by one core and saved in `result.ref`.
   You can change `NUMBEROFPROCESS` in `integrate/general_info` to test with
   multiple cores, and set the executable once via
   `export ABACUS_EXE=/path/to/abacus`.)

## Folders in this directory

| Folder | Description |
| --- | --- |
| `01_PW` | KSDFT calculations in PW basis with multiple k-point setting. |
| `02_NAO_Gamma` | KSDFT calculations in NAO basis with gamma-only k-point setting. |
| `03_NAO_multik` | KSDFT calculations in NAO basis with multiple k-point setting. |
| `04_FF` | Force fields, including Lennard-Jones potentials, Deep Potentials and Neuroevolution Potential. |
| `05_rtTDDFT` | Real-time TDDFT tests. |
| `06_SDFT` | Stochastic DFT tests. |
| `07_OFDFT` | Orbital-free DFT tests. |
| `08_RI` | Resolution-of-identity (RI) based tests, including hybrid functionals, LR-TDDFT, RPA and BSE. |
| `09_DeePKS` | DeePKS tests. |
| `10_others` | Other tests such as LCAO in pw. |
| `11_PW_GPU` | KSDFT calculations in PW basis with multiple k-point setting using GPU. |
| `12_NAO_Gamma_GPU` | KSDFT calculations in NAO basis with gamma-only k-point setting using GPU. |
| `13_NAO_multik_GPU` | KSDFT calculations in NAO basis with multiple k-point setting using GPU. |
| `15_rtTDDFT_GPU` | Real-time TDDFT tests using LCAO basis running on GPU. |
| `16_SDFT_GPU` | Integrate tests for stochastic DFT running on GPU. |
| `CMakeLists.txt` | CMake file for this directory. |
| `integrate` | Stores scripts for integrate tests. |
| `libxc` | Examples related to LibXC, will be refactored soon. |
| `performance` | Examples related to performance of ABACUS, will be refactored soon. |
| `PP_ORB` | Collection of all the used pseudopotentials and numerical atomic orbitals. |
| `README.md` | This file. |

## `CASES_CPU.txt` / `CASES_GPU.txt`

Each category directory (`01_PW`, `02_NAO_Gamma`, ...) contains its own
`CASES_CPU.txt` and/or `CASES_GPU.txt`. These are plain lists of the sub-case
names (the test directories) to run, one per line:

- `CASES_CPU.txt`: cases run on CPU builds.
- `CASES_GPU.txt`: cases run on GPU builds.

`Autotest.sh` reads the list from the file given by `-f <file>` (default
`CASES_CPU.txt` in the current working directory) and runs each named case.
Ctest invokes `Autotest.sh` per category directory with the appropriate `-f`
file, e.g. `tests/11_PW_GPU/CMakeLists.txt` uses `-f CASES_GPU.txt`.

When you add/rename a test case, add/remove its directory name in the relevant
`CASES_*.txt` of its category, otherwise it will not be picked up.

## How to run tests

1. Set the ABACUS executable. The recommended way is a single environment
   variable, which both `Autotest.sh` and `general_info` (via `run_check.sh`)
   read:

   ```bash
   export ABACUS_EXE=/path/to/abacus
   ```

   If `ABACUS_EXE` is not set, both fall back to the bare command name
   `abacus`, which the shell resolves by searching the directories listed in
   `$PATH` (run `which abacus` to see what it resolves to). This PATH fallback
   is what CI uses after `cmake --install`.

   Ways to point at a specific executable, highest priority first:

   a) Command-line flag (`Autotest.sh` only):

      ```bash
      ../integrate/Autotest.sh -a /path/to/abacus -n 4
      ```

      Overrides everything for that invocation. `Single_job.sh`/`run_check.sh`
      do not accept `-a`.
   b) Environment variable (works for `Autotest.sh` and `run_check.sh`):

      ```bash
      export ABACUS_EXE=/path/to/abacus
      ```

      Add it to `~/.bashrc` to make it persistent across sessions.
   c) Local override file (`Autotest.sh` only; not tracked by git):
      create `integrate/general_info.local` containing e.g.

      ```bash
      abacus=/path/to/abacus
      ```

      `Autotest.sh` sources this file at startup when it exists, so the path
      survives new shells while `git status` stays clean.
   d) Edit the files directly (hard-coded; do NOT commit a local absolute path):
      - `integrate/Autotest.sh`: set the `abacus` variable.
      - `integrate/general_info`: set `EXEC` (supports `$VAR`/`${VAR}`
        placeholders, expanded by `integrate/validation_tools/run_check.sh`
        from the environment).

   Effective priority: `-a` flag > `ABACUS_EXE` > `general_info.local` >
   hard-coded value in the files > `abacus` from PATH.

2. Enter each integrate test directory, and run this script for autotests:

   ```bash
   ../integrate/Autotest.sh
   ```

3. If you want to focus on No.xxx example, such as `101_PW_OU`:

   ```bash
   cd 101_PW_OU
   ../../integrate/Single_job.sh $parameter
   ```

   You can choose `$parameter` among `""` (empty), `debug` or `ref`.
   `ref`: generate `result.ref` file (the answer you need).

## How to validate an integrate test case

Each case directory contains a checked-in `result.ref`: a list of
"key value" lines (total energy, forces, stress, matrix-compare flags, ...).
Validation means running ABACUS in the case directory, collecting the same
lines from the fresh run into `result.out`, and comparing the two:

- The collectors live in `integrate/validation_tools/`. The entry point is
  `catch_properties.sh`; it sources `props_common.sh` (shared helpers and
  INPUT parsing) and the per-category modules `props_basic.sh`, `props_mat.sh`,
  `props_cube.sh`, `props_ml.sh`, `props_tddft.sh` and `props_deepks.sh`. Each
  module defines `run_<category>_props()` hooks; add a new kind of check by
  extending the matching module (or dropping in a new `props_<cat>.sh` and
  calling its hook from `catch_properties.sh`).
- `CompareFile.py` compares numerical/text output files at a given digit
  tolerance; `cube_tool.py` handles real-space `.cube` integration and
  wave-function fingerprints.
- `totaltimeref` is the wall time of the run. It is recorded for information
  only and is expected to differ between machines.

Ways to validate, from a whole category down to one case:

1. Whole category (what ctest/CI do), from the category directory
   (e.g. `tests/01_PW`):

   ```bash
   bash ../integrate/Autotest.sh -n 4
   ```

   `Autotest.sh` runs every case in the category's `CASES_CPU.txt`, invokes the collectors to
   build `result.out`, and compares it with `result.ref` using the thresholds
   at the top of `Autotest.sh` (`threshold` / `force_threshold` /
   `stress_threshold` / `descriptor_threshold`; a per-case `threshold` file
   can override them). Useful flags:

   - `-r <regex>`: run only cases whose name matches `<regex>`
   - `-f <file>`: take the case list from `<file>` (e.g. `CASES_GPU.txt`)
   - `-n <N>`: MPI ranks per case; `-o <N>` pins OpenMP threads
   - `-j <N>|auto`: number of cases to run concurrently
   - `-g`: generate `result.ref` instead of checking
   - `-a <exe>`: ABACUS executable for this invocation

   Example:

   ```bash
   bash ../integrate/Autotest.sh -r 035_PW_15_SO
   ```

2. One case by hand, from the case directory (e.g. `tests/01_PW/035_PW_15_SO`):

   ```bash
   bash ../../integrate/Single_job.sh
   ```

   It reads `integrate/general_info` (`NUMBEROFPROCESS`, `CHECKACCURACY`,
   `EXEC`), runs ABACUS, builds `result.out` through `run_check.sh` and prints
   every line that disagrees with `result.ref`. With the `ref` argument it
   regenerates `result.ref` instead; `debug` copies `run_check.sh` and the
   whole `validation_tools/` directory into the case directory and collects
   with that copy, so local edits to the collectors take effect.

3. Collector only (ABACUS output already exists in `OUT.autotest/`):

   ```bash
   bash ../../integrate/validation_tools/catch_properties.sh result.out
   diff result.ref result.out
   ```

   A clean diff (apart from the `totaltimeref` line) means the case still
   validates. This is the fastest check when you only changed the validation
   scripts themselves.

## How to modify an existing test case

When merging or renaming a test case:

1. Move/rename the directory: `git mv <old_dir> <new_dir>`
2. Update the INPUT file with new parameters (keep the original style)
3. Update the case README to describe the new test purpose
4. Update `CASES_CPU.txt` and `CASES_GPU.txt`: remove old name, add new name
5. Remove the old directory: `git rm -rf <old_dir>`
6. Regenerate `result.ref`: `OMP_NUM_THREADS=1 ../../integrate/Single_job.sh ref`
7. Commit the changes

For smearing tests: insulators (e.g., H2) can accept any smearing method
without changing the result. Metals require smearing for convergence.

Marking unimplemented features:
- If a feature is not yet implemented (e.g., a certain code path, GPU backend,
  or functional lacks stress/force support), add the corresponding INPUT
  parameter but comment it out with an explanation, e.g.:

  ```
  #cal_stress        1 # implicit solvation stress is not implemented yet
  ```

  This documents the intended coverage, makes the gap visible to readers, and
  lets the test be enabled later by simply uncommenting once the feature lands.

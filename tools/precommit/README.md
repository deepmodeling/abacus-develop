# ABACUS Precommit

This is a CP2K-style precommit system adapted for ABACUS.  The local driver is
`tools/precommit/precommit.py`; the formatting tools are provided by a small
Flask/gunicorn precommit server.

The design intentionally follows CP2K's usage model:

- one local driver script;
- `obj/precommit/cache.json` mtime cache;
- file-level processing;
- local backup before each tool modifies a file;
- default check-only mode;
- `--allow-modifications` to keep formatter changes;
- remote or local Docker server for external tools.

## Install Git hook

From the ABACUS repository root:

```bash
ln -fs ../../tools/precommit/precommit.py .git/hooks/pre-commit
```

## Run checks

```bash
./tools/precommit/precommit.py
```

The default mode does not modify source files.  If a formatter would change a
file, the driver restores the original content and prints a unified diff.

## Apply modifications

```bash
./tools/precommit/precommit.py --allow-modifications
```

Short option:

```bash
./tools/precommit/precommit.py -m
```

## Ignore cache

```bash
./tools/precommit/precommit.py --no-cache
```

Short option:

```bash
./tools/precommit/precommit.py -a
```

## Process selected files

```bash
./tools/precommit/precommit.py source/module/foo.cpp source/module/foo.h
```

## Local server

Until an ABACUS precommit server is deployed, run the server locally with Docker:

```bash
cd tools/precommit
./start_local_server.sh
```

Then, in another terminal from the repository root:

```bash
export ABACUS_PRECOMMIT_SERVER="http://127.0.0.1:8080"
./tools/precommit/precommit.py
```

## Tools

The first ABACUS version keeps CP2K's lightweight file-level model and maps it to
ABACUS' C++-first source tree:

- C/C++/CUDA/HIP/OpenCL: `clang-format`
- Python: `ast.parse` + `black`
- Shell: `shfmt` + `shellcheck`
- Markdown: `mdformat --wrap=100`
- CMake: `cmake-format -i`
- Makefile: local `format_makefile.py`
- all files: local `check_file_properties.py`

`clang-tidy` and `compile_commands.json` are intentionally not part of this
CP2K-style first version.  They are build-level static-analysis concerns, while
this precommit driver is a lightweight file-level formatting and convention gate.

## Notes for ABACUS maintainers

The only non-CP2K structural change is that C++ is first-class.  CP2K rejects
most C++ files; this ABACUS version formats `.c`, `.cc`, `.cpp`, `.cxx`, `.h`,
`.hh`, `.hpp`, `.hxx`, `.cu`, `.cuh`, `.hip`, and `.cl` files.

`check_file_properties.py` contains the ABACUS-specific convention checks.  The
banner/license policy should be tightened once the exact ABACUS source header is
agreed upon.

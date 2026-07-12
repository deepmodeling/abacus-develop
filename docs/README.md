# ABACUS Documentation

These files are the source for the ABACUS manual published on Read the Docs.

Read the Docs builds a minimal serial ABACUS binary, regenerates
`docs/advanced/input_files/input-main.md` from that binary during the
documentation build, and emits a warning if no ABACUS executable is available.

## Build the Manual Locally

Run the following commands from the repository root.

1. Create a Python virtual environment with `uv`:

   ```bash
   uv venv --python python3 .venv
   ```

1. Install the documentation requirements:

   ```bash
   uv pip install --python .venv/bin/python -r docs/requirements.txt sphinx ninja
   ```

1. Build the reduced ABACUS executable used for INPUT parameter metadata:

   ```bash
   cmake -S . -B build-rtd-docs -G Ninja \
     -DCMAKE_BUILD_TYPE=Release \
     -DENABLE_MPI=OFF \
     -DENABLE_LCAO=OFF \
     -DUSE_OPENMP=OFF \
     -DUSE_ELPA=OFF \
     -DUSE_CUDA=OFF \
     -DUSE_ROCM=OFF \
     -DBUILD_TESTING=OFF \
     -DENABLE_LIBXC=OFF \
     -DENABLE_LIBRI=OFF \
     -DENABLE_DFTD4=OFF \
     -DENABLE_RAPIDJSON=OFF \
     -DENABLE_MLALGO=OFF \
     -DENABLE_FLOAT_FFTW=OFF \
     -DENABLE_CNPY=OFF \
     -DCOMMIT_INFO=OFF \
     -DGIT_SUBMODULE=OFF \
     -DMKLROOT=OFF
   cmake --build build-rtd-docs --target abacus_pw_ser --parallel 2
   ```

1. Build the HTML manual:

   ```bash
   ABACUS_BINARY=./build-rtd-docs/source/abacus_pw_ser \
     .venv/bin/sphinx-build -b html docs build-docs/html
   ```

   If `ABACUS_BINARY` is not set, `docs/conf.py` looks for the Read the Docs
   build binary and then for `abacus` or `abacus_pw_ser` on `PATH`.

1. Open `build-docs/html/index.html` in a browser.

If no ABACUS executable is available, the Sphinx build still continues with the
checked-in `docs/advanced/input_files/input-main.md` and prints a warning that
the INPUT parameter reference may not be up to date.

## Regenerate Only the INPUT Reference

To refresh only the generated INPUT reference after building ABACUS:

```bash
./build-rtd-docs/source/abacus_pw_ser --generate-parameters-yaml \
  | .venv/bin/python docs/generate_input_main.py - \
      --output docs/advanced/input_files/input-main.md
```

## INPUT Parameter Reference

The INPUT parameter reference is generated from metadata registered in
`source/source_io/module_parameter/read_input_item_*.cpp`.

- `abacus --generate-parameters-yaml` produces a transient YAML stream; it is
  not stored in the repository.
- `docs/advanced/input_files/input-main.md` is generated from that metadata.
- `docs/conf.py` refreshes `input-main.md` automatically when it can find an
  ABACUS executable.

Keep `docs/advanced/input_files/input-main.md` in PRs that change INPUT
metadata. The YAML dump is intentionally transient so there is only one
generated reference file to maintain.

## Optional Read the Docs Container Test

To smoke-test the Read the Docs environment locally, use the same image as the
hosted build:

```bash
docker run --rm -v "$PWD:/project" -w /project \
  readthedocs/build:ubuntu-22.04-2024.01.29 \
  /bin/bash -lc 'python -m pip install -r docs/requirements.txt sphinx ninja &&
  cmake -S . -B build-rtd-docs -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_MPI=OFF -DENABLE_LCAO=OFF -DUSE_OPENMP=OFF -DUSE_ELPA=OFF \
    -DUSE_CUDA=OFF -DUSE_ROCM=OFF -DBUILD_TESTING=OFF \
    -DENABLE_LIBXC=OFF -DENABLE_LIBRI=OFF -DENABLE_DFTD4=OFF \
    -DENABLE_RAPIDJSON=OFF -DENABLE_MLALGO=OFF -DENABLE_FLOAT_FFTW=OFF \
    -DENABLE_CNPY=OFF -DCOMMIT_INFO=OFF -DGIT_SUBMODULE=OFF -DMKLROOT=OFF &&
  cmake --build build-rtd-docs --target abacus_pw_ser --parallel 2 &&
  export ABACUS_BINARY=./build-rtd-docs/source/abacus_pw_ser &&
  sphinx-build -b html docs /tmp/abacus-docs-html'
```

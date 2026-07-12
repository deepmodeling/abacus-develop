# ABACUS Documentation

These files are the source for the ABACUS manual published on Read the Docs.

Read the Docs builds a minimal serial ABACUS binary, regenerates
`docs/parameters.yaml` from that binary, and lets Sphinx regenerate
`docs/advanced/input_files/input-main.md` during the documentation build.

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

1. Regenerate the INPUT parameter files:

   ```bash
   ./build-rtd-docs/source/abacus_pw_ser --generate-parameters-yaml > docs/parameters.yaml
   .venv/bin/python docs/generate_input_main.py docs/parameters.yaml \
     --output docs/advanced/input_files/input-main.md
   ```

1. Build the HTML manual:

   ```bash
   .venv/bin/sphinx-build -b html docs build-docs/html
   ```

1. Open `build-docs/html/index.html` in a browser.

## INPUT Parameter Reference

The INPUT parameter reference is generated from metadata registered in
`source/source_io/module_parameter/read_input_item_*.cpp`.

- `docs/parameters.yaml` is an intermediate YAML dump produced by
  `abacus --generate-parameters-yaml`.
- `docs/advanced/input_files/input-main.md` is generated from that YAML file.
- `docs/conf.py` refreshes `input-main.md` automatically when
  `docs/parameters.yaml` exists.

The Read the Docs build no longer depends on a committed `docs/parameters.yaml`
because `.readthedocs.yaml` regenerates it before Sphinx runs. A local Sphinx
build also succeeds without `docs/parameters.yaml` when the checked-in
`input-main.md` is already present, but a fresh INPUT reference still requires
an ABACUS binary.

For now, keep both generated files in PRs that change INPUT metadata, following
the repository contribution guide. Dropping `docs/parameters.yaml` from version
control is technically possible after the project also updates the contributor
guide and PR template to make the generated-file policy explicit.

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
  ./build-rtd-docs/source/abacus_pw_ser --generate-parameters-yaml > docs/parameters.yaml &&
  sphinx-build -b html docs /tmp/abacus-docs-html'
```

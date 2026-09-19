# ABACUS Documentation

These files are the source for the ABACUS manual published on Read the Docs.

Read the Docs builds a minimal serial ABACUS binary, regenerates
`docs/advanced/input_files/input-main.md` from that binary during the
documentation build, and fails the build if the INPUT reference cannot be
refreshed.

## Build the Manual Locally

Run the following commands from the repository root.

1. Create a Python virtual environment:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

1. Install the Python documentation requirements:

   ```bash
   pip3 install -r docs/requirements.txt
   ```

The reduced ABACUS build also requires a C++ compiler, a Fortran compiler,
CMake, Ninja, FFTW3, BLAS, and LAPACK. These are system/build dependencies,
not Python documentation dependencies.

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

1. (Optional) Generate the INPUT reference manually for inspection:

   The ABACUS executable emits the parameter metadata as a transient YAML
   stream. Pipe that stream directly to the Markdown generator; neither
   `docs/parameters.yaml` nor `docs/advanced/input_files/input-main.md` is
   maintained in the repository.

   ```bash
   ./build-rtd-docs/abacus_pw_ser --generate-parameters-yaml \
     | .venv/bin/python docs/generate_input_main.py - \
         --output /tmp/input-main.md
   ```

1. Build the HTML manual:

   ```bash
   ABACUS_BINARY=./build-rtd-docs/abacus_pw_ser \
     .venv/bin/sphinx-build -b html docs build-docs/html
   ```

   If `ABACUS_BINARY` is not set, `docs/conf.py` looks for the Read the Docs
   build binary and then for `abacus` or `abacus_pw_ser` on `PATH`.

1. Open `build-docs/html/index.html` in a browser.

Every documentation build requires an ABACUS executable. If no executable is
available, or if parameter generation fails, Sphinx stops with an error rather
than publishing incomplete or stale INPUT documentation.

## Regenerate Only the INPUT Reference

To refresh only the generated INPUT reference after building ABACUS, pipe the
transient YAML stream directly into the Markdown generator:

```bash
./build-rtd-docs/abacus_pw_ser --generate-parameters-yaml \
  | .venv/bin/python docs/generate_input_main.py - \
      --output /tmp/input-main.md
```

## INPUT Parameter Reference

The INPUT parameter reference is generated from metadata registered in
`source/source_io/module_parameter/read_input_item_*.cpp`.

- `abacus --generate-parameters-yaml` produces a transient YAML stream; it is
  not stored in the repository.
- `docs/advanced/input_files/input-main.md` is generated from that metadata
  during every Sphinx build and is not stored in the repository.
- `docs/conf.py` requires an ABACUS executable and fails the build if the INPUT
  reference cannot be generated.

PRs that change INPUT metadata only update the C++ `Input_Item` registrations.
The YAML stream and Markdown reference are both transient build artifacts.

## Optional Read the Docs Container Test

To smoke-test the Read the Docs environment locally, use the same image as the
hosted build:

```bash
docker run --rm --user root -v "$PWD:/project" -w /project \
  readthedocs/build:ubuntu-22.04-2024.01.29 \
  /bin/bash -lc 'apt-get update &&
  apt-get install -y build-essential cmake gfortran git libfftw3-dev \
    liblapack-dev libopenblas-dev ninja-build pkg-config python3-pip &&
  python3 -m pip install -r docs/requirements.txt &&
  cmake -S . -B /tmp/build-rtd-docs -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_MPI=OFF -DENABLE_LCAO=OFF -DUSE_OPENMP=OFF -DUSE_ELPA=OFF \
    -DUSE_CUDA=OFF -DUSE_ROCM=OFF -DBUILD_TESTING=OFF \
    -DENABLE_LIBXC=OFF -DENABLE_LIBRI=OFF -DENABLE_DFTD4=OFF \
    -DENABLE_RAPIDJSON=OFF -DENABLE_MLALGO=OFF -DENABLE_FLOAT_FFTW=OFF \
    -DENABLE_CNPY=OFF -DCOMMIT_INFO=OFF -DGIT_SUBMODULE=OFF -DMKLROOT=OFF &&
  cmake --build /tmp/build-rtd-docs --target abacus_pw_ser --parallel 2 &&
  export READTHEDOCS=True ABACUS_BINARY=/tmp/build-rtd-docs/abacus_pw_ser &&
  sphinx-build -b html docs /tmp/abacus-docs-html'
```

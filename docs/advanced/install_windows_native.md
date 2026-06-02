# Native Windows Build (experimental)

> **Status:** work in progress. This documents the *native* Windows port of
> ABACUS — a real Windows executable compiled with a Windows toolchain, as
> opposed to the [WSL2 one-click installer](./windows_installer.md), which runs
> the Linux binary inside WSL2 and remains the recommended way to **run**
> full-featured ABACUS on Windows.
>
> The port is staged:
> 1. **Phase 1 — serial, plane-wave (PW) only** ← *current target*
> 2. Phase 2 — serial, add LCAO
> 3. Phase 3 — MPI parallel
>
> Phases 1–3 deliberately exclude ELPA, PEXSI and hybrid functionals (LibRI),
> as well as GPU/DSP backends.

## Toolchain: MinGW-w64 GCC

The native build targets **MinGW-w64 GCC**, not MSVC. Reasons:

- MinGW ships the POSIX headers ABACUS relies on (`unistd.h`, `fcntl.h`,
  `sys/stat.h`, `dirent.h`, `access`, `open/read/write/close`, ...), so most
  I/O code compiles unchanged.
- The codebase has ~hundreds of GCC `__attribute__`/builtin usages (largely in
  vendored container code and CUDA kernels); GCC accepts them as-is, whereas
  MSVC would reject many.
- It pairs cleanly with OpenBLAS + FFTW3, which have good native Windows builds.

MSVC and Intel oneAPI (`icx`) remain possible future targets but are not the
Phase 1 path.

## Prerequisites

The simplest consistent way to get the compiler **and** the math libraries is
[MSYS2](https://www.msys2.org/):

```bash
# in an MSYS2 shell
pacman -S mingw-w64-x86_64-gcc \
          mingw-w64-x86_64-cmake \
          mingw-w64-x86_64-ninja \
          mingw-w64-x86_64-openblas \
          mingw-w64-x86_64-fftw
```

This provides, under `C:\msys64\mingw64`:

| Need            | Package                       |
|-----------------|-------------------------------|
| C++17 compiler  | `mingw-w64-x86_64-gcc`        |
| Build driver    | `mingw-w64-x86_64-cmake` + `ninja` |
| BLAS + LAPACK   | `mingw-w64-x86_64-openblas`   |
| FFTW3 (double)  | `mingw-w64-x86_64-fftw`       |

ScaLAPACK, ELPA, PEXSI and MPI are **not** required for Phase 1.

## Building (Phase 1: serial PW)

From a shell where the MinGW toolchain is on `PATH` (the "MSYS2 MinGW 64-bit"
shell, or PowerShell with `C:\msys64\mingw64\bin` on `PATH`):

```powershell
# PowerShell helper (this repo): tools/windows/build-native-serial.ps1
./tools/windows/build-native-serial.ps1 -PrefixPath "C:\msys64\mingw64"
```

Or invoke CMake directly:

```powershell
cmake -S . -B build_win_serial_pw -G Ninja `
      -DCMAKE_BUILD_TYPE=Release `
      -DENABLE_MPI=OFF -DENABLE_LCAO=OFF -DUSE_OPENMP=OFF `
      -DUSE_ELPA=OFF -DENABLE_PEXSI=OFF -DENABLE_LIBRI=OFF -DENABLE_MLALGO=OFF `
      -DUSE_CUDA=OFF -DBUILD_TESTING=OFF `
      -DCMAKE_PREFIX_PATH="C:\msys64\mingw64"
cmake --build build_win_serial_pw --parallel
```

The resulting executable is `abacus_pw_ser.exe` in the build directory.

### Validate against a Linux baseline

Run a small PW SCF case (e.g. `examples/02_scf/...`) and compare the total
energy / forces with a Linux serial build of the same commit. They should agree
to roughly machine precision (~1e-8 Ry).

## What changed in the source for the port

Phase 1 keeps the Linux build byte-for-byte identical; all changes are guarded
or platform-neutral:

- **`source/source_base/fs_compat.h`** (new): a portable `ModuleBase::make_directory()`
  wrapping `_mkdir` (Windows) / `mkdir(path, 0755)` (POSIX), since the Windows
  CRT `mkdir` takes no permission-mode argument.
- **`source/source_base/global_file.cpp`**, **`global_function.cpp`**: use the
  helper above instead of calling `mkdir(path, 0755)` directly.
- **`CMakeLists.txt`**:
  - `find_package(ScaLAPACK REQUIRED)` is now gated on `ENABLE_MPI` (a serial
    build must not require a distributed-memory library).
  - On Windows, defines `_USE_MATH_DEFINES`, `NOMINMAX`, `_CRT_SECURE_NO_WARNINGS`.
  - The default `-O3 -g` flags and the `-lm` link are skipped for MSVC.
  - The post-install `abacus` symlink step is skipped on Windows.

## Known limitations / not yet ported

- LCAO, MPI, ELPA, PEXSI, hybrid functionals (LibRI/LibComm), DeePKS/ML-KEDF,
  GPU (CUDA/ROCm), DSP — all disabled for Phase 1.
- Expect additional small portability fixes to surface during compilation;
  they are tracked as part of the staged port.

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

Install [MSYS2](https://www.msys2.org/). Everything else (compiler, math
libraries) is installed by the toolchain script below — there is no separate
Windows build script; the native Windows build is just another **toolchain
variant**, alongside `gnu`, `intel`, `gcc-mkl`, …

## Building (Phase 1: serial PW) — via the toolchain

Open the **"MSYS2 MinGW 64-bit"** shell and run the two toolchain scripts, the
same two-step flow as the Linux variants (`toolchain_gnu.sh` →
`build_abacus_gnu.sh`):

```bash
cd toolchain
./toolchain_windows.sh      # pacman-installs gcc/gfortran/openblas/fftw/cmake/ninja
./build_abacus_windows.sh   # configures + builds the serial PW binary
```

`toolchain_windows.sh` is the Windows counterpart of `toolchain_gnu.sh`: on
Linux the dependencies are built from source, while on MSYS2 they come from
`pacman` (under `/mingw64`). `build_abacus_windows.sh` then configures the
serial plane-wave build (`ENABLE_MPI=OFF`, `ENABLE_LCAO=OFF`, OpenBLAS + FFTW)
and builds `build_abacus_windows/abacus_pw_ser.exe`.

A few non-default options the build script sets, and why:
- `BLA_VENDOR=OpenBLAS` — OpenBLAS supplies both BLAS and LAPACK in one library.
- `ENABLE_FLOAT_FFTW=ON` — compiles `fft_cpu_float.cpp` so the `FFT_CPU<float>`
  vtable is fully defined (see *source changes* below); needs `libfftw3f`.
- `CMAKE_CXX_FLAGS="-include cstdint -include cstring -include algorithm"` —
  MSYS2 ships a very new GCC whose libstdc++ dropped several transitive
  standard-header includes; force-including the common ones lets the existing
  sources build unchanged. (Not Windows-specific — tied to GCC ≥ 15. A cleaner
  long-term fix is to add the missing `#include`s per file.)

To run the binary, `source toolchain/abacus_env.sh` (written by the build
script); it puts the binary and the MinGW runtime DLLs (libstdc++, libgcc,
libgfortran, libopenblas, libfftw3) on `PATH`.

## Testing — the existing `01_PW` suite, serial mode

There is **no separate Windows test script and no separate case list**. The PW
test suite is `tests/01_PW`, driven by the standard harness exactly as in CI
(`tests/01_PW/CMakeLists.txt` runs `Autotest.sh` from that directory, which
reads `tests/01_PW/CASES_CPU.txt`). The only addition is a **serial mode** in
`Autotest.sh`: `-n 0` runs the ABACUS binary directly with no MPI launcher.

From the MSYS2 MinGW 64-bit shell:

```bash
cd tests/01_PW
bash ../integrate/Autotest.sh \
     -a "$(pwd)/../../build_abacus_windows/abacus_pw_ser.exe" -n 0
```

This runs the whole `01_PW` suite and compares every case against its
`result.ref` with `tools/catch_properties.sh`, identical to the Linux/MPI run
apart from `-n 0`. (`bc`, used by `catch_properties.sh`, is installed by
`toolchain_windows.sh`.)

Cases requiring features outside the Phase 1 serial-PW build (multi-process
`kpar`, the ScaLAPACK solver, or LibXC functionals) are expected to report
warnings/failures; that is a property of the reduced build, not of the port.

## What changed in the source for the port

Phase 1 keeps the Linux build byte-for-byte identical; all changes are guarded
or platform-neutral:

- **`source/source_base/fs_compat.h`** (new): a portable `ModuleBase::make_directory()`
  wrapping `_mkdir` (Windows) / `mkdir(path, 0755)` (POSIX), since the Windows
  CRT `mkdir` takes no permission-mode argument.
- **`source/source_base/global_file.cpp`**, **`global_function.cpp`**: use the
  helper above instead of calling `mkdir(path, 0755)` directly.
- **`cmake/FindBlas.cmake`**, **`cmake/FindLapack.cmake`**: these wrappers delegate
  to CMake's builtin `FindBLAS`/`FindLAPACK`. On the case-insensitive Windows
  filesystem `FindBlas.cmake` and `FindBLAS.cmake` are the same file, so the
  delegating `find_package(BLAS)`/`find_package(LAPACK)` recursed into the
  wrapper forever ("maximum nesting depth exceeded"). Fixed by temporarily
  dropping our module dir from `CMAKE_MODULE_PATH` around the builtin call.
- **`source/source_base/module_fft/fft_base.h`**, **`fft_cpu.h`**: removed
  `__attribute__((weak))` from the FFT virtual functions. The weak-without-
  definition pattern relied on the ELF linker resolving unbound weak symbols to
  null; on Windows/PE (MinGW) it produced **null vtable slots**, so the first
  FFT dispatch (`FFT_Bundle::setupFFT`) jumped to address 0 and crashed. The
  base virtuals now have trivial default bodies; the float overrides are made
  concrete by building with `ENABLE_FLOAT_FFTW=ON`.
- **`source/source_io/module_parameter/input_conv.h`**: replaced the POSIX
  `<regex.h>` (`regcomp`/`regexec`) expression parser with portable C++
  `<regex>` (`std::regex`). MinGW has no `<regex.h>`.
- **`source/source_base/module_container/base/core/cpu_allocator.cpp`**: replaced
  `posix_memalign` (no Windows CRT equivalent) with `_aligned_malloc`/
  `_aligned_free` on Windows, used consistently across both `allocate` overloads
  and `free`.
- **`source/source_io/module_restart/restart.cpp`**: the POSIX owner-permission
  macros `S_IRUSR`/`S_IWUSR` are undefined in the Windows CRT; mapped them to
  `_S_IREAD`/`_S_IWRITE` and include `<io.h>` for the low-level `open/read/
  write/close`.
- **`source/source_psi/psi_initializer.cpp`**: fixed the seeded (`pw_seed>0`)
  random-wavefunction path in **serial** builds. The per-stick random data was
  only gathered into the working arrays via `stick_to_pool()` under `#ifdef
  __MPI`, so without MPI the wavefunctions stayed all-zero and tripped
  Gram-Schmidt (`psi_norm <= 0.0`). Added the serial direct-copy counterpart.
  (Pre-existing serial bug, not Windows-specific — CI only runs under MPI.)
- **`CMakeLists.txt`**:
  - `find_package(ScaLAPACK REQUIRED)` is now gated on `ENABLE_MPI` (a serial
    build must not require a distributed-memory library).
  - On Windows, defines `_USE_MATH_DEFINES`, `NOMINMAX`, `_CRT_SECURE_NO_WARNINGS`.
  - The default `-O3 -g` flags and the `-lm` link are skipped for MSVC.
  - The post-install `abacus` symlink step is skipped on Windows.
- **`tests/integrate/Autotest.sh`**: added a serial mode (`-n 0`) that runs the
  binary without an MPI launcher, so serial builds reuse the standard harness.
- **`toolchain/toolchain_windows.sh`**, **`toolchain/build_abacus_windows.sh`**
  (new): the native-Windows toolchain variant (MSYS2/MinGW-w64), mirroring the
  `gnu`/`intel`/`gcc-mkl` variants.

## Known limitations / not yet ported

- LCAO, MPI, ELPA, PEXSI, hybrid functionals (LibRI/LibComm), DeePKS/ML-KEDF,
  GPU (CUDA/ROCm), DSP — all disabled for Phase 1.
- Expect additional small portability fixes to surface during compilation;
  they are tracked as part of the staged port.

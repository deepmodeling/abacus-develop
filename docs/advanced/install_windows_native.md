# Native Windows Build (experimental)

> **Status:** work in progress. This documents the *native* Windows port of
> ABACUS — a real Windows executable compiled with a Windows toolchain, as
> opposed to the [WSL2 one-click installer](./windows_installer.md), which runs
> the Linux binary inside WSL2 and remains the recommended way to **run**
> full-featured ABACUS on Windows.
>
> The port was staged, all three phases are now working:
> 1. **Phase 1 — serial, plane-wave (PW)** ✓
> 2. **Phase 2 — add LCAO** ✓ (serial LCAO works for multi-k; the gamma-only
>    serial path has a known bug — see *Known limitations*)
> 3. **Phase 3 — MPI parallel (MS-MPI + ScaLAPACK)** ✓ — the default build
>
> It deliberately excludes ELPA, PEXSI, hybrid functionals (LibRI/LibComm),
> DeePKS/ML-KEDF, LibXC, and GPU/DSP backends — these have no reliable native
> Windows build yet and remain ordinary feature switches.
>
> Validated against the `01_PW`, `02_NAO_Gamma`, and `03_NAO_multik` test
> suites (via the standard harness): under MPI all three pass within the
> expected cross-platform error range; the residual warnings are float noise
> at the harness's strict absolute thresholds or excluded features (e.g. SCAN
> meta-GGA needs LibXC).

## Toolchain: MinGW-w64 GCC

The native build targets **MinGW-w64 GCC**, not MSVC. Reasons:

- MinGW ships the POSIX headers ABACUS relies on (`unistd.h`, `fcntl.h`,
  `sys/stat.h`, `dirent.h`, `access`, `open/read/write/close`, ...), so most
  I/O code compiles unchanged.
- The codebase has ~hundreds of GCC `__attribute__`/builtin usages (largely in
  vendored container code and CUDA kernels); GCC accepts them as-is, whereas
  MSVC would reject many.
- It pairs cleanly with OpenBLAS + FFTW3, which have good native Windows builds.

MSVC and Intel oneAPI (`icx`) remain possible future targets.

## Prerequisites

Install [MSYS2](https://www.msys2.org/). Everything else (compiler, math
libraries) is installed by the toolchain script below — there is no separate
Windows build script; the native Windows build is just another **toolchain
variant**, alongside `gnu`, `intel`, `gcc-mkl`, …

## Building — via the toolchain

Open the **"MSYS2 MinGW 64-bit"** shell and run the two toolchain scripts, the
same two-step flow as the Linux variants (`toolchain_gnu.sh` →
`build_abacus_gnu.sh`):

```bash
cd toolchain
./toolchain_windows.sh      # pacman-installs gcc/gfortran/openblas/fftw/cmake/ninja/
                            #                 cereal/msmpi/scalapack/bc
./build_abacus_windows.sh   # configures + builds the MPI + LCAO binary
```

`toolchain_windows.sh` is the Windows counterpart of `toolchain_gnu.sh`: on
Linux the dependencies are built from source, while on MSYS2 they come from
`pacman` (under `/mingw64`). `build_abacus_windows.sh` then builds the **MPI +
LCAO** configuration by default (`abacus_basic_para.exe`, OpenBLAS + FFTW +
ScaLAPACK). Pick a lighter configuration with environment toggles:

```bash
ENABLE_MPI=OFF ./build_abacus_windows.sh                 # serial LCAO+PW
ENABLE_MPI=OFF ENABLE_LCAO=OFF ./build_abacus_windows.sh # serial PW only
```

The MPI build needs the **MS-MPI runtime** (`msmpi.dll`, `mpiexec`) installed
system-wide — a separate Microsoft redistributable — in addition to the MinGW
`msmpi` package that `toolchain_windows.sh` installs for building.

A few non-default options the build script sets, and why:
- `BLA_VENDOR=OpenBLAS` — OpenBLAS supplies both BLAS and LAPACK in one library.
- `ENABLE_FLOAT_FFTW=ON` — compiles `fft_cpu_float.cpp` so the `FFT_CPU<float>`
  vtable is fully defined (see *source changes* below); needs `libfftw3f`.
- `CMAKE_CXX_FLAGS="-include cstdint -include cstring -include algorithm"` —
  MSYS2 ships a very new GCC whose libstdc++ dropped several transitive
  standard-header includes; force-including the common ones lets the existing
  sources build unchanged. (Not Windows-specific — tied to GCC ≥ 15. A cleaner
  long-term fix is to add the missing `#include`s per file.)

To run it, `source toolchain/abacus_env.sh` and then call `abacus` directly —
exactly like the Linux toolchain:

```bash
source toolchain/abacus_env.sh
abacus --version
```

`abacus_env.sh` puts the binary directory on `PATH` and sets `OMP_NUM_THREADS=1`
(plus `OPENBLAS_NUM_THREADS=1`). The build step also **bundles the dependent
DLLs** (libstdc++, libgcc, libgfortran, libquadmath, libgomp, libwinpthread,
libopenblas, libfftw3, libfftw3f, libscalapack) next to `abacus.exe`; Windows
searches the application directory before `PATH`, so the binary is
self-contained. Because native Windows symlinks need elevation, the build copies
the configured binary to `abacus.exe` (instead of the Linux `abacus` symlink),
so a bare `abacus` resolves in the MSYS2 shell and in cmd/PowerShell. Run in
parallel with MS-MPI:

```bash
mpiexec -n 4 abacus
```

`OMP_NUM_THREADS=1` matters under MPI. MSYS2's OpenBLAS is *OpenMP*-threaded
(it links `libgomp`), so `OMP_NUM_THREADS` — **not** the commonly cited
`OPENBLAS_NUM_THREADS` — is what actually caps its threads. Without it, each MPI
rank spawns a multithreaded BLAS, the ranks oversubscribe the cores, and
OpenBLAS's buffer allocator dies ("Memory allocation still failed after 10
retries"). ABACUS itself is built `USE_OPENMP=OFF`, so pinning the BLAS to one
thread per rank costs nothing — parallelism comes from MPI.

## Testing — the existing harness

There is **no separate Windows test script and no separate case list**. The
suites `tests/01_PW`, `tests/02_NAO_Gamma`, `tests/03_NAO_multik` are driven by
the standard harness exactly as in CI (`tests/<suite>/CMakeLists.txt` runs
`Autotest.sh` from that directory, which reads its `CASES_CPU.txt`).

**Parallel (recommended — matches the MPI references):** just run the harness
normally. MS-MPI's launcher is `mpiexec`, not the `mpirun` the harness invokes,
so the build drops a small **`mpirun` shim** next to the binary (on `PATH` via
`abacus_env.sh`) that forwards to `mpiexec` and pins `OMP_NUM_THREADS=1`. With
that in place the default invocation works unchanged:

```bash
source toolchain/abacus_env.sh
cd tests/02_NAO_Gamma
bash ../integrate/Autotest.sh -a abacus      # default np=4, via the mpirun shim
```

**Serial:** `Autotest.sh` also gained a serial mode — `-n 0` runs the binary
directly with no MPI launcher — for a serial build (PW and multi-k LCAO only;
gamma-only LCAO must use the MPI build, see below):

```bash
cd tests/01_PW
bash ../integrate/Autotest.sh -a abacus -n 0
```

Either way the harness compares every case against its `result.ref` with
`tools/catch_properties.sh`. (`bc`, used by that script, is installed by
`toolchain_windows.sh`.)

Expected residual differences (not bugs): cross-platform/cross-BLAS floating
point that just exceeds the harness's strict absolute thresholds (energies
still match to ~1e-7 eV); gauge-dependent outputs (raw wavefunction values,
Wannier `.amn`); a few file comparisons at ~1e-6; the init-sensitive
`078_PW_S2_elec_add` (see the `pw_seed` note); and excluded features
(SCAN/meta-GGA needs LibXC, DFT+U requires MPI, etc.).

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
- **`source/source_pw/module_pwdft/structure_factor.cpp`**: same family of bug
  in `bspline_sf` (`nbspline>0`) — the real-space plane scatter (`zpiece_to_all`)
  was MPI-only, leaving the structure factor uninitialized in serial → wrong
  total energy/force/stress. Added the serial direct-fill.
- **`source/source_io/module_output/binstream.cpp`**: force binary `fopen` mode;
  on Windows text mode corrupted the binary wavefunction/charge files.
- **`source/source_esolver/esolver_ks_lcao.cpp`**: guard the dereference of the
  DeePKS integrator `overlap_orb_alpha` (null when DeePKS is off) — it is only
  built for DeePKS runs.
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
  `gnu`/`intel`/`gcc-mkl` variants. After the build, `build_abacus_windows.sh`
  bundles the dependent DLLs next to `abacus.exe` (so it runs without `PATH`
  set, including under `mpiexec` redirected to a file) and — for MPI builds —
  drops an `mpirun`→`mpiexec` shim that pins `OMP_NUM_THREADS=1`, letting the
  unmodified harness drive MS-MPI.

## Known limitations / not yet ported

- ELPA, PEXSI, hybrid functionals (LibRI/LibComm), DeePKS/ML-KEDF, LibXC
  (so meta-GGA/SCAN), GPU (CUDA/ROCm), DSP — all disabled. Test cases needing
  them are expected to fail (e.g. `scf_metagga`, `scf_out_chg_tau`).
- **Serial gamma-only LCAO is buggy.** The `gamma_only` LCAO path gives a
  wrong (self-consistently converged) energy in a *serial* build — the same
  serial-only (`#ifndef __MPI`) reduction-gap family as the fixes above, but in
  the gamma H/density assembly and not yet located. The **MPI build is correct**
  (gamma matches the reference to ~1e-11, even on a single rank), so run LCAO
  gamma-only cases under MPI (`mpiexec -n 1 abacus` suffices). Multi-k serial
  LCAO is unaffected.

### `pw_seed` is not bit-reproducible across platforms

The random wavefunction initializer (`pw_seed > 0`) uses C `std::rand()`, whose
sequence and `RAND_MAX` are implementation-defined (e.g. 32767 on the Windows
CRT vs 2^31-1 on glibc). So for a given `pw_seed`, the *initial* wavefunctions
differ between Windows and Linux. For almost all systems the SCF converges to
the same state regardless of initialization, so results still match. But a few
**init-sensitive** cases (near-degenerate / charged / fixed-spin systems, e.g.
`tests/01_PW/078_PW_S2_elec_add`) can settle into a different near-degenerate
solution, so energy/force differ from the Linux-generated `result.ref`. This is
**not a code bug** — both states are valid converged solutions (the reference
state is reachable on Windows with a different seed). A proper cross-platform
fix would replace `std::rand` with a bit-portable generator (e.g. `std::mt19937`)
and regenerate the `pw_seed` references; that is left as a separate, upstream
change because it alters the sequence on Linux too.

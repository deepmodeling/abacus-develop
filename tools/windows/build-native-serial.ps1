<#
.SYNOPSIS
    Configure and build the *native* Windows serial PW-DFT version of ABACUS
    (Phase 1 of the native-Windows port: no MPI, no LCAO, no ELPA/PEXSI/hybrid).

.DESCRIPTION
    This drives a native CMake build with the MinGW-w64 GCC toolchain. MinGW is
    chosen over MSVC because it ships the POSIX headers ABACUS uses
    (unistd.h, fcntl.h, sys/stat.h, ...) and compiles the GCC __attribute__
    code unchanged, so the source needs only minimal portability shims.

    It does NOT replace the WSL2 installer (install-abacus.bat) -- that remains
    the recommended way to *run* full-featured ABACUS on Windows. This script is
    for developing/maintaining the native build.

    Required tools on PATH (or pass paths via parameters):
      - cmake (>= 3.16) and a generator (Ninja recommended, or "MinGW Makefiles")
      - g++ / gcc (MinGW-w64)
    Required libraries (native Windows builds):
      - BLAS + LAPACK   (e.g. OpenBLAS)
      - FFTW3           (double precision)
    The easiest consistent source for all of the above is MSYS2:
      pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake \
                mingw-w64-x86_64-ninja mingw-w64-x86_64-openblas \
                mingw-w64-x86_64-fftw
    Then run this script from a "MSYS2 MinGW 64-bit" shell's environment, or
    point -PrefixPath at the MSYS2 mingw64 prefix (e.g. C:\msys64\mingw64).

.PARAMETER BuildDir
    Out-of-source build directory. Default: build_win_serial_pw

.PARAMETER PrefixPath
    Extra CMAKE_PREFIX_PATH entries (semicolon-separated) where BLAS/LAPACK/FFTW3
    live, e.g. "C:\msys64\mingw64".

.PARAMETER Generator
    CMake generator. Default: "Ninja". Alternative: "MinGW Makefiles".

.PARAMETER Jobs
    Parallel build jobs. Default: number of logical processors.

.EXAMPLE
    ./build-native-serial.ps1 -PrefixPath "C:\msys64\mingw64"
#>
[CmdletBinding()]
param(
    [string]$BuildDir   = "build_win_serial_pw",
    [string]$PrefixPath = "",
    [string]$Generator  = "Ninja",
    [int]   $Jobs       = $env:NUMBER_OF_PROCESSORS
)

$ErrorActionPreference = "Stop"

# Repo root = two levels up from this script (tools/windows/ -> repo root)
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Write-Host "[*] Repo root : $RepoRoot"
Write-Host "[*] Build dir : $BuildDir"
Write-Host "[*] Generator : $Generator"

# --- sanity checks ---------------------------------------------------------
foreach ($tool in @("cmake", "g++")) {
    if (-not (Get-Command $tool -ErrorAction SilentlyContinue)) {
        throw "Required tool '$tool' not found on PATH. See the .DESCRIPTION header for setup (MSYS2 is recommended)."
    }
}

# --- configure -------------------------------------------------------------
$cmakeArgs = @(
    "-S", $RepoRoot,
    "-B", $BuildDir,
    "-G", $Generator,
    "-DCMAKE_BUILD_TYPE=Release",
    # Phase 1 scope: serial, plane-wave only.
    "-DENABLE_MPI=OFF",
    "-DENABLE_LCAO=OFF",
    "-DUSE_OPENMP=OFF",       # start minimal; FFTW3_OMP not required when OFF
    "-DUSE_ELPA=OFF",
    "-DENABLE_PEXSI=OFF",
    "-DENABLE_LIBRI=OFF",
    "-DENABLE_MLALGO=OFF",
    "-DUSE_CUDA=OFF",
    "-DBUILD_TESTING=OFF",
    "-DCMAKE_CXX_COMPILER=g++",
    "-DCMAKE_C_COMPILER=gcc"
)
if ($PrefixPath -ne "") {
    $cmakeArgs += "-DCMAKE_PREFIX_PATH=$PrefixPath"
}

Write-Host "`n[*] Configuring..."
& cmake @cmakeArgs
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed." }

# --- build -----------------------------------------------------------------
Write-Host "`n[*] Building (jobs=$Jobs)..."
& cmake --build $BuildDir --parallel $Jobs
if ($LASTEXITCODE -ne 0) { throw "Build failed." }

Write-Host "`n[+] Build complete. Look for abacus_pw_ser(.exe) under: $BuildDir"

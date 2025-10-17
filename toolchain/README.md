# ABACUS Toolchain

[![Version](https://img.shields.io/badge/version-2025.3-blue.svg)](https://github.com/deepmodeling/abacus-develop/tree/develop/toolchain)
[![License](https://img.shields.io/badge/license-GPL--compatible-green.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)]()

> **Automated dependency management and compilation toolchain for ABACUS**

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Quick Start](#quick-start)
- [Installation Methods](#installation-methods)
- [Supported Toolchains](#supported-toolchains)
- [Dependencies](#dependencies)
- [GPU Support](#gpu-support)
- [Troubleshooting](#troubleshooting)
- [Advanced Usage](#advanced-usage)
- [Developer Guide](#developer-guide)
- [License](#license)
- [Contributing](#contributing)

## Overview

The ABACUS Toolchain is an automated build system inspired by the cp2k-toolchain that simplifies the compilation and installation of ABACUS and its dependencies. It supports both online and offline installation modes, multiple compiler toolchains, and provides a streamlined path from dependency installation to ABACUS compilation.

### Main Developer

**[QuantumMisaka](https://github.com/QuantumMisaka)** (Zhaoqing Liu)  
*Peking University, CCME*

## Features

- ✅ **Multiple Toolchain Support**: GNU, Intel OneAPI, AMD AOCC/AOCL
- ✅ **Flexible Installation**: Online/offline modes with automatic dependency resolution  
- ✅ **GPU Acceleration**: CUDA support for NVIDIA GPUs with ELPA and cuSolverMP
- ✅ **MPI Implementations**: OpenMPI, MPICH, Intel MPI support
- ✅ **Math Libraries**: OpenBLAS, Intel MKL, AMD AOCL integration
- ✅ **Advanced Features**: LibRI, LibComm and MLALGO support
- ✅ **Resumable Installation**: Interrupt and resume capability
- ✅ **Environment Management**: Automatic setup file generation

## Quick Start

### Prerequisites

- **GCC**: Version ≥ 5.0 (recommended ≥ 7.3.0)
- **Internet Connection**: For online installation mode
- **System Libraries**: Basic development tools (see [System Requirements](#system-requirements))

### Basic Installation

For new users, start with one of these pre-configured toolchains:

```bash
# GNU toolchain (GCC + OpenMPI + OpenBLAS)
./toolchain_gnu.sh

# Intel toolchain (Intel compilers + Intel MPI + MKL)
./toolchain_intel.sh

# AMD toolchain options
./toolchain_gcc-aocl.sh    # GCC + AMD AOCL
./toolchain_aocc-aocl.sh   # AMD AOCC + AOCL
```

### Build ABACUS

After successful toolchain installation:

```bash
# For GNU toolchain
./build_abacus_gnu.sh

# For Intel toolchain  
./build_abacus_intel.sh

# For AMD toolchains
./build_abacus_gcc-aocl.sh
./build_abacus_aocc-aocl.sh
```

### Environment Setup

```bash
# Source the generated environment
source install/setup

# Or use the generated ABACUS environment
source abacus_env.sh
```

## Installation Methods

### Online Installation

Downloads packages automatically from official sources:

```bash
./toolchain_gnu.sh  # Uses system package managers and official repositories
```

**Package Sources:**
- **Build Tools:**
  - [GCC](https://mirrors.tuna.tsinghua.edu.cn/gnu/gcc/) - GNU Compiler Collection
  - [CMake](https://cmake.org/download/) - Cross-platform build system
- **MPI Libraries:**
  - [OpenMPI](https://download.open-mpi.org/release/open-mpi/) - Open source MPI implementation
  - [MPICH](https://www.mpich.org/downloads/) - High-performance MPI implementation
- **Math Libraries:**
  - [OpenBLAS](https://github.com/xianyi/OpenBLAS/releases) - Optimized BLAS library
  - [ScaLAPACK](http://www.netlib.org/scalapack/) - Scalable Linear Algebra PACKage
- **Scientific Libraries:**
  - [FFTW](http://www.fftw.org/) - Fast Fourier Transform library
  - [LibXC](https://www.tddft.org/programs/libxc/) - Exchange-correlation functionals
  - [ELPA](https://elpa.mpcdf.mpg.de/) - Eigenvalue solver
- **Advanced Features:**
  - [LibTorch](https://download.pytorch.org/libtorch/cpu/) - PyTorch C++ API
  - [LibNPY](https://github.com/llohse/libnpy) - NumPy I/O for C++
  - [LibRI](https://github.com/abacusmodeling/LibRI) - Resolution of Identity library
  - [LibComm](https://github.com/abacusmodeling/LibComm) - Communication library
  - [NEP](https://github.com/brucefan1983/NEP_CPU) - Neuroevolution Potential
  - [Cereal](https://github.com/USCiLab/cereal) - C++ serialization library
  - [RapidJSON](https://github.com/Tencent/rapidjson) - Fast JSON parser/generator
- **Reference mirror:** [CP2K static downloads](https://www.cp2k.org/static/downloads)

### Offline Installation

For air-gapped systems or unreliable internet:

```bash
# 1. Create build directory and download packages
mkdir build
# Download required packages to build/ directory with proper naming
# e.g., fftw-3.3.10.tar.gz, openmpi-5.0.8.tar.bz2

# 2. Run toolchain (will detect local packages)
./toolchain_gnu.sh
```

### Hybrid Installation

Mix online and offline packages as needed - the toolchain automatically detects locally available packages and downloads missing ones.

## Supported Toolchains

### GNU Toolchain
- **Compilers**: System GCC (≥5.0)
- **MPI**: OpenMPI or MPICH
- **Math**: OpenBLAS + ScaLAPACK
- **Features**: Most stable, widely compatible

### Intel Toolchain
- **Compilers**: Intel OneAPI (icx/icpx/ifx or classic icc/icpc/ifort)
- **MPI**: Intel MPI
- **Math**: Intel MKL
- **Features**: Optimized performance, EXX support

### AMD Toolchain
- **Compilers**: AMD AOCC or GCC
- **Math**: AMD AOCL (Optimized math libraries)
- **Features**: AMD processor optimization

## Dependencies

### Supported Packages

| Package | Version (main/alt) | Purpose | License | Default |
|---------|-------------------|---------|---------|---------|
| **Build Tools** |||||
| CMake | 3.31.7 / 3.30.5 | Build system | BSD-3-Clause | Install |
| GCC | 13.2.0 / 11.4.0 | C/C++ compiler | GPL-3.0-or-later WITH GCC-exception-3.1 | Install |
| **MPI Libraries** |||||
| OpenMPI | 5.0.8 / 4.1.6 | MPI implementation | BSD-3-Clause-Open-MPI | Install |
| MPICH | 4.3.1 / 4.1.0 | Alternative MPI | mpich2 (BSD-like) | Alternative |
| **Math Libraries** |||||
| OpenBLAS | 0.3.30 / 0.3.27 | Linear algebra | BSD-3-Clause | Install |
| ScaLAPACK | 2.2.2 / 2.2.1 | Parallel linear algebra | BSD-3-Clause | Install |
| **Scientific Libraries** |||||
| FFTW | 3.3.10 / 3.3.10 | Fast Fourier Transform | GPL-2.0-or-later | Install |
| LibXC | 7.0.0 / 6.2.2 | Exchange-correlation | MPL-2.0 | Install |
| ELPA | 2025.06.001 / 2024.05.001 | Eigenvalue solver | LGPL-3.0-only | Install |
| **Advanced Features** |||||
| Cereal | master | C++ Serialization | BSD | Install |
| RapidJSON | master | JSON parsing | MIT | Install |
| LibRI | master | EXX calculations | GPL-3.0 | Install |
| LibComm | master | EXX calculations | GPL-3.0 | Install |
| LibTorch | 2.1.2 / 1.12.1 | MLALGO support | BSD-3-Clause | Optional |
| LibNPY | 1.0.1 / 1.0.1 | NumPy I/O | MIT | Optional |
| NEP | main | Neural network potential | MIT | Optional |

### System Requirements

Install system dependencies using provided scripts:

```bash
# Ubuntu/Debian
sudo ./root_requirements/install_requirements_ubuntu.sh

# Fedora/RHEL/CentOS  
sudo ./root_requirements/install_requirements_fedora.sh

# Generic
sudo ./root_requirements/install_requirements.sh
```

## GPU Support

### CUDA Support for NVIDIA GPUs

#### Basic GPU Support

Add to your build script:

```bash
cmake -B $BUILD_DIR \
    -DUSE_CUDA=ON \
    -DCMAKE_CUDA_COMPILER=/path/to/cuda/bin/nvcc \
    # ... other options
```

#### Multi-GPU with ELPA

1. **Configure toolchain with CUDA:**
```bash
export CUDA_PATH=/path/to/CUDA
./toolchain_gnu.sh --enable-cuda --gpu-ver=70  # For V100 (compute capability 7.0)
```

2. **Build with ELPA GPU support:**
```bash
cmake -B $BUILD_DIR \
    -DUSE_CUDA=ON \
    -DUSE_ELPA=ON \
    # ... other options
```

#### Multi-GPU with cuSolverMP

1. **Check or install cuSolverMP manually:**
One may use NVIDIA HPC_SDK as an easy way to install cuSolverMP.

2. **Install dependencies normally:**
```bash
./toolchain_gnu.sh
```

3. **Build with cuSolverMP:**
```bash
cmake -B $BUILD_DIR \
    -DUSE_CUDA=ON \
    -DENABLE_CUSOLVERMP=ON \
    -DCAL_CUSOLVERMP_PATH=/path/to/math_libs/lib \
    # ... other options
```

3. **Set environment variables:**
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hpcx/ucc/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/hpcx/ucx/lib  
export CPATH=$CPATH:/path/to/math_libs/include
```

**Note**: cuSolverMP requires NVIDIA HPC SDK or system installation via package manager.

## Troubleshooting

### Common Issues

#### Intel OneAPI Problems

**OneAPI 2025.0 Compatibility:**
- LibRI compatibility issues ([#6190](https://github.com/deepmodeling/abacus-develop/issues/6190))
- Solution: Using the patch from Cereal and the master version of Cereal to fix the compatibility issue (included in toolchain).

**ELPA on AMD servers with Intel compilers:**
```bash
# Use Intel classic compilers instead
./toolchain_intel.sh --with-intel-classic=yes
```

#### OpenMPI Issues

**Version 5 compatibility problems:**
```bash
# Use OpenMPI v4 instead
./toolchain_gnu.sh --with-openmpi-4th=yes
```

**LibComm compilation with OpenMPI:**
- Fixed in toolchain 2025.2 which downlo the master branch of LibComm
- Alternative: Use MPICH or Intel MPI

#### Shell and Permission Issues

**Line ending problems:**
```bash
./pre_set.sh  # Fixes line endings and permissions
# Or manually:
dos2unix *.sh
chmod +x *.sh
```

#### Library Version Issues

**LibTorch GLIBC errors:**
- Change version from 2.1.2 to 1.12.1
- Use --package-version libtorch:alt when calling toolchain

**DeepMD GLIBC errors:**
- Requires GCC ≥ 11.3.1
- Upgrade system GCC or use newer toolchain

### Getting Help

1. **Check logs**: Look in `build/PKG_NAME/make.log` for compilation errors
2. **Reduce parallelism**: Use `NPROCS_OVERWRITE=N` environment variable to limit parallel processes
3. **System libraries**: Use `--with-PKG=system` for system-installed packages
4. **Clean installation**: Remove `install/` and `build/` directories to restart
5. **Certificate issues**: Use `DOWNLOAD_CERT_POLICY=skip` for download problems

## Advanced Usage

### Package-Specific Options

```bash
# Use Intel MKL instead of installing OpenBLAS
./toolchain_gnu.sh --with-mkl=system

# Use system FFTW instead of installing
./toolchain_gnu.sh --with-fftw=system

# Specify custom package installation path
./toolchain_gnu.sh --with-fftw=/path/to/custom/fftw
```

### Execution Mode Control

```bash
# Test configuration without actual installation (recommended for first run)
./toolchain_gnu.sh --dry-run

# Only download packages without building (useful for offline preparation)
./toolchain_gnu.sh --pack-run
```

### Environment Variable Configuration

The toolchain supports several environment variables for advanced configuration:

#### Download Certificate Verification

Control SSL/TLS certificate verification during package downloads:

```bash
# Strict mode: Always verify certificates (secure)
export DOWNLOAD_CERT_POLICY=strict
./toolchain_gnu.sh

# Smart mode: Try secure first, fallback if needed (default)
export DOWNLOAD_CERT_POLICY=smart  # or leave unset
./toolchain_gnu.sh

# Skip mode: Skip certificate verification (legacy compatibility)
export DOWNLOAD_CERT_POLICY=skip
./toolchain_gnu.sh
```

**Smart Mode Behavior**: The default `smart` mode first attempts secure downloads with certificate verification. If this fails (e.g., due to corporate firewalls or outdated certificates), it automatically falls back to skipping certificate verification while providing clear user feedback.

#### Parallel Compilation Control

Override the automatic CPU core detection for compilation:

```bash
# Use 8 cores for compilation (useful for resource-limited systems)
export NPROCS_OVERWRITE=8
./toolchain_gnu.sh

# Use single core for debugging compilation issues
export NPROCS_OVERWRITE=1
./toolchain_gnu.sh

# Or specify inline
NPROCS_OVERWRITE=4 ./toolchain_gnu.sh --with-gcc --with-openmpi
```

**Use Cases**:
- **Resource-limited systems**: Reduce parallelism to avoid memory exhaustion
- **Shared servers**: Limit resource usage to be considerate of other users
- **CI/CD environments**: Match container resource limits
- **Debugging**: Use single-core compilation for clearer error messages

### Legacy Script Options

The deprecated `install_abacus_toolchain.sh` supports additional options:

| Option | Description | Availability |
|--------|-------------|--------------|
| `--dry-run` | Test configuration without installation | ✅ New & Legacy |
| `--pack-run` | Download packages without building | ✅ New & Legacy |
| `--no-check-certificate` | Skip SSL certificate verification | ⚠️ Legacy only (use `DOWNLOAD_CERT_POLICY=skip`) |
| `-j N` | Limit parallel compilation processes | ⚠️ Legacy only (use `NPROCS_OVERWRITE=N`) |

> **Migration Note**: The new toolchain system (`toolchain_*.sh` scripts) is recommended over the legacy `install_abacus_toolchain.sh`. Legacy options like `--no-check-certificate` and `-j N` are replaced by environment variables `DOWNLOAD_CERT_POLICY` and `NPROCS_OVERWRITE` respectively.

### Environment Management

The toolchain generates several setup files:

- `install/setup`: Main environment setup
- `build/setup_PKG`: Individual package environments  
- `abacus_env.sh`: ABACUS-specific environment (generated by build scripts)

## Developer Guide

### Toolchain Architecture

The toolchain follows a modular design with staged dependency installation:

```
scripts/
├── stage0/          # Compilers and build tools
├── stage1/          # MPI implementations  
├── stage2/          # Math libraries (BLAS, LAPACK)
├── stage3/          # Scientific libraries (FFTW, LibXC, ELPA)
├── stage4/          # Advanced features (LibTorch, LibRI)
└── lib/             # Core toolchain libraries
```

### Key Components

| File | Purpose |
|------|---------|
| `install_abacus_toolchain_new.sh` | Main orchestration script (new version) |
| `install_abacus_toolchain.sh` | Legacy main script (deprecated) |
| `toolchain_*.sh` | Frontend scripts for specific toolchains |
| `scripts/lib/config_manager.sh` | Configuration management |
| `scripts/lib/package_manager.sh` | Package installation logic |
| `scripts/lib/user_interface.sh` | User interaction and output |
| `scripts/common_vars.sh` | Shared variables and defaults |
| `scripts/tool_kit.sh` | Utility functions and macros |
| `scripts/parse_if.py` | Parser for IF_XYZ constructs |
| `checksums.sha256` | Pre-calculated SHA256 checksums for packages |

### Script Structure Details

**Individual Package Scripts**: Each `scripts/stage*/install_PKG.sh` script is relatively independent and should:

1. **Generate setup files**: Write to both `build/setup_PKG` and `install/setup`
   - `build/setup_PKG`: Variables for toolchain compilation and arch file flags
   - `install/setup`: Environment setup for compiling/running ABACUS

2. **Handle dependencies**: May depend on other libraries being installed with correct environment variables

3. **Use toolkit macros**: Leverage functionality from `scripts/tool_kit.sh` for common operations

### Package Installation Scripts

Each `scripts/stage*/install_PKG.sh` script:

1. **Downloads** the package (if not available locally)
2. **Configures** build with appropriate flags
3. **Compiles** with error handling and logging
4. **Installs** to the toolchain directory
5. **Generates** setup files for environment configuration

### Configuration System

#### Package Control Options (`--with-PKG`)

The `--with-PKG` options control how a package is going to be installed:

- `--with-PKG=install` (or `--with-PKG` alone): Compile and install from source downloaded (default)
- `--with-PKG=system`: Link to locations provided by system search paths
- `--with-PKG=/path/to/pkg`: Link to locations provided by the user (custom path)
- `--with-PKG=no`: Skip package installation entirely

**System Search Paths**: When using `system` mode, the installation script searches in:
- `LD_LIBRARY_PATH`, `LD_RUN_PATH`, `LIBRARY_PATH`
- `/usr/local/lib64`, `/usr/local/lib`, `/usr/lib64`, `/usr/lib`
- For MKL libraries: `MKLROOT` environment variable

**Troubleshooting System Libraries**: If `--with-PKG=system` cannot find the library:
1. Use `module show PKG` to see module-defined paths
2. Find the root installation directory manually
3. Use `--with-PKG=/path/to/pkg` to specify exact location

#### Feature Control Options (`--enable-FEATURE`)

The `--enable-FEATURE` options control whether optional features are enabled:

- `--enable-FEATURE=yes` (or `--enable-FEATURE` alone): Enable the feature
- `--enable-FEATURE=no`: Disable the feature

#### Mode Selection (`PKG_MODE` Variables)

For packages serving the same purpose, mode variables act as selectors:

- `--mpi-mode=openmpi|mpich|intelmpi`: Choose MPI implementation
- `--math-mode=openblas|mkl|aocl`: Choose math library

**Note**: While `--with-PKG` controls the installation method, the `PKG_MODE` variable picks which package to actually use, providing maximum flexibility.

### Adding New Packages

1. **Create installation script**: `scripts/stageN/install_newpkg.sh`
2. **Add to stage script**: Include in `scripts/stageN/install_stageN.sh`
3. **Update configuration**: Add options to `config_manager.sh`
4. **Add version info**: Update `scripts/package_versions.sh`
5. **Test thoroughly**: Verify with different toolchain combinations

### Advanced Developer Features

#### The IF_XYZ Constructs

The toolchain uses a special syntax construct for conditional compilation flags:

```shell
IF_XYZ(A | B)
```

This construct is parsed by `scripts/parse_if.py`:
- Evaluates to *A* if *XYZ* is passed as command line option
- Evaluates to *B* if *XYZ* is not passed

**Nested Constructs**: The `IF_XYZ(A|B)` construct can be nested:

```shell
IF_XYZ(IF_ABC(flag1|flag2) | flag3)
```

This parses to:
- *flag1* if both *XYZ* and *ABC* are present
- *flag2* if only *XYZ* is present  
- *flag3* if neither is present

#### Portability Requirements

**Compiler Flag Filtering**: Always pass compiler flags through compatibility filters:

```shell
# Filter flags for GCC compatibility
CFLAGS="$(allowed_gcc_flags $CFLAGS)"
FCFLAGS="$(allowed_gfortran_flags $FCFLAGS)"
```

**IF_XYZ with Flag Filtering**: Since filters don't work with IF_XYZ constructs, break them down:

```shell
# Instead of: FCFLAGS="IF_XYZ(flag1 flag2 | flag3 flag4)"
XYZ_TRUE_FLAGS="flag1 flag2"
XYZ_FALSE_FLAGS="flag3 flag4"
# Apply filtering
XYZ_TRUE_FLAGS="$(allowed_gcc_flags $XYZ_TRUE_FLAGS)"
XYZ_FALSE_FLAGS="$(allowed_gcc_flags $XYZ_FALSE_FLAGS)"
# Reconstruct
FCFLAGS="IF_XYZ($XYZ_TRUE_FLAGS | $XYZ_FALSE_FLAGS)"
```

**Fortran Module Checking**: Check intrinsic Fortran modules with:

```shell
check_gfortran_module module_name
```

**Avoid Hard Coding**: Use common variables instead of hard-coded paths:

```shell
# Good practice
./configure --prefix=some_dir CC=${MPICC} FC=${MPIFC}
# Avoid
./configure --prefix=some_dir CC=mpicc FC=mpif90
```

### Best Practices

- **Reuse toolkit functions**: Use macros from `scripts/tool_kit.sh`
- **Modular functionality**: Add new functionality as macros in `scripts/tool_kit.sh` rather than inline code
- **Portable compiler flags**: Filter through `allowed_gcc_flags` and `allowed_gfortran_flags`
- **Environment variables**: Use `${VAR:-default}` pattern for configurable defaults
- **Lock files**: Create completion markers for resumable installation
- **Separate directories**: Install each package in its own directory
- **Error handling**: Provide clear error messages and recovery suggestions

## License

The ABACUS Toolchain downloads and installs only [GPL-compatible](https://www.gnu.org/licenses/gpl-faq.html#WhatDoesCompatMean) packages. All included packages maintain their original licenses as listed in the Dependencies section above.

**License Compatibility**: All packages use GPL-compatible licenses including BSD, MIT, LGPL, MPL-2.0, and GPL variants, ensuring seamless integration with GPL-licensed software.

**Note**: Proprietary packages like Intel OneAPI (MKL/Compiler/MPI) and AMD AOCC/AOCL are supported but must be installed separately by the user.

## Contributing

We welcome contributions to improve the ABACUS Toolchain! Here's how you can help:

### Reporting Issues

1. **Search existing issues** before creating new ones
2. **Provide detailed information**:
   - Operating system and version
   - Compiler versions
   - Complete error messages and logs
   - Steps to reproduce

### Contributing Code

1. **Fork the repository** and create a feature branch
2. **Follow coding standards**:
   - Use consistent shell scripting style
   - Add comments for complex logic
   - Test with multiple toolchain combinations
3. **Update documentation** for new features
4. **Submit pull request** with clear description

### Development Setup

```bash
# Clone the repository
git clone https://github.com/deepmodeling/abacus-develop.git
cd abacus-develop/toolchain

# Test your changes
./toolchain_gnu.sh --dry-run
```

### Areas for Contribution

- 🔧 **New package support**: Add support for additional scientific libraries
- 🐛 **Bug fixes**: Resolve compatibility issues and installation problems  
- 📚 **Documentation**: Improve guides and troubleshooting information
- 🧪 **Testing**: Expand test coverage for different systems and configurations
- 🚀 **Performance**: Optimize installation speed and resource usage

---

**For questions, issues, or contributions, please visit the [ABACUS GitHub repository](https://github.com/deepmodeling/abacus-develop).**

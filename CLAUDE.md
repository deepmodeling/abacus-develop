# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About ABACUS

ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is an open-source electronic structure package based on density functional theory (DFT). It supports both plane wave (PW) and linear combination of atomic orbitals (LCAO) basis sets with norm-conserving pseudopotentials. The package provides various DFT functionals (LDA, GGA, meta-GGA, hybrid), geometry optimization, ab-initio molecular dynamics, and advanced features like DFT+U, VdW corrections, and machine-learning-assisted DFT methods.

Documentation: https://abacus.deepmodeling.com/

## Build System

ABACUS uses CMake as its build system. The main build configuration is in `CMakeLists.txt` at the root.

### Basic Build Commands

```bash
# Configure with default options (LCAO enabled, MPI enabled)
cmake -B build

# Build with specific number of processors
cmake --build build -j${number_of_processors}

# Install (required before running tests)
cmake --install build
```

### Common Build Options

- `-DBUILD_TESTING=ON` - Enable unit tests (required for testing)
- `-DENABLE_LCAO=ON/OFF` - Enable/disable LCAO calculations (default: ON)
- `-DUSE_CUDA=ON` - Enable CUDA support
- `-DUSE_ROCM=ON` - Enable ROCm/HIP support
- `-DUSE_OPENMP=ON/OFF` - Enable OpenMP (default: ON)
- `-DENABLE_COVERAGE=ON` - Enable code coverage (sets Debug build)
- `-DENABLE_ASAN=ON` - Enable AddressSanitizer (sets RelWithDebInfo build)
- `-DCMAKE_BUILD_TYPE=Debug/Release/RelWithDebInfo` - Build type
- `-DENABLE_DEEPKS=ON` - Enable DeePKS functionality
- `-DENABLE_LIBXC=ON` - Enable LibXC functionality
- `-DUSE_ELPA=ON/OFF` - Enable ELPA diagonalization (default: ON)
- `-DENABLE_LIBRI=ON` - Enable EXX with LibRI
- `-DENABLE_PAW=ON` - Enable PAW calculations
- `-DMKLROOT=$MKLROOT` - Use Intel MKL for math libraries

### Build Variants

The executable name depends on build configuration:
- `abacus` - LCAO enabled, MPI enabled (default)
- `abacus_pw` - LCAO disabled, MPI enabled
- `abacus_serial` - LCAO enabled, MPI disabled
- `abacus_pw_serial` - LCAO disabled, MPI disabled

## Testing

### Unit Tests

Unit tests use GoogleTest framework and are located in `source/*/test` directories.

```bash
# Build with unit tests enabled
cmake -B build -DBUILD_TESTING=ON
cmake --build build -j${nproc}
cmake --install build  # REQUIRED before running tests

# Run all tests
cd build
ctest -V

# Run specific test by name
ctest -R <test-name>

# Run tests matching pattern
ctest -R <pattern>

# Build specific unit test
cmake --build build -j${nproc} --target ${unit_test_name}
```

Unit test executables are located in `build/source/${module_name}/test` (or `test_parallel`, `test_pw` for specialized tests).

### Integration Tests

Integration tests are in `tests/integrate/` directory. Each test case is a subdirectory with input files (INPUT, STRU, KPT) and reference results (result.ref).

```bash
# Run all integration tests
cd tests/integrate
bash Autotest.sh

# Run specific test case
cd tests/integrate/<test_case_name>
bash ../Single_job.sh

# Run with custom parameters
bash Autotest.sh -a /path/to/abacus -n 4 -t 0.0000001

# Generate reference results for new test
bash Autotest.sh -g -r <test_case_name>
```

Key Autotest.sh options:
- `-a <path>` - ABACUS executable path (default: abacus)
- `-n <num>` - Number of MPI processes (default: 4)
- `-t <threshold>` - Energy threshold in eV (default: 0.0000001)
- `-c <accuracy>` - Check accuracy (default: 8)
- `-g` - Generate reference results
- `-r <regex>` - Test case name regex filter
- `-f <file>` - Test cases file (default: CASES_CPU.txt)

## Code Architecture

### Module Structure

The source code is organized into modules under `source/`:

**Core Infrastructure:**
- `module_base/` - Mathematical library interfaces (BLAS, LAPACK, ScaLAPACK), custom data structures (matrix, vector), parallelization (MPI, OpenMP), utilities (timer, memory), and global parameters
- `module_container/` - Container module for data storage and operations across different architectures
- `module_parameter/` - Input parameters and global variables

**Basis Sets:**
- `module_basis/module_pw/` - Plane wave basis data structures and methods
- `module_basis/module_nao/` - Numerical atomic orbital basis for LCAO two-center integrals
- `module_basis/module_ao/` - Legacy atomic orbital basis (being refactored)

**Cell and Structure:**
- `module_cell/` - Unit cell definition, operations, and pseudopotential reading
- `module_cell/module_neighbor/` - Neighbor finding algorithms
- `module_cell/module_symmetry/` - Symmetry operations
- `module_cell/module_paw/` - PAW calculations

**Electronic Structure:**
- `module_elecstate/` - Electronic state definition and operations
- `module_elecstate/module_charge/` - Charge density calculation and mixing
- `module_elecstate/potentials/` - Potential calculations (Hartree, XC, local pseudopotential)
- `module_psi/` - Wave function definition and operations

**Hamiltonians:**
- `module_hamilt_general/` - General Hamiltonian components used in both PW and LCAO
  - `module_ewald/` - Ewald summation
  - `module_surchem/` - Surface charge correction
  - `module_vdw/` - van der Waals corrections
  - `module_xc/` - Exchange-correlation energy and potential
- `module_hamilt_pw/` - PW-specific Hamiltonians
  - `hamilt_pwdft/operator_pw/` - PW-DFT operators
  - `hamilt_ofdft/` - Orbital-free DFT
  - `hamilt_stodft/` - Stochastic DFT
- `module_hamilt_lcao/` - LCAO-specific Hamiltonians
  - `hamilt_lcaodft/operator_lcao/` - LCAO-DFT operators
  - `module_hcontainer/` - Hamiltonian matrix storage
  - `module_gint/` - Grid integration
  - `module_dftu/` - DFT+U implementation
  - `module_deepks/` - DeepKS integration
  - `module_tddft/` - Time-dependent DFT

**Solvers and Drivers:**
- `module_hsolver/` - Hamiltonian diagonalization methods (CG, Davidson for PW; ScaLAPACK, ELPA for LCAO)
- `module_esolver/` - Task-specific workflow drivers (KS-DFT, SDFT, OFDFT, LJ, DP, etc.)
- `module_md/` - Molecular dynamics
- `module_relax/` - Structural optimization (relax_new for simultaneous cell+ion, relax_old for separate)
- `module_lr/` - Linear response calculations

**I/O:**
- `module_io/` - INPUT file reading and property output (band structure, DOS, charge density, etc.)

### Program Flow

1. **main.cpp** - Entry point
   - Parse command-line arguments
   - Initialize MPI environment
   - Initialize FFTW threads (if OpenMP enabled)
   - Create Driver object and call init()

2. **Driver::init()** (driver.cpp)
   - Read input parameters via Driver::reading()
   - Call atomic_world() for main computation
   - Output timing and close logs
   - Generate JSON output

3. **Driver::reading()**
   - Read INPUT file parameters
   - Create output directory and log files
   - Split MPI communicators for diagonalization and k-point parallelization
   - Read wannier function parameters

4. **Driver::atomic_world()**
   - Call driver_run() for actual calculations
   - Output timing and memory statistics

5. **driver_run()** (driver_run.cpp)
   - Initialize ESolver (energy solver) based on calculation type
   - Run the calculation workflow (SCF, relaxation, MD, etc.)

### Key Design Patterns

**Basis Set Abstraction:** The code supports multiple basis sets (PW, LCAO) through polymorphism. Hamiltonian and solver modules have separate implementations for each basis type.

**Operator Pattern:** Hamiltonians are built from operator components (kinetic, local potential, nonlocal potential, etc.) that can be combined flexibly.

**MPI Parallelization:** Multi-level parallelization over k-points (pools), bands, and real-space grids. Communicators are split in Driver::reading().

**Module Independence:** Modules are designed to be relatively independent with clear interfaces, though some global state exists (being refactored).

## Code Style and Formatting

- Use `clang-format` with the `.clang-format` file in root directory
- Doxygen comments for documentation (Javadoc style preferred)
- Comment only in `.h` files for Doxygen visibility
- Use `@param` for parameters, `\f$...\f$` for inline formulas
- Pre-commit hooks available but not required locally (runs on CI)

## Development Workflow

### Adding New Features

1. Read relevant source files first before proposing changes
2. Follow existing code patterns in the module
3. Add unit tests in `source/${module}/test/` using GoogleTest
4. Add integration tests in `tests/integrate/` if needed
5. Update CMakeLists.txt if adding new source files or tests
6. Ensure code follows clang-format style

### Debugging

```bash
# Build with debug symbols
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j${nproc}
cmake --install build

# Use GDB
cd <input_directory>
gdb /path/to/abacus
# Set breakpoints and run

# Debug with MPI (attach GDB to rank 0)
mpirun -n 1 gdb abacus : -n 3 abacus

# Use AddressSanitizer for memory issues
cmake -B build -DENABLE_ASAN=1
cmake --build build -j${nproc}
cmake --install build
# Run normally - ASan will report issues
```

### Code Coverage

```bash
# Requires GCC compiler
cmake -B build -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build -j${nproc}
cmake --install build
cmake --build build --target test ARGS="-V --timeout 21600"
cd build
make lcov
# View build/lcov/html/all_targets/index.html
```

## Commit Message Format

Follow Conventional Commits specification:

```
<type>[optional scope]: <description>

[optional body]

[optional footer]
```

Types: `Feature`, `Fix`, `Docs`, `Style`, `Refactor`, `Perf`, `Test`, `Build`, `CI`, `Revert`

Example:
```
Fix(lcao): use correct scalapack interface

`pzgemv_` and `pzgemm_` used `double*` for alpha and beta parameters
but not `complex*`, this would cause error in GNU compiler.

Fix #753.
```

## Important Notes

- Always run `cmake --install build` after building and before running tests
- For PW calculations, set `pw_seed 1` in INPUT file for reproducible tests
- Integration tests should run in < 20 seconds (reduce atoms, k-points, ecutwfc, steps)
- Pseudopotential and orbital files go in `tests/PP_ORB/`
- Use relative paths in INPUT for `pseudo_dir` and `orb_dir`
- Never use `-uall` flag with `git status` (can cause memory issues on large repos)
- The main development branch is `develop`, not `main` or `master`

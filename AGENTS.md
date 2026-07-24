# AGENTS.md

This file provides guidance to agentic coding agents working with the ABACUS codebase.

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

## Code Style and Formatting

### Formatting

- Use `clang-format` with the `.clang-format` file in root directory
- Based on Microsoft style with customizations
- 4-space indentation, no tabs
- Left-aligned pointers (`int* ptr` not `int *ptr`)
- Sort includes and using declarations

### Documentation

- Doxygen comments for documentation (Javadoc style preferred)
- Comment only in `.h` files for Doxygen visibility
- Use `@param` for parameters, `\f$...\f$` for inline formulas

### Code Conventions

**Naming:**
- Classes: PascalCase (e.g., `HSolver`, `Matrix`)
- Functions: snake_case (e.g., `solve_hamiltonian`, `calculate_energy`)
- Variables: snake_case (e.g., `num_bands`, `ecutwfc`)
- Constants: UPPER_SNAKE_CASE (e.g., `MAX_ITERATIONS`)
- Private members: trailing underscore (e.g., `num_bands_`)

**Imports and Includes:**
- System headers first, then project headers
- Use angle brackets for system headers, quotes for project headers
- Group includes: C standard library, C++ standard library, external libraries, internal headers
- Sort includes alphabetically within groups

**Types:**
- Use `double` for floating-point by default, `float` only when memory is critical
- Use `std::complex<double>` for complex numbers
- Prefer `std::vector` over C-style arrays
- Use `size_t` for sizes and indices
- Use `const` and references where appropriate

**Error Handling:**
- Use assertions (`assert()`) for internal invariants
- Return error codes or use exceptions for recoverable errors
- Log errors using the existing logging infrastructure
- Validate input parameters at function boundaries

**Templates:**
- Use templates for generic algorithms and data structures
- Provide explicit instantiations for common types
- Use `template <>` for specializations

**Memory Management:**
- Use RAII principles with smart pointers (`std::unique_ptr`, `std::shared_ptr`)
- Avoid raw `new`/`delete` when possible
- Follow the Rule of Three/Five for classes managing resources

**Parallelization:**
- Use MPI for distributed memory parallelization
- Use OpenMP for shared memory parallelization
- Follow existing patterns for communicator splitting
- Be careful with global variables in parallel contexts

### Module Structure

The source code is organized into modules under `source/`:

**Core Infrastructure:**
- `module_base/` - Mathematical library interfaces, data structures, parallelization, utilities
- `module_container/` - Container module for data storage and operations
- `module_parameter/` - Input parameters and global variables

**Basis Sets:**
- `module_basis/module_pw/` - Plane wave basis
- `module_basis/module_nao/` - Numerical atomic orbital basis
- `module_basis/module_ao/` - Legacy atomic orbital basis

**Cell and Structure:**
- `module_cell/` - Unit cell definition and operations
- `module_cell/module_neighbor/` - Neighbor finding
- `module_cell/module_symmetry/` - Symmetry operations

**Electronic Structure:**
- `module_elecstate/` - Electronic state definition and operations
- `module_elecstate/module_charge/` - Charge density calculation and mixing
- `module_elecstate/potentials/` - Potential calculations
- `module_psi/` - Wave function definition and operations

**Hamiltonians:**
- `module_hamilt_general/` - General Hamiltonian components
- `module_hamilt_pw/` - PW-specific Hamiltonians
- `module_hamilt_lcao/` - LCAO-specific Hamiltonians

**Solvers and Drivers:**
- `module_hsolver/` - Hamiltonian diagonalization methods
- `module_esolver/` - Task-specific workflow drivers
- `module_md/` - Molecular dynamics
- `module_relax/` - Structural optimization

**I/O:**
- `module_io/` - INPUT file reading and property output

## Development Workflow

### Adding New Features

1. Read relevant source files first before proposing changes
2. Follow existing code patterns in the module
3. Add unit tests in `source/${module}/test/` using GoogleTest
4. Add integration tests in `tests/integrate/` if needed
5. Update CMakeLists.txt if adding new source files or tests
6. Ensure code follows clang-format style

### Important Notes

- Always run `cmake --install build` after building and before running tests
- For PW calculations, set `pw_seed 1` in INPUT file for reproducible tests
- Integration tests should run in < 20 seconds (reduce atoms, k-points, ecutwfc, steps)
- Pseudopotential and orbital files go in `tests/PP_ORB/`
- Use relative paths in INPUT for `pseudo_dir` and `orb_dir`
- The main development branch is `develop`, not `main` or `master`

### Commit Message Format

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
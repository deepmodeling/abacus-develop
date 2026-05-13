# ABACUS DFT+U / DeltaSpin PW Port

## Project

ABACUS electronic structure calculation software. The current branch `feat/dftu-pw-port-v2` ports DFT+U and DeltaSpin functionality to the PW (plane wave) basis set.

- Language: C++ (C++17)
- Build: CMake + make
- Parallel: MPI + OpenMP (optional CUDA)
- Testing: Integration tests in `tests/`, unit tests with GoogleTest

## Build

```bash
cmake --build build -j$(nproc)
```

Build artifact: `build/abacus_pw_para`

Current build configuration: Debug, MPI=ON, CUDA=ON, OpenMP=OFF

## Running Test Cases

```bash
cd tests/17_DS_DFTU/13_PW_DS_S4_XY
mpirun -np 4 /root/abacus-dftu-pw-port/build/abacus_pw_para
# or single process:
mpirun -np 1 /root/abacus-dftu-pw-port/build/abacus_pw_para
```

## Current Task

Debug the DeltaSpin feature causing ~1000 eV energy deviation under the PW basis set. See `DELTA_SPIN_DEBUG.md` for details.

### Fixed
- `nproc_in_pool` parallel reduction issue (locale double-counting when kpar>1)
- Pure DFT+U (without DeltaSpin) verified correct

### Pending Investigation
- Changes in `update_psi_charge_pw_cpu` and underlying functions between v3.7.1→v3.11.0
- SCF restart logic and charge mixing parameters
- HSolver diagonalization solver differences
- DFT+U and DeltaSpin interaction

## Code Structure (relevant parts)

```
source/source_lcao/module_deltaspin/   # DeltaSpin core (lambda loop, spin constraint)
source/source_pw/module_pwdft/         # PW-layer DeltaSpin/DFT+U interface
source/source_lcao/module_dftu/        # DFT+U module
source/source_esolver/                 # ESolver (SCF main loop)
source/source_hsolver/                 # Diagonalization solver
source/source_estate/                  # Electronic state (ElecState, Charge)
```

## Baseline Comparison

- Baseline branch: `tmp` (zdy), ABACUS v3.7.1
- Baseline build: `/root/abacus-zdy-tmp/build_zdy/abacus_pw`
- Baseline output: `tests/17_DS_DFTU/13_PW_DS_S4_XY/log2-1`

## Development Guidelines

- Compile and verify after each logical change
- Read relevant source files before fixing
- Retry the same error at most 2 times, then output root cause analysis
- Do not arbitrarily disable features (CUDA, MPI, etc.)
- Ask immediately when uncertain
- Use `[DS-DIAG]` tag for diagnostic output, for easy grep

# DFPT-PW Module

## Overview

This module implements Density Functional Perturbation Theory (DFPT) for 
plane-wave basis set in ABACUS. It allows calculation of phonon frequencies,
dielectric tensor, Born effective charges, and related properties.

**Note:** This code is currently in the design phase and has not been 
put into production yet. It may change in the future. Please use with caution.

## Directory Structure

```
module_dfpt/
├── README.md           # This file
├── dfpt_pw.h           # DFPT-PW main interface class
├── dfpt_pw.cpp         # DFPT-PW implementation
├── dfpt_pw_data.h      # PW-specific DFPT data container
├── dfpt_pw_data.cpp
├── dfpt_pert.h         # Perturbation construction
├── dfpt_pert.cpp
├── dfpt_stern.h        # Sternheimer equation solver
├── dfpt_stern.cpp
├── dfpt_rho.h          # First-order density handling
├── dfpt_rho.cpp
├── dfpt_phon.h         # Phonon/dynamical matrix
├── dfpt_phon.cpp
├── dfpt_q0.h           # q=0 special handling
├── dfpt_q0.cpp
├── dfpt_metal.h        # Metal system handling
├── dfpt_metal.cpp
└── CMakeLists.txt      # Build configuration
```

## Design Philosophy

### 1. Separation of Concerns
- **Data Layer**: `DFPT_PW_Data` stores all DFPT-related data
- **Algorithm Layer**: Individual classes handle specific algorithms
- **Interface Layer**: `DFPT_PW` provides a clean API to ESolver

### 2. Encapsulation
- All data members in `DFPT_PW_Data` are private
- Access is through getter/setter methods
- Pimpl idiom used to hide implementation details from ESolver

### 3. Reusability
- Uses existing ABACUS components: Psi, Charge_Mixing, Monkhorst-Pack
- q-point management via `ModuleCell::QList`
- Conjugate gradient via existing HSolver

### 4. KISS Principle
- Short, descriptive function and variable names
- Minimal dependencies between components
- Clear separation of PW-specific code

## Module Dependencies

```
DFPT_PW
    ├── DFPT_PW_Data     # Data container
    ├── DFPT_Pert        # Perturbation construction
    ├── DFPT_Stern       # Sternheimer solver
    ├── DFPT_Rho         # Density handling
    ├── DFPT_Phon        # Phonon calculation
    ├── DFPT_Q0          # q=0 special handling
    ├── DFPT_Metal       # Metal system handling
    └── ModuleCell::QList # q-point management
```

## Key Features

1. **Monochromatic Perturbation**: Handles q≠0 perturbations
2. **q=0 Specialization**: Computes dielectric tensor, Born charges, LO-TO splitting
3. **Metal System Support**: Handles smearing, Fermi level correction
4. **Phonon Calculation**: Assembles dynamical matrix, computes frequencies
5. **Symmetry Support**: Uses irreducible representations for efficiency

## Usage

```cpp
// In ESolver
ModuleDFPT::DFPT_PW dfpt;
dfpt.init(ucell, psi, nelec, ecutwfc);
dfpt.set_qmesh(4, 4, 4);
dfpt.set_conv_thr(1e-8);
dfpt.run();

// Get results
std::vector<double> freq = dfpt.get_phonon_freq(q_idx);
ModuleBase::matrix eps = dfpt.get_dielectric_tensor();
```

## Development Status

- **Phase**: Design phase
- **Author**: Mohan Chen
- **Date**: 2026-05-18
- **Status**: Not yet production-ready

## Future Work

1. Implement core Sternheimer solver
2. Add proper error handling
3. Complete unit tests
4. Optimize parallelization
5. Add LCAO support (separate module)
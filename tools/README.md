# ABACUS Tools

This directory contains various auxiliary tools for ABACUS calculations.

## Directory Structure

```
tools/
├── README.md
│
├── 01_NAO_generation/                # Numerical atomic orbital generation tools
│   ├── SIAB/                         # Simulated Annealing method (C++)
│   ├── pytorch/                      # PyTorch gradient method V1
│   ├── pytorch_dpsi/                 # PyTorch gradient method V2 (with dpsi)
│   ├── pytorch_gradient_source/      # Original PyTorch gradient implementation
│   ├── lcao_bash/                    # LCAO basis set bash tools
│   ├── abfs_bash/                    # ABFS basis set bash tools
│   ├── qo/                           # Quasiatomic orbital (QO) generation
│   ├── Generate_Orbital_AllInOne.sh  # Main orbital generation script
│   └── examples/                     # Usage examples
│
├── 02_postprocessing/                # Post-processing and visualization tools
│   ├── rt-tddft-tools/               # Real-time TDDFT analysis
│   ├── stm/                          # STM image generation
│   ├── average_pot/                  # Average electrostatic potential
│   ├── selective_dynamics/           # ABACUS + Phonopy phonon calculation
│   └── plot-tools/                   # Band structure, DOS, dipole and absorption
│
├── 03_code_analysis/                 # Source code analysis tools
│   └── generate_include_analysis.py
│
└── 04_windows_installation/          # Windows one-click installer via WSL2
    ├── install-abacus.bat
    ├── uninstall-abacus.bat
    ├── provision.sh
    └── README.md
```

## Quick Start

### Generate Numerical Atomic Orbitals

```bash
cd 01_NAO_generation
./Generate_Orbital_AllInOne.sh ORBITAL_INPUT
```

### RT-TDDFT Post-processing

```bash
cd 02_postprocessing/rt-tddft-tools
python plot_absorption.py --help
```

### STM Image Generation

```bash
cd 02_postprocessing/stm
python stm.py --help
```

### Average Electrostatic Potential

```bash
cd 02_postprocessing/average_pot
python aveElecStatPot.py --help
```

### Windows Installation

Run as administrator:
```bash
04_windows_installation/install-abacus.bat
```

## Tool Descriptions

### 01_NAO_generation/
- **SIAB/** - Simulated Annealing method for NAO optimization (C++)
- **pytorch/** - PyTorch gradient method V1 for NAO optimization
- **pytorch_dpsi/** - PyTorch gradient method V2 with dpsi calculation
- **pytorch_gradient_source/** - Original PyTorch gradient implementation
- **lcao_bash/** - Bash tools for LCAO basis set generation
- **abfs_bash/** - Bash tools for ABFS basis set generation
- **qo/** - Quasiatomic orbital (QO) generation tool
- **examples/** - Example input files and test cases

### 02_postprocessing/
- **rt-tddft-tools/** - Tools for real-time time-dependent density functional theory analysis
- **stm/** - Tools for generating STM images from LDOS cube files
- **average_pot/** - Python script to calculate and plot average electrostatic potential
- **selective_dynamics/** - Tools for selective dynamics with ABACUS + Phonopy
- **plot-tools/** - Band structure, DOS, dipole and absorption plotting tools

### 03_code_analysis/
- **generate_include_analysis.py** - Header file dependency depth analysis tool

### 04_windows_installation/
- Windows one-click installer via WSL2 + conda-forge
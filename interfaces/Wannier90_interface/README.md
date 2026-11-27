# A Brief Introduction to ABACUS Wannier90 Interface
This package provides a user-friendly Python interface: abacusw90 to bridge **ABACUS** (Atomic-scale Simulation Package) with **Wannier90**. It automates the workflow of generating Maximally Localized Wannier Functions (MLWFs) and tight-binding models from ABACUS calculations.

## What is Wannier90?

Wannier90 is an open-source code that calculates maximally localized Wannier functions (MLWFs) from first-principles calculations. It is designed to:

- Generate maximally localized Wannier functions
- Calculate band structures and density of states
- Compute Berry phases and orbital magnetization
- Provide a basis for tight-binding models
- Support various first-principles calculation codes through interfaces

Wannier functions are particularly useful for:
- Electronic structure calculations
- Transport properties
- Spectroscopy calculations
- Model Hamiltonian construction

## ABACUS-Wannier90 Interface: abacusw90

The abacusw90 package allows ABACUS to generate the necessary files for Wannier90, including:

- `*.amn` files: Overlap matrix between Bloch functions and Wannier functions
- `*.mmn` files: Overlap matrix between Bloch functions at neighboring k-points
- `UNK*` files: Bloch wavefunctions

## Key features
- **Various Basis Sets**: Support for both plane wave (PW) and LCAO basis sets.
- **Automated Workflow**: Handles the core coupling pipeline (Steps 3-5 of the standard tutorial workflow).
- **Input Generation**: Automatically generates `wannier90.win`, `INPUT`, `KPT`, and `STRU` files.
- **Method Support**: Supports the recommended `wannier_method = 2` for efficient overlap matrix calculation.
- **Spin-Orbit Coupling**: Full support for SOC calculations (`nspin=4`, `lspinorb=1`).

## Examples

This directory contains three examples demonstrating different use cases of the ABACUS-Wannier90 interface:

### 1. 01_lcao
- **System**: Diamond (C)
- **Basis**: LCAO (Linear Combination of Atomic Orbitals)
- **Purpose**: Demonstrates Wannier90 calculation using LCAO basis set
- **Input Files**:
  - `INPUT-scf`: ABACUS input file for SCF calculation
  - `INPUT-nscf`: ABACUS input file for NSCF calculation
  - `KPT-scf`: k-point sampling file for SCF calculation
  - `KPT-nscf`: k-point sampling file for NSCF calculation
  - `STRU`: Crystal structure file for diamond
  - `diamond.win`: Wannier90 input file
  - `diamond.nnkp`: Wannier90 preprocessing file

### 2. 02_pw
- **System**: Diamond (C)
- **Basis**: Plane wave (PW)
- **Purpose**: Demonstrates Wannier90 calculation using plane wave basis set
- **Input Files**: Similar to 01_lcao, but configured for plane wave basis

### 3. 03_lcao_in_pw
- **System**: Diamond (C)
- **Basis**: LCAO in plane wave mode
- **Purpose**: Demonstrates Wannier90 calculation using LCAO basis set in plane wave mode
- **Input Files**: Similar to 01_lcao, but configured for LCAO in plane wave mode

# How to Use abacusw90
## Installation
```bash
pip install .
# Or for development
pip install -e .
```
## Workflow Scope
This interface automates the technical coupling steps between ABACUS and Wannier90. In the context of the standard tutorial workflow, it covers the following stages:
| Step | Description | Responsibility |
| :--- | :--- | :--- |
| **Automated** | **Step 1**: ABACUS SCF Calculation | **Interface Step 0** |
| **Prerequisite** | **Step 2**: Determine Energy Windows | User provides `dis_win` parameters |
| **Automated** | **Step 3**: Generate `wannier90.win` & Run `-pp` | **Interface Step 1** |
| **Automated** | **Step 4**: ABACUS NSCF (Interface Mode) | **Interface Step 2 & 3** |
| **Automated** | **Step 5**: Wannier90 Minimization | **Interface Step 4** |
| **Post-process** | **Step 6**: WannierTools Analysis | User (Downstream tool) |
## Quick Start
Here is an example of generating Wannier functions for Bi2Se3:
```python
from abacusw90 import ABACUSWannier90
# 1. Initialize
# Assumes 'scf_dir' contains results from Step 1 (CHG, HR files)
job = ABACUSWannier90(work_dir="./Bi2Se3_wannier", scf_dir="./Bi2Se3_scf")
# 2. Define Structure
lattice = [[-2.069, -3.583614, 0.0], [2.069, -3.583614, 0.0], [0.0, 2.389075, 9.546667]]
atoms = [
    {"name": "Bi", "pos": [0.399, 0.399, 0.697]},
    {"name": "Bi", "pos": [0.601, 0.601, 0.303]},
    # ... (other atoms)
]
job.set_structure(lattice, atoms)
# 3. Configure Wannier90
# Parameters usually determined in Step 2 (Band structure analysis)
job.set_wannier_parameters(
    num_wann=30,
    num_bands=100,
    projections=["Bi : pz; px; py", "Se : pz; px; py"],
    dis_win_min=3.0,
    dis_win_max=18.0,
    dis_froz_min=3.0,
    dis_froz_max=14.8,
    mp_grid=[4, 4, 4],
    kpath=[
        {"start_label": "G", "start_pos": [0,0,0], "end_label": "Z", "end_pos": [0,0,0.5]}
    ]
)
# 4. Configure ABACUS
job.set_abacus_parameters(ecutwfc=100, nbands=100, lspinorb=1)
# 5. Run Automation (Covers Tutorial Steps 3, 4, 5)
job.run()
```

## Detailed Workflow Steps
The `run()` method executes the following automated sequence:
1.  **Generate Inputs & Preprocess**: Write `wannier90.win` and execute `wannier90 -pp` to generate `.nnkp`.
2.  **Prepare ABACUS**: Parse `.nnkp` to generate ABACUS `KPT`, `INPUT`, and `STRU` files. Copy SCF charge densities.
3.  **Run ABACUS Interface**: Execute ABACUS in NSCF mode with `towannier90=1`. This generates `mmn`, `amn`, `eig` files.
4.  **Run Wannier90**: Execute `wannier90.x` to compute MLWFs and output `wannier90_hr.dat`.

## Requirements
- **ABACUS**: v3.0 or higher (with Wannier90 interface support).
- **Wannier90**: v3.0 or higher.
- **Python**: 3.8+

# Important Notes

- The k-point grid in the ABACUS NSCF calculation must match the one in the Wannier90 input file
- Set `wvfn_formatted = .true.` in the Wannier90 input file to ensure compatibility with ABACUS output
- For LCAO calculations, ensure that the orbital basis set is appropriate for the Wannier functions you want to generate

# Troubleshooting

- **Files not found**: Ensure ABACUS is generating `*.amn`, `*.mmn`, and `UNK*` files in the OUT.* directory
- **Wannier90 cannot read ABACUS output**: Check that `wvfn_formatted = .true.` is set in the Wannier90 input file
- **Convergence issues**: Ensure the SCF and NSCF calculations are properly converged

# References

- **Wannier90 website**: [http://www.wannier.org/](http://www.wannier.org/)
- **Wannier90 paper**: A. A. Mostofi et al., *Comput. Phys. Commun.* **185**, 2309 (2014)
- **ABACUS documentation**: Refer to the ABACUS user manual for more details on input parameters

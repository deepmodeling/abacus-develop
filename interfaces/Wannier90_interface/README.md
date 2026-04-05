# ABACUS Wannier90 Interface
This package provides a user-friendly Python interface to bridge **ABACUS** (Atomic-scale Simulation Package) with **Wannier90**. It automates the workflow of generating Maximally Localized Wannier Functions (MLWFs) and tight-binding models from ABACUS calculations.
## Features
- **Automated Workflow**: Handles the core coupling pipeline (Steps 3-5 of the standard tutorial workflow).
- **Input Generation**: Automatically generates `wannier90.win`, `INPUT`, `KPT`, and `STRU` files.
- **Method Support**: Supports the recommended `wannier_method = 2` for efficient overlap matrix calculation.
- **Spin-Orbit Coupling**: Full support for SOC calculations (`nspin=4`, `lspinorb=1`).
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
| **Prerequisite** | **Step 1**: ABACUS SCF Calculation | User provides `scf_dir` |
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
```


# Hefei-NAMD Interface

This directory contains the interface between ABACUS and Hefei-NAMD, an open-source, ab initio non-adiabatic molecular dynamics (NAMD) package designed for simulating ultrafast excited-state carrier dynamics in condensed-matter systems.

## What is Hefei-NAMD?

Hefei-NAMD is an open-source, ab initio non-adiabatic molecular dynamics (NAMD) package designed for simulating ultrafast excited-state carrier dynamics in condensed-matter systems (e.g., semiconductors, 2D materials, oxides, and interfaces) across real/momentum space, energy, and time domains. Developed primarily by Prof. Qijing Zheng (University of Science and Technology of China, USTC) and collaborators, it is hosted publicly on GitHub: [github.com/QijingZheng/Hefei-NAMD](https://github.com/QijingZheng/Hefei-NAMD).

### Core Methodology

- **Theoretical foundation**: Combines time-dependent density-functional theory (TDDFT) with the fewest-switches surface hopping (FSSH) algorithm and decoherence-induced surface hopping (DISH) to treat electron–phonon coupling and non-adiabatic transitions between electronic states.
- **Key quantities**: Computes non-adiabatic coupling (NAC) matrices and electron–phonon coupling (EPC) from first-principles wavefunctions.
- **Implementation**: Mainly written in Fortran (≈80%) with Python pre/post-processing tools; minimal external dependencies, easy to compile via make.

### Key Features

- **Interfaces to major DFT codes**: Works natively with VASP and ABACUS (via abacus-namd scripts) to read MD trajectories, wavefunctions, and band structures.
- **Real- and momentum-space dynamics**: Supports carrier relaxation, recombination, hot-carrier transport, and exciton dynamics.
- **Advanced capabilities**:
  - Spin–orbit coupling (SOC) (via NAMDwithSOC extension).
  - Light–matter interaction (NAMD-LMI) for photoexcitation and stimulated emission.
  - Finite-temperature MD and trajectory sampling.
- **Parallelization**: MPI support for high-throughput calculations.

### Workflow

1. **Geometry optimization** → ab initio MD → trajectory snapshots.
2. **DFT single-points (SCF)** → wavefunction output.
3. **NAMD simulation** → analysis of carrier lifetime, relaxation pathways, etc.

### Applications

- **Carrier relaxation/recombination**: e.g., in TiO₂, perovskites, 2D van der Waals heterostructures.
- **Hot-carrier dynamics and transport** in photovoltaics and photocatalysts.
- **Exciton dissociation and charge transfer** at interfaces.
- **Ultrafast spectroscopy** (transient absorption, photoluminescence) interpretation.

### Availability & Use

- **Open-source**: Free under academic license; source code, manual, tutorials, and examples are fully public.
- **Citation**: Users are requested to cite core publications (e.g., Zheng et al., WIREs Comput. Mol. Sci., 2021).

In short, Hefei-NAMD is a widely used, community-driven tool for first-principles excited-state dynamics in materials science, valued for its efficiency, flexibility, and compatibility with standard DFT workflows.

## ABACUS-Hefei-NAMD Interface

The ABACUS-Hefei-NAMD interface allows ABACUS to generate the necessary files for Hefei-NAMD calculations, including overlap matrices, wavefunctions, and Hamiltonian matrices.

### Key Parameters for ABACUS-Hefei-NAMD Interface:

- `cal_syns`: Set to 1 to calculate asynchronous overlap matrix
- `dmax`: Maximum displacement of all atoms in one step (in bohr) for calculating asynchronous overlap matrix
- `out_wfc_lcao`: Set to 1 to output wavefunction files
- `out_mat_hs`: Set to 1 to output Hamiltonian and overlap matrix files

### Examples

The `example01` directory contains a simple example demonstrating how to set up an ABACUS calculation that generates the necessary files for Hefei-NAMD.

### How to Use the Interface

1. Set up an ABACUS calculation with the appropriate parameters for Hefei-NAMD (see examples)
2. Run ABACUS to generate the necessary output files
3. Use the generated files as input for Hefei-NAMD calculations

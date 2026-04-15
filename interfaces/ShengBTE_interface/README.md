# ShengBTE Interface

This directory contains the interface between ABACUS and ShengBTE, an open-source code for calculating lattice thermal conductivity using the Boltzmann Transport Equation (BTE).

## What is ShengBTE?

ShengBTE is an open-source code that calculates lattice thermal conductivity from first principles using the Boltzmann Transport Equation (BTE). It is designed to:

- Calculate lattice thermal conductivity
- Compute phonon dispersion relations
- Calculate phonon-phonon scattering rates
- Analyze thermal transport properties
- Support various first-principles calculation codes through interfaces

ShengBTE is particularly useful for:
- Studying thermal transport in materials
- Designing materials with specific thermal properties
- Understanding phonon scattering mechanisms
- Predicting thermal conductivity of new materials

## ABACUS-ShengBTE Interface

The ABACUS-ShengBTE interface allows ABACUS to generate the necessary files for ShengBTE, including:

- Second-order interatomic force constants (IFCs)
- Third-order interatomic force constants (IFCs)
- Phonon dispersion relations
- Crystal structure information

### Key Features of the Interface

- **Support for both plane wave (PW) and LCAO basis sets**
- **Compatibility with ShengBTE's input file format**
- **Automated generation of force constants**
- **Seamless integration with ShengBTE's workflow**

## Examples

This directory contains two examples demonstrating different use cases of the ABACUS-ShengBTE interface:

### 1. 01_pw
- **System**: Silicon (Si)
- **Basis**: Plane wave (PW)
- **Purpose**: Demonstrates ShengBTE calculation using plane wave basis set
- **Input Files**:
  - ABACUS input files for force constant calculations
  - ShengBTE input files
  - Scripts for running the calculation

### 2. 02_lcao
- **System**: Silicon (Si)
- **Basis**: LCAO (Linear Combination of Atomic Orbitals)
- **Purpose**: Demonstrates ShengBTE calculation using LCAO basis set
- **Input Files**:
  - ABACUS input files for force constant calculations
  - ShengBTE input files
  - Scripts for running the calculation

## How to Use the Interface

### Prerequisites

1. **Install ShengBTE**: Follow the installation instructions on the [ShengBTE website](https://shengbte.org/)
2. **Set up ABACUS**: Ensure ABACUS is compiled with ShengBTE support
3. **Prepare input files**: Create ABACUS input files and ShengBTE input files

### Basic Workflow

1. **Run ABACUS calculations to generate force constants**:
   - Calculate second-order force constants
   - Calculate third-order force constants
   - This will generate files containing the force constant information

2. **Prepare ShengBTE input files**:
   - Create `control.f90` file with calculation parameters
   - Create `structure.f90` file with crystal structure information
   - Create `BTE.out` file with initial conditions

3. **Run ShengBTE**:
   - Use the force constants from ABACUS as input
   - This will calculate the lattice thermal conductivity

### Important Notes

- The accuracy of the thermal conductivity calculation depends on the quality of the force constants
- Ensure that the supercell size is sufficient for accurate force constant calculation
- Convergence tests should be performed for supercell size and k-point sampling

## Troubleshooting

- **Force constants not found**: Ensure ABACUS is generating the necessary force constant files
- **ShengBTE cannot read ABACUS output**: Check that the force constant files are in the correct format
- **Convergence issues**: Ensure the force constant calculations are properly converged

## References

- **ShengBTE website**: [https://shengbte.org/](https://shengbte.org/)
- **ShengBTE paper**: L. Cheng et al., *Phys. Rev. B* **84**, 214301 (2011)
- **ABACUS documentation**: Refer to the ABACUS user manual for more details on input parameters

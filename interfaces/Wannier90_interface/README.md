# Wannier90 Interface

This directory contains the interface between ABACUS and Wannier90, an open-source code for generating maximally localized Wannier functions (MLWFs) and using them for various electronic structure calculations.

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

## ABACUS-Wannier90 Interface

The ABACUS-Wannier90 interface allows ABACUS to generate the necessary files for Wannier90, including:

- `*.amn` files: Overlap matrix between Bloch functions and Wannier functions
- `*.mmn` files: Overlap matrix between Bloch functions at neighboring k-points
- `UNK*` files: Bloch wavefunctions

### Key Features of the Interface

- **Support for both plane wave (PW) and LCAO basis sets**
- **Compatibility with Wannier90's input file format**
- **Automated generation of Wannier90 input files**
- **Seamless integration with Wannier90's workflow**

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

## How to Use the Interface

### Prerequisites

1. **Install Wannier90**: Follow the installation instructions on the [Wannier90 website](http://www.wannier.org/)
2. **Set up ABACUS**: Ensure ABACUS is compiled with Wannier90 support
3. **Prepare input files**: Create ABACUS input files and Wannier90 input files

### Basic Workflow

1. **Prepare Wannier90 input file** (`diamond.win`)
2. **Run Wannier90 preprocessing**:
   ```bash
   wannier90 -pp diamond.win
   ```
   This will generate `diamond.nnkp` file

3. **Run ABACUS SCF calculation**:
   - No need for `diamond.nnkp` file in this step
   - This will generate the converged charge density

4. **Run ABACUS NSCF calculation**:
   - Include `diamond.nnkp` file in the calculation directory
   - Use a k-point grid similar to what's defined in `diamond.win`
   - This will generate `diamond.amn`, `diamond.mmn`, and `UNK*` files in the OUT.* directory

5. **Run Wannier90**:
   - Copy `diamond.amn`, `diamond.mmn`, and `UNK*` files to the Wannier90 directory
   - Ensure `wvfn_formatted = .true.` is set in `diamond.win`
   - Run:
     ```bash
     wannier90 diamond.win
     ```
   - This will generate the maximally localized Wannier functions

### Important Notes

- The k-point grid in the ABACUS NSCF calculation must match the one in the Wannier90 input file
- Set `wvfn_formatted = .true.` in the Wannier90 input file to ensure compatibility with ABACUS output
- For LCAO calculations, ensure that the orbital basis set is appropriate for the Wannier functions you want to generate

## Troubleshooting

- **Files not found**: Ensure ABACUS is generating `*.amn`, `*.mmn`, and `UNK*` files in the OUT.* directory
- **Wannier90 cannot read ABACUS output**: Check that `wvfn_formatted = .true.` is set in the Wannier90 input file
- **Convergence issues**: Ensure the SCF and NSCF calculations are properly converged

## References

- **Wannier90 website**: [http://www.wannier.org/](http://www.wannier.org/)
- **Wannier90 paper**: A. A. Mostofi et al., *Comput. Phys. Commun.* **185**, 2309 (2014)
- **ABACUS documentation**: Refer to the ABACUS user manual for more details on input parameters

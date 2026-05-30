# DPGEN Interface

This directory contains the interface between ABACUS and DPGEN (Deep Potential GENerator), an open-source tool for generating deep learning potentials for molecular dynamics simulations.

## What is DPGEN?

DPGEN is an open-source tool that automates the generation of deep learning potentials (DPLs) for molecular dynamics simulations. It is designed to:

- Generate training data through ab initio calculations
- Train deep learning models for potential energy surfaces
- Validate and test the trained models
- Provide a framework for high-throughput materials simulations

DPGEN is particularly useful for:
- Accelerating molecular dynamics simulations
- Studying large systems that are computationally expensive with ab initio methods
- Simulating long-time dynamics
- Exploring complex phase spaces

## ABACUS-DPGEN Interface

The ABACUS-DPGEN interface allows ABACUS to be used as the ab initio calculator for DPGEN. This enables DPGEN to generate training data using ABACUS's various basis sets (plane wave and LCAO) and functional choices.

### Key Features of the Interface

- **Support for both plane wave (PW) and LCAO basis sets**
- **Compatibility with DPGEN's input file format**
- **Automated generation of training data**
- **Seamless integration with DPGEN's workflow**

## Examples

This directory contains two examples demonstrating different use cases of the ABACUS-DPGEN interface:

### 1. autotest
- **Purpose**: Demonstrates automated testing of the ABACUS-DPGEN interface
- **Input Files**:
  - DPGEN input files
  - ABACUS input files
  - Scripts for running the calculation

### 2. init_and_run
- **Purpose**: Demonstrates the complete workflow of initializing and running DPGEN with ABACUS
- **Input Files**:
  - DPGEN input files
  - ABACUS input files
  - Scripts for running the calculation

## How to Use the Interface

### Prerequisites

1. **Install DPGEN**: Follow the installation instructions on the [DPGEN GitHub repository](https://github.com/deepmodeling/dpgen)
2. **Set up ABACUS**: Ensure ABACUS is compiled and accessible in your PATH
3. **Prepare input files**: Create ABACUS input files and DPGEN input files

### Basic Workflow

1. **Prepare DPGEN input files**:
   - Create `input.json` file with calculation parameters
   - Define the structure and configuration space
   - Set up the ABACUS calculator settings

2. **Run DPGEN**:
   - DPGEN will automatically generate initial configurations
   - It will run ABACUS calculations to generate training data
   - It will train deep learning models based on the data
   - It will validate and test the trained models

3. **Use the trained model**:
   - The trained deep learning potential can be used for molecular dynamics simulations
   - It can be integrated with molecular dynamics packages like LAMMPS

### Important Notes

- The quality of the deep learning potential depends on the quality and diversity of the training data
- Ensure that the ABACUS calculations are properly converged
- Convergence tests should be performed for basis set, k-point sampling, and other parameters

## Troubleshooting

- **ABACUS calculations failing**: Check that ABACUS is properly installed and accessible
- **DPGEN cannot read ABACUS output**: Ensure that ABACUS is generating the necessary output files
- **Training not converging**: Consider increasing the amount of training data or adjusting the model parameters

## References

- **DPGEN GitHub repository**: [https://github.com/deepmodeling/dpgen](https://github.com/deepmodeling/dpgen)
- **DPGEN paper**: L. Zhang et al., *npj Comput. Mater.* **6**, 1 (2020)
- **ABACUS documentation**: Refer to the ABACUS user manual for more details on input parameters

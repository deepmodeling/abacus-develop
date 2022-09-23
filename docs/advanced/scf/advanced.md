# SCF in Complex Environments

## Implicit Solvation Model

## External Electric Field

##  Compensating Charge

## Van-der-Waals Correction

Conventional DFT functionals often suffer from an inadequate treatment of long-range dispersion, or Van der Waals (VdW) interactions. In order to describe materials where VdW interactions are prominent, one simple and popular approach is to add a Lennard-Jones type term. The resulting VdW-corrected DFT has been proved to be a very effective method for description of both short-range chemical bonding and long-range dispersive interactions.

Currently ABACUS provides three Grimme DFT-D methods, including [D2](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20495), [D3(0)](https://aip.scitation.org/doi/10.1063/1.3382344) and [D3(BJ)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21759), to describe Van der Waals interactions. Among them, the D3 method has been implemented in ABACUS based on the
dftd3 [program](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3) written by Stefan Grimme, Stephan Ehrlich and Helge Krieg.

To use VdW-correction, users need to supply value to the `vdw_method` keyword in the `INPUT` file:

   - (Default) none: no VdW correction
   - d2: DFT-D2 method
   - d3_0: DFT-D3(0) method
   - d3_bj: DFT-D3(BJ) method

Furthermore, ABACUS also provides a [list of keywords](../input_files/input-main.md#vdw-correction) to control relevant parmeters used in calculating the VdW correction, such as the scale factor (s6) term. Recommended values of such parameters can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of the parameters in ABACUS are set to be the recommended values for PBE.

Examples of VdW-corrected DFT calculations are provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/vdw/si2). There are two input files, where `INPUT1` shows how to apply D2 correction with user-specified $C_6$ parameter, and `INPUT2` shows how to apply D3(BJ) correction with default VdW parameters.

To run any of the two cases, users may enter the directory, copy the corresponding input file to `INPUT`, and run ABACUS.

## Dipole Correction

# SCF in Complex Environments

## Implicit Solvation Model

Solid-liquid interfaces are ubiquitous in nature and frequently encountered and employed in materials simulation. The solvation effect should be taken into account in first-principles calculations of such systems so as to obtain accurate predictions.  

Implicit solvation model is a well-developed method to deal with solvation effects, which has been widely used in finite and periodic systems. This approach treats the solvent as a continuous medium instead of individual “explicit” solvent molecules, which means that the solute embedded in an implicit solvent, and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath. Compared to the “explicit” method, such implicit solvation model can provide qualitatively correct results with much less computational cost, which is particularly suited for large and complex systems. The implicit solvation model implemented in ABACUS follows the [methodology](https://aip.scitation.org/doi/10.1063/1.4865107) developed by Mathew, Sundararaman, Letchworth-Weaver, Arias, and Hennig in 2014. 

Input parameters that control the implicit solvation model are listed as follows with detailed explaination and recommended values provided on this [webpage](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#implicit-solvation-model):

```
INPUT_PARAMETERS
imp_sol                 1
eb_k                    80
tau                     0.000010798
sigma_k                 0.6
nc_k                    0.00037
```

Example of running DFT calculation with the implicit solvation model is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/implicit_solvation_model/Pt-slab).

## External Electric Field 

A saw-like potential simulating an electric field
can be added to the bare ionic potential, which is a simplified simulation to the field-effect measurements, in which the system is separated from the gate electrode by a dielectric such as silicon oxide.

Whether to apply the external field is controlled via the keyword `efield_flag` in `INPUT` (setting to 1 to turn on the field). Related keywords that control the external field are listed as follows with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#electric-field-and-dipole-correction):
```
INPUT_PARAMETERS
efield_flag        1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0.001
```

Example of running DFT calculation with added external electric field is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/electric_field/Pt-slab). 


## Dipole Correction
A dipole correction can be added to the bare ionic potential, which can compensate for the artificial dipole field within the context of a periodic supercell calculation. The dipole correction implemented in ABACUS follows the [methodology](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.59.12301) proposed by Bengtsson in 1999. This correction must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE. Note that the common input parameters shared between the external electric field and dipole correction, with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#electric-field-and-dipole-correction). The following keywords settings add dipole correction only without applying any external electric field:
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0
```

While The external electric field and dipole correction can also be added together to the bare ionic potential as follows: 
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0.001
```

Examples of running DFT calculations with dipole correction are provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dipole_correction/Pt-slab). There are two input files, where `INPUT1` considers only the dipole correction without no applied external field, while `INPUT2` considers the dipole correction under an applied external field.

To run any of the two cases, users may enter the directory, copy the corresponding input file to `INPUT`, and run ABACUS.


## Compensating Charge

Modeling a constant-potential electronchemcial surface reaction requires adjustment of electron numbers in a simulation cell. At the mean time, we need to maintain the supercell's neutrality due to the periodic boundary condition. A distribution of compensating charge thus needs to be implemented in the vacuum region of surface models when extra electrons are added/extracted from the system.

The compensating charge implemented in ABACUS follows the [methodology](http://dx.doi.org/10.1103/PhysRevB.89.245406) developed by Brumme, Calandra, and Mauri in 2014. Input parameters that control the compensating charge are listed as follows with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#gate-field-compensating-charge): 

```
INPUT_PARAMETERS
gate_field         1
efield_dir         2
zgate              0.5
block              1
block_down         0.45
block_up           0.55
block_height       0.1
```

Example of running DFT calculation with the compensating charge is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/compensating_charge/Pt-slab). 

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


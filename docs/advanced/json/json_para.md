# JSON Configuration Parameters Documentation

- [JSON Configuration Parameters Documentation](#json-configuration-parameters-documentation)
  - [Overview](#overview)
  - [General Information](#general-information)
  - [Input](#input)
  - [Init](#init)
  - [Output](#output)
  - [Final Structure](#final-structure)

## Overview

This JSON template provides input and output configurations for ABACUS. It contains parameters for the execution of the program and the output of results, primarily for recording computational processes and outcomes for post processing. 

## General Information

- `version` - [str] The version number of ABACUS.
- `commit` - [str] The commit hash of ABACUS code at the time of computation.
- `start_time` - [str] The start time of the computation.
- `end_time` - [str] The end time of the computation.
- `device` - [str] The name of the hardware device on which the computation was run.
- `omp_number` - [int] The number of OpenMP threads.
- `mpi_number` - [int] The number of MPI processes.
- `out_dir` - [str] The output directory, e.g., "OUT.ABACUS".
- `log_file` - [str] The name of the log file, e.g., "running_scf.log".
- `pseudo_dir` - [str] The directory where pseudopotential files are stored.
- `orbital_dir` - [str] The directory where atomic orbital files are stored.
- `stru_file` - [str] The name of the structure file.
- `kpt_file` - [str] The name of the k-point file.

## Input
- A dictionary of parameters and their values as defined by the user in the INPUT file. (This part of the content will not be output in the current version yet.)


## Init


- `Input` - Lists the value of all input parameters. (This part of the content will not be output in the current version yet.)
- `point_group` - [str] the Schoenflies name of the point group.
- `point_group_in_space` - [str] the Schoenflies name of the point group in the space group.
- `nkstot`, `nkstot_ibz` - [int] Total number of k-points and total number of irreducible k-points.
- `nelectron_each_type` - [object(str-int)] The number of electrons for each atom type, e.g., `{"C": 2, "H":1}`.
- `nelectron` - [int] Total number of electrons.
- `nband` - [int] Number of bands.
- `natom` - [int] Total number of atoms.
- `label` - [array(str)] An array of atomic labels.
- `element` - [array(object(str:str))] The element of each atom type.
- `cell` - [array(array(double))] The lattice vector. Unit in Angstrom.
- `coordinate` - [array(array(double))] The cartesian coordinates of each atom. Unit in Angstrom.
- `mag` - [array(double)] The magnetic moments for each atom. 
- `pp` - [object(str-str)] The pseudopotential file of each atom type.
- `orb` - [object(str-str)] The orbital file of each atom type.


## Output

An array of dicts, including information about each self-consistent field (SCF) step, such as energy, convergence, and configuration:

- `energy`, `e_fermi` - [double] The total energy and Fermi energy. Unit in eV.
- `force` -  [array(array(double))] The forces calculated on each atom. Unit in eV/Angstrom.
- `stress` - [array(array(double))] The stress tensor. Unit in Kbar.
- `cell` - [array(array(double))] The cell parameters. Unit in Angstrom.
- `coordinate` - [array(array(double))] The coordinates of the atoms in the box after the simulation.
- `total_mag` , `absolute_mag` , `mag` - [double] The total magnetic moment; total absolute magnetic moment; and a list of magnetic moments for each atom, respectively.
- `scf_converge` - [bool] A boolean indicating whether the scf optimization has converged.
- `scf` - [array(object(str:double)] A list of each scf step, each item contains:
  - `energy` - [double] The total energy. Unit in eV. 
  - `ediff` - [double] The energy difference between the current and previous step. Unit in eV.
  - `drho` - [double] The charge density difference between the current and previous step.
  - `time` - [double] The time used for the current step. Unit in seconds.

## Final Structure
Parameters regarding the final converged results and the optimized geometry:

- `energy` - [double] The final energy.
- `label` - [str] An array of atomic labels.
- `cell` - [array(array(double))] The resulting cell parameters.
- `coordinate` - [array(array(double))] The final atomic coordinates.
- `relax_converge` - [bool] A boolean indicating whether the geometry optimization has converged.
- `dos` - [array(array(array(double)))] The state energy, and the dimension is NSPIN\*NKPOINT\*NBAND. 
- `dos_weight` - [array(array(array(double)))] The weight of each state, and the dimension is same as `dos`.
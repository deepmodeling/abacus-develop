# JSON Output Documentation

- [JSON Output Documentation](#json-output-documentation)
  - [Overview](#overview)
  - [General Information](#general-information)
  - [Initialization Information](#initialization-information)
  - [Output](#output)

## Overview

When JSON support is enabled, ABACUS writes calculation metadata and results to `abacus.json` for post-processing.

The current top-level JSON members are `comment`, `init`, `output`, and `general_info`. Some fields are populated only when the corresponding calculation data are available.

## General Information

The `general_info` object records basic build and runtime metadata:

- `version` - [string] ABACUS version.
- `commit` - [string] Git commit information when available at build time.
- `device` - [string] Hardware device selected for the calculation.
- `mpi_num` - [int] Number of MPI processes.
- `omp_num` - [int] Number of OpenMP threads.
- `pseudo_dir` - [string] Pseudopotential directory.
- `orbital_dir` - [string] Numerical atomic orbital directory.
- `stru_file` - [string] Structure input file.
- `kpt_file` - [string] K-point input file.
- `start_time` - [string] Calculation start time.
- `end_time` - [string] Time at which the JSON output is finalized.

## Initialization Information

The top-level `comment` describes the default units used by the JSON output. The `init` object records the initial structure and calculation settings. Depending on the calculation path, it can contain:

- `element` - [object(string:string)] Element/pseudopotential element information keyed by atom label.
- `orb` - [object(string:string/null)] Numerical orbital file for each atom type; `null` when no orbital file is used.
- `pp` - [object(string:string)] Pseudopotential file for each atom type.
- `coordinate` - [array(array(double))] Initial Cartesian coordinates in Angstrom.
- `mag` - [array(double)] Initial magnetic moment for each atom.
- `label` - [array(string)] Atomic labels.
- `cell` - [array(array(double))] Initial lattice vectors in Angstrom.
- `point_group` - [string] Schoenflies name of the point group.
- `point_group_in_space` - [string] Schoenflies name of the point group in the space group.
- `natom` - [int] Total number of atoms.
- `nband` - [int] Number of bands.
- `natom_each_type` - [object(string:int)] Number of atoms of each type.
- `nelectron_each_type` - [object(string:double)] Number of valence electrons for each atom type.
- `nelectron` - [int] Total number of electrons.
- `ecutwfc` - [double] Wavefunction energy cutoff.
- `ecutwfc_unit` - [string] Unit of `ecutwfc`, currently `Ry`.
- `smearing_method` - [string] Smearing method.
- `smearing_sigma` - [double] Smearing width.
- `smearing_sigma_unit` - [string] Unit of `smearing_sigma`, currently `Ry`.
- `kmesh_type` - [string] K-point mesh type.
- `kspacing` - [array(double)] K-point spacing parameters.
- `koffset` - [array(double)] K-point mesh offsets.
- `nkstot` - [int] Total number of k-points, when available on the calculation path.

## Output

`output` is an array. Each element represents one calculation/ionic-step output record. Fields are filled as the corresponding results become available:

- `energy` - [double/null] Total energy in eV.
- `e_fermi` - [double/null] Fermi energy in eV.
- `scf_converge` - [bool/null] Whether the SCF calculation converged.
- `force` - [array(array(double))/null] Atomic forces in eV/Angstrom when force calculation is enabled.
- `stress` - [array(array(double))/null] Stress tensor in kbar when stress calculation is enabled.
- `coordinate` - [array(array(double))] Cartesian coordinates in Angstrom.
- `mag` - [array(double)] Magnetic moment for each atom.
- `cell` - [array(array(double))] Lattice vectors in Angstrom.
- `total_mag` - [double] Total magnetic moment when available.
- `absolute_mag` - [double] Absolute magnetic moment when available.
- `scf` - [array(object)] SCF iteration history. Each entry contains:
  - `energy` - [double] Total energy in eV.
  - `ediff` - [double] Energy change from the previous SCF step in eV.
  - `drho` - [double] Charge-density difference.
  - `time` - [double] Time used by the SCF step in seconds.

JSON numbers are intended to be consumed as numeric values. Their textual representation (for example, decimal versus scientific notation) is not part of the output schema.

# Developer Guide

Abacuslite has the following file structure:
```
.
├── __init__.py
├── core.py
├── io
│   ├── __init__.py
│   ├── generalio.py
│   ├── latestio.py
│   ├── legacyio.py
│   └── testfiles
├── utils
│   ├── __init__.py
│   └── ksampling.py
└── xtest.sh
```

## core.py

This file contains the implementation of the Atomic Simulation Environment (ASE) calculator.

## io

This directory contains the input/output (I/O) functions for extracting information from ABACUS dumped files. A summary is here:

### Long-term-support (LTS) version

|Item|Regular expression|Example|
|----|----|----|
|esolver_type|'The esolver type has been set to : (\S+)'|`The esolver type has been set to : ksdft_pw`|
|nspin|'nspin\s+=\s+(\d+)'|`nspin = 4`|
|number of bands|'NBANDS\s+=\s+(\d+)'|`NBANDS = 40`|
|number of atoms|'TOTAL ATOM NUMBER\s*=\s*(\d+)'|`TOTAL ATOM NUMBER = 2`|
|lattice constant|'lattice constant \(Angstrom\)\s*=\s*(-?\d+(\.\d+)?)'|`lattice constant (Angstrom) = 0.529177`|
|lattice vectors|'^Lattice vectors: \(Cartesian coordinate: in unit of a_0\)$'|` Lattice vectors: (Cartesian coordinate: in unit of a_0)`|
|coordinate system|'^(DIRECT\|CARTESIAN) COORDINATES'|`DIRECT COORDINATES`|
|atomic positions|'^tau[c\|d]_([A-Z][a-z]?)\d+\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)\s+(-?\d+(\.\d+)?)*'|`taud_As1            0.2500000000        0.2500000000        0.2500000000  1.7321        0.0000000000        0.0000000000        0.0000000000`|
|eigenvalues|'\d+/\d+\s+kpoint\s+\(Cartesian\)\s+=\s+(-?\d(\.\d+)?(e-\d+)?)\s+(-?\d(\.\d+)?(e-\d+)?)\s+(-?\d(\.\d+)?(e-\d\|+)?)\s+\(\d+\s+pws\)'|` 1/1 kpoint (Cartesian) = 0.0000 0.0000 0.0000 (230 pws)`|
|atomic forces|'\s\*TOTAL\-FORCE\s\*\(eV\s*/Angstrom\)'|` TOTAL-FORCE (eV/Angstrom)     `|
|total stress|'\s\*TOTAL\-STRESS\s\*\(KBAR\)'|` TOTAL-STRESS (KBAR)          `|
|kpoints and weights|'\s\*(IBZ\|KPOINTS)\s+(DIRECT\|CARTESIAN)_X\s+(DIRECT\|CARTESIAN)_Y\s+(DIRECT\|CARTESIAN)_Z\s+WEIGHT(\s+ibz2bz)?'|` KPOINTS    DIRECT_X    DIRECT_Y    DIRECT_Z  WEIGHT`|
|energy|'\s*ENERGY\s+Rydberg\s+eV'|`     Energy           Rydberg                 eV      `|
|total magnetism|'\s\*Total\sMagnetism\s\(uB\)(\s+x\s+y\s+z)?\s\*'|` Total Magnetism (uB)    `|   

### Latest version

|Item|Regular expression|Example|Notes|
|----|----|----|----|
|esolver_type|'#ENERGY SOLVER#\s+(\S+)'|` #ENERGY SOLVER# ksdft_pw`||
|nspin|'nspin\s+=\s+(\d+)'|`nspin = 4`||
|number of bands|||has been removed|
|number of atoms|'TOTAL ATOM NUMBER\s*=\s*(\d+)'|`TOTAL ATOM NUMBER = 2`||
|lattice constant|'lattice constant \(Angstrom\)\s*=\s*(-?\d+(\.\d+)?)'|`lattice constant (Angstrom) = 0.529177`||
|lattice vectors|'^Lattice vectors: \(Cartesian coordinate: in unit of a_0\)$'|` Lattice vectors: (Cartesian coordinate: in unit of a_0)`||
|coordinate system|'^(DIRECT\|CARTESIAN) COORDINATES'|`DIRECT COORDINATES`||
|atomic positions|'atom\s+x\s+y\s+z\s+mag‘|` atom                  x                  y                  z     mag`|"tauc/taud" suffix has been removed, therefore the way to read the coordinates changes to find the head of the table, then read the following `number of atoms` lines. On the other hand, if the `calculation` is set to `md`, ABACUS will not dump the atomic positions to running log anymore, instead, will read from the MD_dump file|
|eigenvalues|'spin=(\d)\s+k-point=(\d+)/(\d+)\s+Cartesian=\s*(-?\d(\.\d+)?(e-\d+)?)\s+(-?\d(\.\d+)?(e-\d+)?)\s+(-?\d(\.\d\|+)?(e-\d+)?)\s+\(\d+\s+plane wave\)'|  spin=1 k-point=1/1 Cartesian=0.0000000 0.0000000 0.0000000 (1837 plane wave)|eigenvalues information has been removed from the running log, the file istate.info is renamed as eig_occ.txt, where the eigenvalues are read|
|atomic forces|'#\s\*TOTAL\-FORCE\s*\(eV\s*/Angstrom\)\s\*#'|`  #TOTAL-FORCE (eV/Angstrom)#`||
|total stress|'#\s\*TOTAL\-STRESS\s*\(kbar\)\s\*#'|` #TOTAL-STRESS (kbar)#    `||
|kpoints and weights|'\s*(IBZ\|KPOINTS)\s+(DIRECT\|CARTESIAN)_X\s+(DIRECT\|CARTESIAN)_Y\s+(DIRECT\|CARTESIAN)_Z\s+WEIGHT(\s+ibz2bz)?'|` KPOINTS    DIRECT_X    DIRECT_Y    DIRECT_Z  WEIGHT`||
|energy|'\s*ENERGY\s+Rydberg\s+eV'|`     Energy           Rydberg                 eV      `||
|total magnetism|'\s\*Total\sMagnetism\s\(uB\)(\s+x\s+y\s+z)?\s\*'|` Total Magnetism (uB)    `||

Please look at detailed implementations in the following files.

### generalio.py

This module contains the general I/O functions for ABACUS.

### latestio.py

This module contains the I/O functions for the latest version of ABACUS.

### legacyio.py

This module contains the I/O functions for the Long-Term-Support (LTS) version of ABACUS.

## utils

This directory contains the utility modules for Abacuslite.

### ksampling.py

This module contains the wrapper of k-point sampling functions and helper functions.

## xtest.sh

This script can run all the unittests that programmed in all Python source files.
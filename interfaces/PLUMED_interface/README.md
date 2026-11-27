# ABACUS + PLUMED interface

[PLUMED](https://www.plumed.org) is a library for enhanced-sampling methods
(metadynamics, umbrella sampling, ...). During an ABACUS molecular-dynamics run
it can compute collective variables and apply biasing potentials, controlled by
the `plumed` and `plumed_file` INPUT keywords.

## Building ABACUS with PLUMED

Install PLUMED 2.x first, for example:

    ./configure --prefix=/path/to/plumed && make -j && make install

Then configure ABACUS with `-DENABLE_PLUMED=ON`:

    cmake -B build -DENABLE_PLUMED=ON ...
    cmake --build build -j

The PLUMED installation is located with `pkg-config` (through the `plumed.pc`
file shipped by PLUMED); if PLUMED is installed in a non-standard prefix, add
`-DPLUMED_ROOT=/path/to/plumed` so that it can be found without pkg-config.
When ABACUS is built without `-DENABLE_PLUMED=ON`, setting `plumed` in INPUT
stops the run with an explicit error message.

## Usage

In the `INPUT` file:

    plumed          1            # enable the PLUMED interface
    plumed_file     plumed.dat   # PLUMED input file (default: plumed.dat)

The PLUMED input file defines the collective variables and the biasing
potential, e.g. for a dihedral angle:

    UNITS LENGTH=A TIME=fs
    t: TORSION ATOMS=3,1,2,7
    metad: METAD ARG=t PACE=50 HEIGHT=2.0 SIGMA=0.15 FILE=HILLS BIASFACTOR=20 TEMP=600
    PRINT ARG=t,metad.bias FILE=COLVAR STRIDE=1

ABACUS hands the atomic positions (Bohr), masses (amu), cell, potential energy
and virial over to PLUMED after every MD step, and PLUMED adds the biasing
forces in place.

Note: the interface currently supports a single MPI rank only, because PLUMED
must see the whole system while ABACUS distributes the atoms of an MD run over
the ranks. OpenMP threads inside that single rank are fine, e.g.
`OMP_NUM_THREADS=32 abacus > log`.

## Example

`example01` runs a well-tempered metadynamics of the ethane (C2H6) H-C-C-H
torsion at 600 K (LCAO). Running ABACUS in that directory produces

  - `COLVAR`: time, the torsion angle and the metadynamics bias at every MD step;
  - `HILLS`: the Gaussian hills deposited by metadynamics.

The free-energy surface can be reconstructed with the tools shipped with PLUMED:

    plumed sum_hills --hills HILLS --mintozero --outfile fes.dat

`fes.dat` contains the torsion angle (rad) in the first column and the free
energy (kJ/mol) in the second one.

## References

  - PLUMED manual: <https://www.plumed.org/doc>
  - ABACUS `plumed` / `plumed_file` keywords:
    <https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html>

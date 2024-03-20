# CANDELA

[CANDELA](https://github.com/MCresearch/Candela) is short for Collection of ANalysis DEsigned for Large-scale Atomic simulations. 
It is developed by MCresearch to conduct analyses on MD trajectory in different formats. 
Right now the program only supports analysis of pair distribution function (PDF), static structure factor (SSF) and mean square displacement (MSD).
The minimum supported version of ABACUS is 3.2.0.

## Requirements for using CANDELA 
For Detailed usage of CANDELA, please refer to the [official document](https://candela-docs.readthedocs.io/en/latest/).
There are two things which need special attention in using CANDELA with ABACUS.
First, the input file of CANDELA only takes the name of `INPUT`, the same as ABACUS input file, so you should not run CANDELA in the same folder where you run ABACUS.
Second, to use CANDELA to postprocess ABACUS MD trajectory, the following parameters have to be specified in the `INPUT` file of CANDELA in addition to other required parameters:
1. `geo_in_type` has to be set to `ABACUS`;
2. `msd_dt` has to be specified in unit of picosecond, especially in the case of `calculation` = `msd`;
3. `geo_directory` has to be set to the path to the `MD_dump` file in the `OUT.xxx` folder.
As a result, a CANDELA `INPUT` file for calculating PDF from ABACUS should be something like this:
```
calculation  pdf # Pair Distribution Function.
system Al
geo_in_type  LAMMPS
geo_directory ../geo/Al64.dump
geo_1        0
geo_2        20
geo_interval 1
geo_ignore   4

geo_out      pdf.txt # output pdf name.

ntype        1        # number of different types of atoms.
natom        64   # total number of atoms.
natom1       64
rcut         6
dr           0.01     # delta r in real space

struf_dgx   0.05
struf_ng    200
```

More examples of CANDELA `INPUT` with ABACUS can be found in the [test](https://github.com/MCresearch/Candela/tree/main/test) directory.
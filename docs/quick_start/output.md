# Brief Introduction of the Output Files

The following files are the central output files for ABACUS. After executing the program, you can obtain a log file containing the screen output, more detailed outputs are stored in the working directory `OUT.suffix` (Default one is `OUT.ABACUS`). Here we give some simple descriptions.

## *INPUT*

Different from `INPUT` given by the users, `OUT.suffix/INPUT` contains all parameters in ABACUS.

> **Note:** `OUT.suffix/INPUT` contain the initial default of ABACUS instead of the real parameters used in calculations. If you want to figure out the real parameters used in calculations, you can open `OUT.suffix/runing_scf.log` and research corresponding parameter you are interested.

For a complete list of input parameters, please consult this [instruction](../advanced/input_files/input-main.md).

## *running_scf.log*

`running_scf.log` contains information on nearly all function calls made during the execution of ABACUS.

## *KPT*

This file contains the information of all generated k-points, as well as the list of k-points actually used for calculations after considering symmetry.

## *istate.info*

This file includes the energy levels computed for all k-points. From left to right, the columns represent: energy level index, eigenenergy, and occupancy number.

Below is an example `istate.info`:

```
BAND               Energy(ev)               Occupation                Kpoint = 1                        (0 0 0)
      1                 -5.33892                  0.03125
      2                  6.68535                0.0312006
      3                  6.68535                0.0312006
      4                  6.68535                0.0312006
      5                  9.41058                        0
```

## *STRU_SIMPLE.cif*

ABACUS generates a `.cif` format structure file based on the input file `STRU`, facilitating users to visualize with commonly used software. `STRU_READIN_ADJUST.cif` is the structure after considering symmetry.

## *warning.log*

The file contains all the warning messages generated during the ABACUS run.
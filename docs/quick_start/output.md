# Brief Introduction of the Output Files

The following files are the central output files for ABACUS. After executing the program, you can obtain a log file containing the screen output, more detailed outputs are stored in the working directory `OUT.suffix` (Default one is `OUT.ABACUS`). Here we give some simple descriptions.

## *INPUT*

Different from `INPUT` given by the users, `OUT.suffix/INPUT` contains all parameters in ABACUS.

> **Note:** `OUT.suffix/INPUT` contains the **actual parameters used in the calculation**, including:
> 1. **User-specified parameters** (explicitly defined in your input file or command-line arguments, overriding default parameters).
> 2. **System default parameters** (automatically applied when not explicitly provided by the user).


This file ensures calculations can be fully reproduced, even if default values change in future ABACUS versions.
Also notice that in rare cases, a small number of parameters may be dynamically reset to appropriate values during runtime.

For a complete list of input parameters, please consult this [instruction](../advanced/input_files/input-main.md).

## *running_scf.log*

`running_scf.log` contains information on nearly all function calls made during the execution of ABACUS.

## *KPT.info*

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

## *STRU.cif*

ABACUS generates a `.cif` format structure file based on the input file `STRU`, facilitating users to visualize with commonly used software.

## *warning.log*

The file contains all the warning messages generated during the ABACUS run.

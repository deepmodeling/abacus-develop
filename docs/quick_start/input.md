# Brief Introduction of the Input Files

The following files are the central input files for ABACUS. Before executing the program, please make sure these files are prepared and stored in the working directory. Here we give some simple descriptions XXX. For more details, users should consult the Advanced session.

## *INPUT*

The `INPUT` file contains parameters that control the type of calculation as well as a variety of settings.

Below is an example `INPUT` file with some of the most important parameters that need to be set:

```
INPUT_PARAMETERS
#Parameters (General)
ntype 1
nbands 4

#Parameters (Accuracy)
ecutwfc 60
```

The parameter list always starts with key word `INPUT_PARAMETERS`. Any content before `INPUT_PARAMETERS` will be ignored.

Any line starting with `#` or `/` will also be ignored.

Each parameter value is provided by specifying the name of the input variable
and then putting the value after the name, separated by one or more blank characters(space or tab). The following characters (≤ 150) in the same line will be neglected.

Depending on the input variable, the value may be an integer, a real number or a string. The parameters can be given in any order, but only one parameter should be given per line.

Furthermore, if a given parameter name appeared more than once in the input file, only the last value will be taken.

> Note: if a parameter name is not recognized by the program, the program will stop with an error message.

In the above example, the meanings of the parameters are:

- `ntype` : how many types of elements in the unit cell
- `nbands` : the number of bands to be calculated
- `ecutwfc` : the plane-wave energy cutoff for the wave function expansion (UNIT: Rydberg)    

For a complete list of input parameters, please consult this [instruction](docs/input-main.md)

> Note: Users cannot change the filename “INPUT” to other names.

## *STRU*

The structure file contains structural information about the system, e.g., lattice constant, lattice vectors, and positions of the atoms within a unit cell. The positions can be given either in direct or Cartesian coordinates. 

An example of the `STRU` file is given as follows :

XXXXXXX

> Note : users may choose a different name for their structure file using the keyword XXXXX

For a more detailed description of STRU file, please consult XXXXX

## *KPT*

This file contains information of the kpoint grid setting for the Brillouin zone sampling.
    
An example of the `KPT` file is given by XXXXXXX

> Note : users may choose a different name for their k-point file using keyword XXXXX


For a more detailed description, please consult XXXXX

- The pseudopotential files

    Norm-conserving pseudopotentials are used in ABACUS, in the UPF file format.The filename of each element’s pseudopotential needs to be specified in the STRU file, together with the directory of the pseudopotential files unless they are already present in the working directory.

    More information on pseudopotentials is given [here](docs/features.md#pseudopotentials).

- The numerical orbital files

    This part is only required in LCAO calculations.
    The filename for each element’s numerical orbital basis needs to be specified in the STRU file, together with the directory of the orbital files unless they are already present in the working directory.
    ABACUS provides numerical atomic basis sets of different accuracy levels for most elements commonly used. Users can download these basis sets from the [website](http://abacus.ustc.edu.cn/pseudo/list.htm). Moreover, users can generate numerical atomic orbitals by themselves, and the procedure is provided in this [short introduction](docs/generate-basis.md).
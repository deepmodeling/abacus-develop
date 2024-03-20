<!--
 * @Date: 2021-08-21 21:58:06
 * @LastEditors: jiyuyang
 * @LastEditTime: 2022-01-03 17:21:08
 * @Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
-->

# Plotting tool for ABACUS

-   Band Structure and Fat Bands
-   DOS and PDOS
-   Dipole and Absorption

## Requirements

-   [NumPy](https://numpy.org/)
-   [Matplotlib](https://matplotlib.org/)
-   [lxml](https://lxml.de/)
-   [setuptools](https://setuptools.pypa.io/en/latest/index.html)

## Installation

Use the following command to install this tool:

```shell
python setup.py install
```

## Usage

There are two ways to use this tool:

1. Specify parameters in `band.py` (`dos.py`, `dipole.py`) directly, and then `python **.py`. And you can also import module in your own script e.g. `from abacus_plot.band import Band`. (Recommend)
2. Command-line tools are also supported in this tool. In this way, you need prepare an input file and execute some commands (see below). You can use `abacus-plot -h` to check command-line information

### Band Structure

First, prepare a file named 'config.json' in json format:

```json
{
    "bandfile": "BANDS_1.dat",
    "efermi": 6.585653952007503,
    "energy_range": [-1.5, 6],
    "kptfile": "KPT"
}
```

|    Property    |                                    Type                                     |                                                                                                Note                                                                                                |
| :------------: | :-------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|   _bandfile_   |                            `str` or `List[str]`                             |                                         Bands data file output from ABACUS. 'BANDS*\*.dat' for band structure and 'PBANDS*\*' for projected band structure                                         |
|    _efermi_    |                          `float` or `List[float]`                           |                                                                                         Fermi level in eV                                                                                          |
| _energy_range_ |                                   `list`                                    |                                                                                       Range of energy in eV                                                                                        |
|    _shift_     |                                   `bool`                                    |                                              If set `'true'`, it will evaluate band gap (only for semiconductors and insulators). Default: `'false'`                                               |
|    _index_     | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` |                     Extract PBANDS of each atom from 'PBANDS\_\*' file e.g. [index_0, index_1...] or {index_0:[l_0, l_1, ...], ...}, [index_0:{l_0:[m_0, m_1, ...], ...}, ...]                     |
|  _atom_index_  | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` | Extract PBANDS of each atom with same atom*index from 'PBANDS*\*' file e.g. [atom_index_0, atom_index_1...] or {atom_index_0:[l_0, l_1, ...], ...}, [atom_index_0:{l_0:[m_0, m_1, ...], ...}, ...] |
|   _species_    | `List[str]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` |              Extract PBANDS of each atom with same species from 'PBANDS\_\*' file e.g. [elem_0, elem_1...] or {elem_0:[l_0, l_1, ...], ...}, [elem_0:{l_0:[m_0, m_1, ...], ...}, ...]              |
|   _kptfile_    |                                    `str`                                    |                                                                          K-points file with `Line` mode in ABACUS format                                                                           |
|    _label_     |                            `str` or `List[str]`                             |                                                                                      Label of band structure                                                                                       |
|    _color_     |                            `str` or `List[str]`                             |                                                                                      Color of band structure                                                                                       |

If you want to plot band structure of `nspin=2` or compare two band structure on same k-path, you can set _bandfile_, _efermi_, _label_, _color_ in `list` type.

The _kptfile_ should be as follows, and notes after `#` will be set as k-points label automatically.

```shell
K_POINTS
7
Line
0.5     0.5     0.5     20      #L
0.0     0.0     0.0     20      #G
0.5     0.0     0.5     10      #X
0.5     0.25    0.75    15      #W
0.5     0.5     0.5     10      #L
0.375   0.375   0.75    20      #K
0.0     0.0     0.0     1       #G
```

Then, use the following command to plot band structure:

```shell
abacus-plot -b
```

Then, the following command will plot projected band structure and figures output to directory `PBAND*_FIG` ($*=1$ for $nspin=1$, $*=2$ for $nspin=2$):

```shell
abacus-plot -d -p
```

Then, the following command will output parsed partial DOS to directory `PBAND*_FILE` ($*=1$ for $nspin=1$, $*=2$ for $nspin=2$):

```shell
abacus-plot -d -o
```

### DOS

First, prepare a file named 'config.json' in json format:

```json
{
    "pdosfile": "PDOS",
    "efermi": 6.585653952007503,
    "energy_range": [-5, 7],
    "dos_range": [0, 5],
    "species": {
        "C": {
            "0": [0],
            "1": [0, 1]
        },
        "Si": [0, 1]
    },
    "pdosfig": "pdos.png"
}
```

If you only want to plot total DOS, you can modify `pdosfile` to `tdosfile` and do not set `species` and `pdosfig`.
| Property | Type | Note |
| :------------: | :-------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| _tdosfile_ | `str` | Total DOS data file output from ABACUS |
| _pdosfile_ | `str` | Partial DOS data file output from ABACUS in xml format |
| _efermi_ | `float` | Fermi level in eV |
| _energy_range_ | `list` | Range of energy in eV |
| _shift_ | `bool` | If set `'true'`, it will evaluate band gap. Default: `'false'` |
| _index_ | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` | Extract PDOS of each atom e.g. [index_0, index_1...] or {index*0:[l_0, l_1, ...], ...}, [index_0:{l_0:[m_0, m_1, ...], ...}, ...] |
| \_atom_index* | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` | Extract PDOS of each atom with same atom*index e.g. [atom_index_0, atom_index_1...] or {atom_index_0:[l_0, l_1, ...], ...}, [atom_index_0:{l_0:[m_0, m_1, ...], ...}, ...] |
| \_species* | `List[str]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` | Extract PDOS of each atom with same species e.g. [elem_0, elem_1...] or {elem*0:[l_0, l_1, ...], ...}, [elem_0:{l_0:[m_0, m_1, ...], ...}, ...] |
| \_tdosfig* | `str` | Output picture of total DOS |
| _pdosfig_ | `str` | Output picture of partial DOS |

Then, the following command will plot total DOS:

```shell
abacus-plot -d
```

Then, the following command will plot partial DOS:

```shell
abacus-plot -d -p
```

Then, the following command will output parsed partial DOS to directory `PDOS_FILE`:

```shell
abacus-plot -d -o
```

### Dipole and Absorption

```python

from abacus_plot.dipole import Dipole
from abacus_plot.dipole import Absorption
import matplotlib.pyplot as plt

dipolefile = './SPIN1_DIPOLE'
dipole = Dipole(dipolefile, dt=0.0024)
Efile=[["efield_0.dat"],["efield_1.dat"],["efield_2.dat"]]
#Efile is a 2D list, Efile[i][j] is the jth Efield data file in ith direction
Abs = Absorption(dipolefile, Efile, dt=0.0024)

fig1, ax1 = plt.subplots()
fig1, ax1 = dipole.plot_dipole(fig1, ax1)
fig1.savefig('dipole.png')

fig2, ax2 = plt.subplots()
x_range = [0, 20]
fig2, ax2 = Abs.plot_abs(
    fig2, ax2, x_range=x_range, unit='eV')
fig2.savefig('abs.png')
```

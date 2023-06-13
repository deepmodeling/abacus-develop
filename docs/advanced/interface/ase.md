# ASE

## Introduction

[ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment) provides a set of Python tools for setting, running, and analysing atomic simulations. We have developed an ABACUS calculator ([ase-abacus](https://gitlab.com/1041176461/ase-abacus )) to be used together with the ASE tools, which exists as an external project with respect to ASE and is maintained by ABACUS developers.

## Installation

```bash
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
python3 setup.py install
```

## Environment variables

[ABACUS](http://abacus.ustc.edu.cn) supports two types of basis sets: PW, LCAO. The path of pseudopotential and numerical orbital files can be set throught the environment variables `ABACUS_PP_PATH` and `ABACUS_ORBITAL_PATH`, respectively, e.g.:

```bash
  PP=${HOME}/pseudopotentials
  ORB=${HOME}/orbitals
  export ABACUS_PP_PATH=${PP}
  export ABACUS_ORBITAL_PATH=${ORB}
```
 
For PW calculations, only `ABACUS_PP_PATH` is needed. For LCAO calculations, both `ABACUS_PP_PATH` and `ABACUS_ORBITAL_PATH` should be set.

## ABACUS Calculator

The default initialization command for the ABACUS calculator is

```python
from ase.calculators.abacus import Abacus
```

In order to run a calculation, you have to ensure that at least the following parameters are specified, either in the initialization or as environment variables:

|keyword         |description
|:---------------|:----------------------------------------------------------
|`pp`            |dict of pseudopotentials for involved elememts, <br> such as `pp={'Al':'Al_ONCV_PBE-1.0.upf',...}`.
|`pseudo_dir`    |directory where the pseudopotential are located, <br> Can also be specified with the `ABACUS_PP_PATH` <br> environment variable. Default: `pseudo_dir=./`.
|`basis`         |dict of orbital files for involved elememts, such as <br> `basis={'Al':'Al_gga_10au_100Ry_4s4p1d.orb'}`.<br> It must be set if you want to do LCAO <br> calculations. But for pw calculations, it can be omitted.
|`basis_dir`     |directory where the orbital files are located, <br> Can also be specified with the `ABACUS_ORBITAL_PATH`<br> environment variable. Default: `basis_dir=./`.
|`xc`            |which exchange-correlation functional is used.<br> An alternative way to set this parameter is via <br> seting `dft_functional` which is an ABACUS <br> parameter used to specify exchange-correlation <br> functional
|`kpts`          |a tuple (or list) of 3 integers `kpts=(int, int, int)`, <br>it is interpreted as the dimensions of a Monkhorst-Pack <br>  grid, when `kmode` is `Gamma` or `MP`. It is <br>  interpreted as k-points, when `kmode` is `Direct`,<br>  `Cartesian` or `Line`, and `knumber` should also<br>  be set in these modes to denote the number of k-points.<br>  Some other parameters for k-grid settings:<br>  including `koffset` and `kspacing`.

For more information on pseudopotentials and numerical orbitals, please visit [ABACUS]. The elaboration of input parameters can be found [here](../input_files/input-main.md).


The input parameters can be set like::
```python
  calc = Abacus(profile=profile, ntype=1, ecutwfc=50, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.01, basis_type='pw', ks_solver='cg', calculation='scf' pp=pp, basis=basis, kpts=kpts)
```

The command to run jobs can be set by specifying `AbacusProfile`::

```python
  from ase.calculators.abacus import AbacusProfile
  abacus = '/usr/local/bin/abacus'
  profile = AbacusProfile(argv=['mpirun','-n','2',abacus])
```

in which `abacus` sets the absolute path of the `abacus` executable.

## MD Analysis
After molecular dynamics calculations, the log file `running_md.log` can be read. If the 'STRU_MD_*' files are not continuous (e.g. 'STRU_MD_0', 'STRU_MD_5', 'STRU_MD_10'...), the index parameter of read should be as a slice object. For example, when using the command `read('running_md.log', index=slice(0, 15, 5), format='abacus-out')` to parse 'running_md.log', 'STRU_MD_0', 'STRU_MD_5' and 'STRU_MD_10' will be read.


## SPAP Analysis

[SPAP](https://github.com/chuanxun/StructurePrototypeAnalysisPackage) (Structure Prototype Analysis Package) is written by Dr. Chuanxun Su to analyze symmetry and compare similarity of large amount of atomic structures. The coordination characterization function (CCF) is used to 
measure structural similarity. An unique and advanced clustering method is developed to automatically classify structures into groups. 


If you use this program and method in your research, please read and cite the publication:

`Su C, Lv J, Li Q, Wang H, Zhang L, Wang Y, Ma Y. Construction of crystal structure prototype database: methods and applications. J Phys Condens Matter. 2017 Apr 26;29(16):165901.`

and you should install it first with command `pip install spap`.

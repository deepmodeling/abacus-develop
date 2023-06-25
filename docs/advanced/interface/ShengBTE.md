# ShengBTE
## Introduction

This tutorial aims to introduce the process of performing density functional theory calculations using ABACUS and calculating lattice thermal conductivity using the ShengBTE software. During the entire calculation process, the following external software are also used: 1) Phonopy for calculating second-order force constants, 2) ASE for converting atomic structures, 3) ShengBTE's thirdorder program for calculating third-order force constants, and 4) finally using ShengBTE to calculate the material's lattice thermal conductivity.

Here is the announcement of ShengBTE with ABACUS: [ShengBTE - The ABACUS DFT software can now be used with ShengBTE](https://www.shengbte.org/announcements/the-abacus-dft-software-can-now-be-used-with-shengbte)

Some external packages that need to be combined are mentioned above, and here it is recommended to read the relevant documentation and instructions of these packages:

ShengBTE：[https://bitbucket.org/sousaw/shengbte/src/master/](https://bitbucket.org/sousaw/shengbte/src/master/)

phonopy：[http://abacus.deepmodeling.com/en/latest/advanced/interface/phonopy.html](http://abacus.deepmodeling.com/en/latest/advanced/interface/phonopy.html)

ASE：[http://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html](http://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html)

thirdorder: [https://bitbucket.org/sousaw/thirdorder/src/master/](https://bitbucket.org/sousaw/thirdorder/src/master/)

## Prepare

The ABACUS software package provides an example of using ABACUS+ShengBTE to calculate the lattice thermal conductivity in the `examples/interface_ShengBTE` folder. The example includes two folders: `LCAO` (Linear Combination of Atomic Orbitals) which uses numerical atomic orbitals and `PW` (Plane wave) which uses plane wave basis vectors. Each folder contains three subfolders: `2nd`, `3rd`, and `shengbte`, which respectively store the relevant files for calculating second-order force constants using Phonopy (`2nd`), third-order force constants using the thirdorder program (`3rd`), and lattice thermal conductivity using ShengBTE (`shengbte`).

## How to use

Taking the `LCAO` folder as an example, we provide the test case of a diamond structure Si structure containing 2 atoms with the norm-conserving pseudopotential `Si_ONCV_PBE-1.0.upf` and the atomic orbital file `Si_gga_7au_100Ry_2s2p1d.orb` (GGA functional, 7 a.u. cut-off radius, 100 Ry energy cut-off, and DZP orbitals containing 2s2p1d).

### 1. Calculating the second-order force constants

It would be best to combine Phonopy and ASE with ABACUS to calculate the second-order force constants. First, enter the `2nd` folder.

### 1.1 Structure optimization

Before performing lattice thermal conductivity calculations, it is necessary to optimize the atomic configuration of the simulated material system. The following is the atomic configuration file `STRU` obtained after structure optimization (relax) using ABACUS. In this example, for simplicity, a 2\*2\*2 Brillouin zone k-point sampling was used in the structure optimization process, with an energy cutoff value of 100 Ry for plane waves (the plane wave basis vector is also used in LCAO). Note that the actual calculation should use more convergent k-point sampling.

```bash
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.88972612546

LATTICE_VECTORS
0 2.81594778072 2.81594778072 #latvec1
2.81594778072 0 2.81594778072 #latvec2
2.81594778072 2.81594778072 0 #latvec3

ATOMIC_POSITIONS
Direct # direct coordinate

Si #label
0 #magnetism
2 #number of atoms
0.875  0.875  0.875  m  0  0  0
0.125  0.125  0.125  m  0  0  0
```

### 1.2 Calculating the second-order force constants

The Phonopy software is called to generate multiple atomic configurations of the supercell and corresponding perturbations needed for calculation with the following command:

```bash
phonopy setting.conf --abacus -d
```

where the `setting.conf` file reads:

```bash
DIM = 2 2 2
ATOM_NAME = Si
```

In this Si example, we only need to generate one perturbed atomic configuration, `STRU-001`. Perform SCF calculations (SCF stands for Self-Consistent Field and represents the iterative self-consistent calculation of density functional theory) on all perturbed configurations (in this case, there is only one for Si) to obtain the forces on the atoms. Afterward, use the following command to generate the `FORCE_SET` file:

```bash
phonopy -f OUT.DIA-50/running_scf.log
```

Tip: In the input file `INPUT` of ABACUS, you can set the variable `stru_file`, which corresponds to the atomic configuration file `STRU-001`, and ABACUS will read the structure file directly.

Next, set the `band.conf` file to calculate the phonon spectrum and the second-order force constants:

```bash
phonopy -p band.conf --abacus
```

The `band.conf` file mentioned here contains the following contents (you can refer to the Phonopy documentation for specific parameter meanings):

```bash
ATOM_NAME = Si
DIM = 2 2 2
MESH = 8 8 8
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.5  0.625  0.25  0.625, 0.375 0.375 0.75  00 0.0 0.0  0.5 0.5 0.5
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
FORCE_CONSTANTS = WRITE
FULL_FORCE_CONSTANTS = .TRUE.
```

After this step, the Phonopy software will generate `band.yaml` (for plotting the phonon spectrum) and the `FORCE_CONSTANTS` file. The data contained in the `FORCE_CONSTANTS` file is the second-order force constants. It is important to set `FULL_FORCE_CONSTANTS = .TRUE.`, which outputs all the second-order force constants. Otherwise, there may be errors when ShengBTE reads the data.

In addition, you can use the following command to output the gnuplot format of the phonon spectrum for plotting:

```bash
phonopy-bandplot --gnuplot > pho.dat
```

### 1.3 Post-processing

Note that ShengBTE software requires the unit of the data in the `FORCE_CONSTANTS_2ND` file to be eV/Å^2, but the unit of the `FORCE_CONSTANTS` calculated by ABACUS combined with Phonopy is eV/(Å*au), where au is the atomic unit system and 1 au = 0.52918 Å. You can use the provided `au2si.py` script in the `2nd` directory to convert the units and generate the `FORCE_CONSTANTS_2ND` file. The command is as follows:

```python
python au2si.py
```

The `FORCE_CONSTANTS_2ND` file is provided in the shengbte folder for reference to the calculation results.

### 2. Calculating the third-order force constants

To calculate the third-order force constants, you need to combine with the thirdorder program and output the third-order force constant file `FORCE_CONSTANTS_3RD`. However, thirdorder currently only supports reading input and output files from VASP and QE. Therefore, we are using thirdorder by converting ABACUS's structure and output force files to `POSCAR` and `vasprun.xml`, respectively. Please enter the `3rd` folder first, and the specific steps will be described below.

### 2.1 Obtaining perturbed configurations

First, convert the optimized `STRU` file from ABACUS software to `POSCAR` (the converted `POSCAR` file is already provided in the directory, or you can do this conversion by yourself).

Then, run the `thirdorder_vasp.py` program to generate a series of atomic configuration files `3RD.POSCAR.*` after perturbation. For example, in this example, a total of 40 configurations were generated:

```bash
thirdorder_vasp.py sow 2 2 2 -2
```

Run `pos2stru.py` to convert the above `POSCAR` to `STRU` file. Note that this script calls functions from the ASE software package (ASE needs to be installed in advance):

```python
python pos2stru.py
```

Note: The dpdata software cannot be called here to perform the conversion. This is because the dpdata software forces the lattice to change into a lower triangular matrix, which is equivalent to rotating the lattice and leads to a corresponding rotation in the direction of the interatomic forces, which will cause errors.

### 2.2 Calculation of atomic forces for perturbation configurations

You can refer to the `run_stru.sh` script provided in the directory to batch generate `SCF-*` folders and submit calculations. Here, ABACUS needs to perform SCF calculations on 40 atomic configurations, which may take some time. It is recommended to run each SCF separately in the `SCF-*` folder. The `scf_thr` parameter in INPUT file should be set to at least 1e-8 to obtain converged results.

After the calculations are complete, run `aba2vasp.py` to package the atomic forces calculated by ABACUS into the `vasprun.xml` format and place them in each `SCF-\*` folder with the following command:

```python
python aba2vasp.py
```

The `vasprun.xml` format is illustrated as follows:

```xml
<modeling>
    <calculation>
        <varray name="forces">
            <v>1.865e-05 -0.04644196 -0.00153852</v>
            <v>-1.77e-05 -0.00037715 -0.00149635</v>
            <v>1.973e-05 0.002213 -0.00149461</v>
            <v>-1.976e-05 0.00065303 -0.0014804</v>
            <v>8.31e-06 -0.0003306 -0.00024288</v>
            <v>-8.25e-06 -0.00038306 -0.00025385</v>
            <v>1.071e-05 0.00060621 -0.00025797</v>
            <v>-1.05e-05 -0.00014553 -0.00027532</v>
            <v>0.00668053 0.00645634 -0.04642593</v>
            <v>-0.00668085 0.00645595 -0.00040122</v>
            <v>-0.00650454 0.00628877 -0.00025123</v>
            <v>0.00650504 0.00628892 -0.00028948</v>
            <v>-0.00039591 2.479e-05 0.00223371</v>
            <v>0.00039608 2.426e-05 0.0006732</v>
            <v>0.0003264 3.122e-05 0.00052874</v>
            <v>-0.00032589 3.415e-05 -0.00023577</v>
            <v>-2.908e-05 -0.00832477 0.00635709</v>
            <v>3.737e-05 -0.00125057 -7.444e-05</v>
            <v>-2.582e-05 0.00656076 0.00636285</v>
            <v>2.566e-05 -0.00049974 -6.661e-05</v>
            <v>-5.431e-05 0.00502637 0.00639077</v>
            <v>4.553e-05 -0.00180978 0.0001325</v>
            <v>-3.609e-05 -0.00676473 0.00638092</v>
            <v>3.806e-05 5.503e-05 0.00012759</v>
            <v>-0.00670704 0.00646596 0.01310437</v>
            <v>0.00670119 3.673e-05 0.00602948</v>
            <v>0.00036366 0.00627899 -0.00657272</v>
            <v>-0.00036508 2.288e-05 0.00026009</v>
            <v>0.00648649 0.0064463 -0.00036521</v>
            <v>-0.00648098 1.594e-05 0.00671469</v>
            <v>-0.00034493 0.00630074 0.00662932</v>
            <v>0.00034331 4.157e-05 -0.0002028</v>
        </varray>
    </calculation>
</modeling>
```

Finally, execute the following command:

```bash
find SCF-* -name vasprun.xml|sort -n|thirdorder_vasp.py reap 2 2 2 -2
```

Then, the third-order force constant file `FORCE_CONSTANTS_3RD` can be obtained by running the above command. The `FORCE_CONSTANTS_3RD` file is provided in the shengbte folder for reference in calculating the results.

### 3. Run ShengBTE to obtain lattice thermal conductivity

Enter the `shengbte` folder, in which the three files `CONTROL` (parameter file of ShengBTE), `FORCE_CONSTANTS_2ND` (second-order force constant file), and `FORCE_CONSTANTS_3RD` (third-order force constant file) have been prepared. Next, run ShengBTE with the following command to obtain the lattice thermal conductivity, where the calculation results are given in the Ref folder for reference:

```bash
mpirun -n 10 ShengBTE
```

## Conclusion

For using plane wave (PW) in ABACUS to perform ShengBTE calculations, similar procedures should be followed. However, the `scf_thr` parameter in the `INPUT` file for calculating the third-order force constant needs to be set to at least 1e-12. The experimental lattice thermal conductivity of Si at 300 K is around 150 W/(m K), while the calculated thermal conductivity of Si at 300 K is around 100 W/(m K) by using the provided example. This is because, as a demo, a 2\*2\*2 expanded cell and a 2\*2\*2 K-point are used in the example, but the results are not converged yet with respect to the given system size and k-points. In actual research, the size of the supercell and the sampling scheme of K-points need to be tested to obtain converged results.

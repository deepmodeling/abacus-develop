# Calculating DOS and PDOS

## DOS

ABACUS can calculate the density of states (DOS) of the system, and the examples can be found in [examples/dos](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dos).
We first, do a ground-state energy calculation ***with one additional keyword "[out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg)" in the INPUT file***:

```
out_chg              1
```

this will produce the converged charge density, which is contained in the file SPIN1_CHG.cube.
Then, use the same `STRU` file, pseudopotential file and atomic orbital file (and the local density matrix file onsite.dm if DFT+U is used) to do a non-self-consistent calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:

```
INPUT_PARAMETERS
#Parameters (General)
suffix Si2_diamond
ntype 1
nbands 8
calculation nscf
basis_type lcao
read_file_dir   ./

#Parameters (Accuracy)
ecutwfc 60
symmetry 1
scf_nmax 50
scf_thr 1.0e-9
pw_diag_thr 1.0e-7

#Parameters (File)
init_chg file
out_dos 1
dos_sigma 0.07
```

Some parameters in the INPUT file are explained:

- calculation

  choose which kind of calculation: scf calculation, nscf calculation, structure relaxation or Molecular Dynamics. Now we need to do one step of nscf calculation.
  Attention: This is a main variable of ABACUS, and for its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#calculation).
- pw_diag_thr

  threshold for the CG method which diagonalizes the Hamiltonian to get eigenvalues and eigen wave functions. If one wants to do nscf calculation, pw_diag_thr needs to be changed to a smaller account, typically smaller than 1.0e-3. Note that this parameter only apply to plane-wave calculations that employ the CG or Davidson method to diagonalize the Hamiltonian. For its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#pw_diag_thr).

  For LCAO calculations, this parameter will be neglected !
- init_chg

  the type of starting density. When doing scf calculation, this variable can be set ”atomic”. When doing nscf calculation, the charge density already exists(eg. in SPIN1_CHG.cube), and the variable should be set as ”file”. It means the density will be read from the existing file SPIN1_CHG.cube. For its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#init_chg).
- out_dos

  output density of state(DOS). The unit of DOS is `(number of states)/(eV * unitcell)`. For its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#out_dos).
- dos_sigma

  the gaussian smearing parameter(DOS), in unit of eV. For its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#dos_sigma).
- read_file_dir

  the location of electron density file. For its more information please see the [here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#read_file_dir).

To have an accurate DOS, one needs to have a denser k-point mesh. For example, the KPT file can be set as:

```
K_POINTS
0
Gamma
8 8 8 0 0 0
```

Run the program, and you will see a file named DOS1_smearing.dat in the output directory. The first two columns in the file are the energy and DOS, respectively, and the third column is the sum of DOS. Plot file DOS1_smearing.dat with graphing software, and you’ll get the DOS.

```
            -5.49311           0.0518133           0.0518133
            -5.48311           0.0641955            0.116009
            -5.47311           0.0779299            0.193939
            -5.46311           0.0926918             0.28663
            -5.45311            0.108023            0.394653
            -5.44311            0.123346            0.517999
            ...
```

## PDOS

Along with the DOS1_smearing.dat file, we also produce the projected density of states (PDOS) in a file called PDOS.

The PDOS file starts with number of atomic orbitals in the system, then a list of energy values, such as:

```
<pdos>
<nspin>1</nspin>
<norbitals>26</norbitals>
<energy_values units="eV">
            -5.50311
            -5.49311
            -5.48311
            -5.47311
...

```

The rest of the fileis arranged in sections, each section with a header such as below:

```
<orbital
 index="                                       1"
 atom_index="                                       1"
 species="Si"
 l="                                       0"
 m="                                       0"
 z="                                       1"
>
<data>
...
</data>

```

which tells the atom and symmetry of the current atomic orbital, and followed by the PDOS values. The values can thus be plotted against the energies. The unit of PDOS is also `(number of states)/(eV * unitcell)`.

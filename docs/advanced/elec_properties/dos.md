# Calculating DOS and PDOS

## DOS

ABACUS can calculate the density of states (DOS) of the system, and the examples can be found in [examples/dos](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dos).
We first, do a ground-state energy calculation ***with one additional keyword "[out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg)" in the INPUT file***:

```
out_chg              1
```

this will produce the converged charge density, which is contained in the file SPIN1_CHG.cube.
Then, use the same `STRU` file, pseudopotential file and atomic orbital file (and the local density matrix file dm_onsite.txt if DFT+U is used) to do a non-self-consistent calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:

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

Run the program, and you will see a file named doss1g1_nao.txt in the output directory. The columns are: energy(eV), dos(1/eV), dos_int (integrated DOS), dos_smear(1/eV), dos_smear_int (integrated smeared DOS). Plot the file with graphing software, and you'll get the DOS.

```
#   energy(eV)      dos(1/eV)        dos_int  dos_smear(1/eV)  dos_smear_int
      -5.49311       0.0518133       0.0518133       0.0518133       0.0518133
      -5.48311       0.0641955        0.116009       0.0641955        0.116009
      -5.47311       0.0779299        0.193939       0.0779299        0.193939
      ...
```

## PDOS

Along with the DOS files, we also produce the projected density of states (PDOS) in files named pdoss{spin}g{geom}_{basis}.txt (e.g., pdoss1g1_nao.txt).

The PDOS file uses a plain-text format. Each row corresponds to one energy point and one atom. Columns: energy(eV), atom (1-based), species, then pdos values ordered as s(1), p(3), d(5), f(7), etc. Zeta components are summed, and values below 1e-6 are zeroed out.

```
# istep: 1
# npoints: 2736
# energy(eV)  atom  species  pdos(1/eV), columns: s(m=0) p(m=0,+1,-1) d(m=0,+1,-1,+2,-2) f(m=0,+1,-1,+2,-2,+3,-3)
  -55.607730    1     Fe    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  ...
```

The `# istep:` line is only present for geometry-optimization (istep) output. Within each angular momentum `l`, the `2l+1` columns follow the magnetic-quantum-number order `m = 0, +1, -1, +2, -2, ...`.

For nspin=2, two files are written (pdoss1* and pdoss2*), one per spin channel. For nspin=4, the two spinor components are summed into a single file. The unit of PDOS is also `(number of states)/(eV * unitcell)`.

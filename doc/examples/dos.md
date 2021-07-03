# DOS, wave functions

[back to main page](../../README.md)

The main task of this example is to calculate the density of states (DOS) of the system. At first, do a ground-state energy calculation as in [this example](#basic-lcao.md) ***with one additional keyword in the INPUT file***:

```
out_charge              1
```

this will produce the converged charge density, which is contained in the file SPIN1_CHG. Copy the file along with the `STRU` file, the pseudopotential file and the atomic orbital file to the new working directory where we will do a non-self-consistent calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:
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
niter 50
dr2 1.0e-9
ethr 1.0e-7

#Parameters (File)
start_charge file
out_dos 1
dos_sigma 0.07

#Parameters (Smearing)
smearing gaussian
sigma 0.02

```

Some parameters in the INPUT file are explained:
- calculation

    choose which kind of calculation: scf calculation, nscf calculation, structure relaxation or Molecular Dynamics. Now we need to do one step of nscf calculation.
    Attention: This is a main variable of ABACUS, and for its more information please see the [list of input variables](../input-main.md).
- ethr

    threshold for the CG method which diagonalizes the Hamiltonian to get eigenvalues and eigen wave functions. If one wants to do nscf calculation, ethr needs to be changed to a smaller account, typically smaller than 1.0e-3. Note that this parameter only apply to plane-wave calculations that employ the CG method to diagonalize the Hamiltonian.
    
    For LCAO calculations, this parameter will be neglected !
- start_charge

    the type of starting density. When doing scf calculation, this variable can be set ”atomic”. When doing nscf calculation, the charge density already exists(eg. in SPIN1_CHG), and the variable should be set as ”file”. It means the density will be read from the existing file SPIN1_CHG. For more information please see the [list of input variables](../input-main.md).

- out_dos

    output density of state(DOS).

- dos_sigma

    the gaussian smearing parameter(DOS).

- read_file_dir

    the location of electron density file.

To have an accurate DOS, one needs to have a denser k-point mesh. For example, the KPT file can be set as:
```
K_POINTS
0
Gamma
8 8 8 0 0 0
```
Run the program, and you will see a file named DOS1_smearing.dat in the output directory.

Next, plot file DOS1_smearing.dat with graphing software, and you’ll get the DOS.

[back to top](#dos,-wave-functions)
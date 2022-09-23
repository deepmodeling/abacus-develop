# Two Quick Examples

## Running SCF Calculation

Here, we take FCC MgO as a quick start example. Firstly we need to build the structure file. The default name of a structure file in ABACUS is `STRU`, e.g.:

```
#This is the atom file containing all the information
#about the lattice structure.

ATOMIC_SPECIES
Mg 24.305  Mg_ONCV_PBE-1.0.upf  # element name, atomic mass, pseudopotential file
O  15.999 O_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
8.04411131              # 3.01/sqrt(2.)/0.52918

LATTICE_VECTORS
1.00000000 0.00000000 0.00000000                # 1/0.52918 = 1.88971616
0.00000000 1.00000000 0.00000000
0.00000000 0.00000000 1.00000000

ATOMIC_POSITIONS
Direct                  #Cartesian(Unit is LATTICE_CONSTANT)
Mg                      #Name of element        
0.0                     #Magnetic for this element.
4                       #Number of atoms
0.0  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
O                       #Name of element        
0.0                     #Magnetic for this element.
4                       #Number of atoms
0.5  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
```

Secondly, one needs to set the `INPUT` file, which sets all key parameters to direct ABACUS how to calculte and what to output:
```
INPUT_PARAMETERS
#Parameters     (System)
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
calculation             scf		# this is the key parameter telling abacus to do a scf calculation
#Parameters (PW)
ecutwfc                 100             # Rydberg
scf_thr                 0.01		# Rydberg

#Parameters (electronic opt)
basis_type              pw              # or lcao, lcao_in_pw
ks_solver               cg              # & david for pw; genelpa, hpseps & lapack (single cpu) for lcao
smearing_method         fixed           # or gauss, mp
smearing_sigma          0.01
mixing_type             pulay           # or kerker, plain, pulay-kerker
mixing_beta             0.7
```
The pseudopotential files of `Mg_ONCV_PBE-1.0.upf` and `O_ONCV_PBE-1.0.upf` should be provided under the directory of `pseudo_dir`.

The next mandatory input file is called `KPT`, which sets the k-mesh. Below is an example:

```
K_POINTS
0 
Gamma
4 4 4 0 0 0
```

Ok, now that all key input files have been set, we should be able to run the first quick example. The simplest way is to use the command line, e.g.:

```
mpirun -np 4 abacus
```

The main output information appears in `OUT.MgO/running_scf.log`, which starts with

```

                             WELCOME TO ABACUS

               'Atomic-orbital Based Ab-initio Computation at UStc'

                     Website: http://abacus.ustc.edu.cn/

    Version: Parallel, in development
    Processor Number is 4
    Start Time is Tue Sep 22 18:23:54 2022

 ------------------------------------------------------------------------------------

 READING GENERAL INFORMATION
                           global_out_dir = OUT.MgO/
                           global_in_card = INPUT
                               pseudo_dir =
                              orbital_dir =
                              pseudo_type = auto
                                    DRANK = 1
                                    DSIZE = 8
                                   DCOLOR = 1
                                    GRANK = 1
                                    GSIZE = 1




 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 |                                                                    |
 | Reading atom information in unitcell:                              |
 | From the input file and the structure file we know the number of   |
 | different elments in this unitcell, then we list the detail        |
 | information for each element, especially the zeta and polar atomic |
 | orbital number for each element. The total atom number is counted. |
 | We calculate the nearest atom distance for each atom and show the  |
 | Cartesian and Direct coordinates for each atom. We list the file   |
 | address for atomic orbitals. The volume and the lattice vectors    |
 | in real and reciprocal space is also shown.                        |
 |                                                                    |
 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

```

If ABAUCS finishes successfully, you will be able to obtain the total energy in `OUT.MgO/running_scf.log`:

```
 --------------------------------------------
 !FINAL_ETOT_IS -1908.456553781335 eV
 --------------------------------------------
```

Congratulations! ABACUS has come under you command from now on!


## Running Geometry Optimization

Let's consider how to run geometry optimization in ABACUS. As you may have guessed, one have only to set `calculation` to `cell-relax`. But wait a bit, how should we set the convergence criteria for atomics force and cell stress? What's the maximum number of ionc steps in one optimization trail? Well, one contact `INPUT` reads like:

```
INPUT_PARAMETERS
#Parameters     (System)
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
calculation             cell-relax	# this is the key parameter telling abacus to do a optimization calculation
force_thr_ev		0.01		# the threshold of the force convergence, in unit of eV/Angstrom
stress_thr		1		# the threshold of the stress convergence, in unit of kBar
relax_nmax		100		# the maximal number of ionic iteration steps
#Parameters (PW)
ecutwfc                 100             # Rydberg
scf_thr                 0.01		# Rydberg

#Parameters (electronic opt)
basis_type              pw              # or lcao, lcao_in_pw
ks_solver               cg              # & david for pw; genelpa, hpseps & lapack (single cpu) for lcao
smearing_method         fixed           # or gauss, mp
smearing_sigma          0.01
mixing_type             pulay           # or kerker, plain, pulay-kerker
mixing_beta             0.7
```

Now we can use the same `KPT` and `STRU` file as above. Run this optimization example, the final optimized structure will appear under the directory `OUT.MgO` with the name of `STRU_NOW.cif`.

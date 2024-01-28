# Two Quick Examples

## Running SCF Calculation

### A quick LCAO example

ABACUS is well known for its support of LCAO (Linear Combination of Atomic Orbital) basis set in calculating periodic condensed matter systems, so it's a good choice to start from a LCAO example of self-consistent field (SCF) calculation. Here, FCC MgO has been chosen as a quick start example. The default name of a structure file in ABACUS is `STRU`. The `STRU` file for FCC MgO in a LCAO calculation is shown below:

```
#This is the atom file containing all the information
#about the lattice structure.

ATOMIC_SPECIES
Mg 24.305  Mg_ONCV_PBE-1.0.upf  # element name, atomic mass, pseudopotential file
O  15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897259886 		# 1.8897259886 Bohr =  1.0 Angstrom

LATTICE_VECTORS
4.25648 0.00000 0.00000  
0.00000 4.25648 0.00000
0.00000 0.00000 4.25648

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

Next, the `INPUT` file is required, which sets all key parameters to direct ABACUS how to calculte and what to output:
```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
pseudo_dir              ./
orbital_dir		./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              lcao            
calculation             scf		# this is the key parameter telling abacus to do a scf calculation
```

The pseudopotential files of `Mg_ONCV_PBE-1.0.upf` and `O_ONCV_PBE-1.0.upf` should be provided under the directory of `pseudo_dir`, and the orbital files `Mg_gga_8au_100Ry_4s2p1d.orb` and `O_gga_8au_100Ry_2s2p1d.orb` under the directory of `orbital_dir`. The pseudopotential and orbital files can be downloaded from the [ABACUS website](http://abacus.ustc.edu.cn/pseudo/list.htm).

The final mandatory input file is called `KPT`, which sets the reciprocal space k-mesh. Below is an example:

```
K_POINTS
0 
Gamma
4 4 4 0 0 0
```

After all the above input files have been set, one should be able to run the first quick example. The simplest way is to use the command line, e.g.:

```
mpirun -np 2 abacus
```

The main output information is stored in the file `OUT.MgO/running_scf.log`, which starts with

```
                             WELCOME TO ABACUS v3.2

               'Atomic-orbital Based Ab-initio Computation at UStc'

                     Website: http://abacus.ustc.edu.cn/

    Version: Parallel, in development
    Processor Number is 2
    Start Time is Mon Oct 24 01:47:54 2022

 ------------------------------------------------------------------------------------

 READING GENERAL INFORMATION
                           global_out_dir = OUT.MgO/
                           global_in_card = INPUT
                               pseudo_dir =
                              orbital_dir =
                                    DRANK = 1
                                    DSIZE = 2
                                   DCOLOR = 1
                                    GRANK = 1
                                    GSIZE = 1
 The esolver type has been set to : ksdft_lcao




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

......

```

If ABAUCS finishes successfully, the total energy will be output in `OUT.MgO/running_scf.log`:

```
 --------------------------------------------
 !FINAL_ETOT_IS -7663.897267807250 eV
 --------------------------------------------
```

### A quick PW example

In order to run a SCF calculation with PW (Plane Wave) basis set, one has only to change the tag `basis_type` from `lcao` to `pw` in the `INPUT` file, and no longer needs to provide orbital files under `NUMERICAL_ORBITAL` in the `STRU` file.

The `INPUT` file follows as:
```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
pseudo_dir              ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              pw              # changes the type of basis set
calculation             scf		# this is the key parameter telling abacus to do a scf calculation
```

And the `STRU` file will be:

```
#This is the atom file containing all the information
#about the lattice structure.

ATOMIC_SPECIES
Mg 24.305  Mg_ONCV_PBE-1.0.upf  # element name, atomic mass, pseudopotential file
O  15.999 O_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897259886 		# 1.8897259886 Bohr =  1.0 Angstrom

LATTICE_VECTORS
4.25648 0.00000 0.00000  
0.00000 4.25648 0.00000
0.00000 0.00000 4.25648

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

Use the same pseudopotential and `KPT` files as the above LCAO example. The final total energy will be output:

```
 --------------------------------------------
 !FINAL_ETOT_IS -7665.688319476949 eV
 --------------------------------------------
```

## Running Geometry Optimization

In order to run a full geometry optimization in ABACUS, the tag `calculation` in `INPUT` should be set to `cell-relax`. In addition, the convergence criteria for atomics force and cell stress can be set through the tags `force_thr_ev` and `stress_thr`, respectively. The maximum number of ionc steps is controlled by `relax_nmax`.

### A quick LCAO example

The `INPUT` is provided as follows:

```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
orbital_dir             ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              lcao 
calculation             cell-relax	# this is the key parameter telling abacus to do a optimization calculation
force_thr_ev		0.01		# the threshold of the force convergence, in unit of eV/Angstrom
stress_thr		5		# the threshold of the stress convergence, in unit of kBar
relax_nmax		100		# the maximal number of ionic iteration steps
out_stru		1
```
Use the same `KPT`, `STRU`, pseudopotential, and orbital files as in the above SCF-LCAO example. The final optimized structure can be found in `STRU_NOW.cif` and `OUT.MgO/running_cell-relax.log`.

### A quick PW example

The `INPUT` is provided as follows:

```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              pw
calculation             cell-relax	# this is the key parameter telling abacus to do a optimization calculation
force_thr_ev		0.01		# the threshold of the force convergence, in unit of eV/Angstrom
stress_thr		5		# the threshold of the stress convergence, in unit of kBar
relax_nmax		100		# the maximal number of ionic iteration steps
out_stru		1
```

Use the same `KPT`, `STRU`, and pseudopotential files as in the above SCF-PW examples. The final optimized structure can be found in `STRU_NOW.cif` and `OUT.MgO/running_cell-relax.log`.

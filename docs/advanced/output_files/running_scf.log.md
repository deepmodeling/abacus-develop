# The running_scf.log file

- [The running\_scf.log file](#the-running_scf.log-file)
  - [Reading information](#reading-information)
    - [Reading version information](#reading-version-information)
    - [Reading general information](#reading-general-information)
    - [Reading unitcell](#reading-unitcell)
    - [Reading pseudopotentials files](#reading-pseudopotentials-files)
  - [Setup Tasks](#setup-tasks)
    - [Setup plane waves of charge/potential](#setup-plane-waves-of-chargepotential)
    - [Doing symmetry analysis](#doing-symmetry-analysis)
    - [Setup K-points](#setup-k-points)
    - [Setup plane waves of wave functions](#setup-plane-waves-of-wave-functions)
  - [Running scf processes](#running-scf-processes)
    - [Search adjacent atoms (init)](#search-adjacent-atoms-init)
    - [Scf iteration](#scf-iteration)
    - [Scf results](#scf-results)
  - [Results summary](#results-summary)
    - [TOTAL-FORCE](#total-force)
    - [TOTAL-STRESS](#total-stress)
    - [FINAL\_ETOT\_IS](#final_etot_is)
    - [TIME STATISTICS](#time-statistics)
    - [MEMORY STATISTICS](#memory-statistics)
    - [Start and end times](#start-and-end-times)

  

  
## Reading information
### Reading version information
```
ABACUS v3.7.0
Atomic-orbital Based Ab-initio Computation at UStc 
// ABACUS version and description.

Website: http://abacus.ustc.edu.cn/                              
// Official website for more information about ABACUS.
               
Documentation: https://abacus.deepmodeling.com/                        
// Documentation provides detailed user guides and references.

Repository: https://github.com/abacusmodeling/abacus-develop        
// GitHub repository where the source code is hosted.
             https://github.com/deepmodeling/abacus-develop          
 // Another link to the same repository, ensuring access.

Commit: a339356 (Thu Jun 27 16:40:42 2024 +0800)
// The specific commit of the ABACUS source code used for this run.

Start Time is Thu Jul 18 11:34:56 2024
// The start time of the ABACUS calculation.
```



------------------------------------------------------------------------------------
### Reading general information
```
READING GENERAL INFORMATION
// The following section is reading general settings and preparing the computation.

global_out_dir = OUT.CeAl/
// The directory where output files will be stored.

global_in_card = INPUT
// The main input file for the ABACUS calculation.

pseudo_dir = 
orbital_dir = 
// Directories for pseudopotential and orbital files, if used.

DRANK = 1
DSIZE = 8
DCOLOR = 1
GRANK = 1
GSIZE = 1
// Parameters related to the domain decomposition and parallel execution.

The esolver type has been set to : ksdft_lcao
// The type of electronic structure solver being used in this calculation.

RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
// The device and CPU model used for the calculation.
```

### Reading unitcell
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>=
|                                                                    |
| Reading atom information in unitcell:                              |
| From the input file and the structure file we know the number of   |
| different elements in this unitcell, then we list the detail       |
| information for each element, especially the zeta and polar atomic |
| orbital number for each element. The total atom number is counted. |
| We calculate the nearest atom distance for each atom and show the  |
| Cartesian and Direct coordinates for each atom. We list the file   |
| address for atomic orbitals. The volume and the lattice vectors    |
| in real and reciprocal space is also shown.                        |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

READING UNITCELL INFORMATION
// The following section reads the unit cell information necessary for the calculation.

ntype = 2
// There are two different types of atoms in the unit cell.

lattice constant (Bohr) = 1.88973
// The lattice constant of the unit cell in Bohr units.

lattice constant (Angstrom) = 1
// The lattice constant of the unit cell in Angstrom units.

READING ATOM TYPE 1
// The section starts reading details about the first type of atom in the unit cell.

atom label = Ce
// The label for the first type of atom, which is Cerium (Ce).

L=0, number of zeta = 6
// For angular momentum quantum number L=0, there are 6 zeta functions defined for the atomic orbitals.

L=1, number of zeta = 3
// For L=1, there are 3 zeta functions.

L=2, number of zeta = 3
// For L=2, there are also 3 zeta functions.

L=3, number of zeta = 3
// Similarly, for L=3, there are 3 zeta functions.

L=4, number of zeta = 2
// And for L=4, there are 2 zeta functions.

number of atom for this type = 8
// There are 8 atoms of this type in the unit cell.

READING ATOM TYPE 2
// The section now reads details about the second type of atom in the unit cell.

atom label = Al
// The label for the second type of atom, which is Aluminum (Al).

L=0, number of zeta = 4
// For L=0, there are 4 zeta functions for Al.

L=1, number of zeta = 4
// For L=1, there are 4 zeta functions.

L=2, number of zeta = 1
// And for L=2, there is 1 zeta function.

number of atom for this type = 16
// There are 16 atoms of this type in the unit cell.

TOTAL ATOM NUMBER = 24
// The total number of atoms in the unit cell is 24.

DIRECT COORDINATES
    atom                   x                   y                   z     mag                  vx                  vy                  vz
// The following table lists the direct coordinates (in units of the lattice constant), magnetic quantum number, and velocities (vx, vy, vz) for each atom.

[...]
// The actual coordinates and velocities for each atom are listed here.

                          Volume (Bohr^3) = 3395.25
// The volume of the unit cell in Bohr^3.

                             Volume (A^3) = 503.123
// The volume of the unit cell in cubic Angstroms.

Lattice vectors: (Cartesian coordinate: in unit of a_0)
// The lattice vectors of the unit cell in Cartesian coordinates, in units of the Bohr radius (a_0).

Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)
// The reciprocal lattice vectors, also in Cartesian coordinates, in units of 2π/a_0.
```

### Reading pseudopotentials files
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Reading pseudopotentials files:                                    |
| The pseudopotential file is in UPF format. The 'NC' indicates that |
| the type of pseudopotential is 'norm conserving'. Functional of    |
| exchange and correlation is decided by 4 given parameters in UPF   |
| file.  We also read in the 'core correction' if there exists.      |
| Also we can read the valence electrons number and the maximal      |
| angular momentum used in this pseudopotential. We also read in the |
| trail wave function, trail atomic density and local-pseudopotential|
| on logrithmic grid. The non-local pseudopotential projector is also|
| read in if there is any.                                           |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

PAO radial cut off (Bohr) = 15
// The radial cutoff for pseudo atomic orbitals is given in Bohr units.

Read in pseudopotential file is Ce-sp.PD04.PBE.UPF
// The filename of the pseudopotential file used for Cerium (Ce) in UPF format.

                     pseudopotential type = NC
// Indicates that the pseudopotential is of 'norm conserving' (NC) type.

          exchange-correlation functional = PBE
// Specifies that the Perdew-Burke-Ernzerhof (PBE) functional is used for exchange-correlation.

                 nonlocal core correction = 1
// Indicates that there is a nonlocal core correction included in the pseudopotential.

                        valence electrons = 12
// The number of valence electrons considered in the pseudopotential, which affects the accuracy and efficiency of the calculation.

                                     lmax = 3
// The maximum angular momentum quantum number l used in the pseudopotential.

                           number of zeta = 5
// The number of zeta functions used in the pseudopotential.

                     number of projectors = 8
// The number of projectors in the pseudopotential, which is related to the nonlocal part.

                           L of projector = 0, 0, 1, 1, 2, 2, 3, 3
// The angular momentum quantum numbers for each projector in the pseudopotential.

                PAO radial cut off (Bohr) = 15
// The radial cutoff for pseudo atomic orbitals (PAOs) in Bohr units, which defines the range of the orbitals.

Read in pseudopotential file is Al_ONCV_PBE-1.0.upf
// The filename of the pseudopotential file used for Aluminum (Al) in UPF format.

                     pseudopotential type = NC
// Again, indicates that the pseudopotential is of 'norm conserving' (NC) type.

          exchange-correlation functional = PBE
// The PBE functional is used for exchange-correlation for Aluminum as well.

                 nonlocal core correction = 0
// Indicates that there is no nonlocal core correction for Aluminum.

                        valence electrons = 11
// The number of valence electrons considered in the pseudopotential for Aluminum.

                                     lmax = 1
// The maximum angular momentum quantum number l for Aluminum is 1.

                           number of zeta = 0
// There are no additional zeta functions for Aluminum in this pseudopotential.

                     number of projectors = 4
// The number of projectors for the Aluminum pseudopotential.

                           L of projector = 0, 0, 1, 1
// The angular momentum quantum numbers for each projector in the Aluminum pseudopotential.

     initial pseudo atomic orbital number = 136
// The initial number of pseudo atomic orbitals used in the calculation.

                                   NLOCAL = 888
// The number of nonlocal operations, which may be related to the nonlocal pseudopotential projectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Warning: the number of valence electrons in pseudopotential > 4 for Ce: [Xe] 4f1 5d1 6s2
 Warning: the number of valence electrons in pseudopotential > 3 for Al: [Ne] 3s2 3p1
 Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
 If you‘re confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// A warning message indicating that the pseudopotentials include more valence electrons than typical, which can improve accuracy but at the cost of computational efficiency.

Warning_Memory_Consuming allocated:  FFT::grid 6.5918 MB
// A warning about the memory consumption for the Fast Fourier Transform (FFT) grid.

```
## Setup Tasks
### Setup plane waves of charge/potential
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Setup plane waves of charge/potential:                             |
| Use the energy cutoff and the lattice vectors to generate the      |
| dimensions of FFT grid. The number of FFT grid on each processor   |
| is 'nrxx'. The number of plane wave basis in reciprocal space is   |
| different for charge/potential and wave functions. We also set    |
| the 'sticks' for the parallel of FFT. The number of plane waves    |
| is 'npw' in each processor.                                        |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SETUP THE PLANE WAVE BASIS
energy cutoff for charge/potential (unit:Ry) = 600
// The energy cutoff used for charge and potential calculations in Rydberg units.

            fft grid for charge/potential = [ 120, 120, 120 ]
// The dimensions of the FFT grid used for charge and potential calculations.

                        fft grid division = [ 3, 3, 3 ]
// The division of the FFT grid among processors.

        big fft grid for charge/potential = [ 40, 40, 40 ]
// The dimensions of the 'big' FFT grid, which might be used for parallelization.

                                     nbxx = 8000
// The total number of plane waves in the charge/potential basis.

                                     nrxx = 216000
// The number of FFT grid points on each processor.

SETUP PLANE WAVES FOR CHARGE/POTENTIAL
                    number of plane waves = 842641
// The total number of plane waves used for charge/potential.

                         number of sticks = 10781
// The number of 'sticks' used in the FFT parallelization.

PARALLEL PW FOR CHARGE/POTENTIAL
// The distribution of plane waves among processors for charge/potential calculations.

[...]
// The actual distribution of plane waves and the number of plane waves per processor.

--------------- sum -------------------
        8          10781         842641
// The sum of plane waves and the total number of processors used.

                            number of |g| = 2844
// The total number of G-vectors in the reciprocal space.

                                  max |g| = 54.2697
// The maximum magnitude of the G-vectors.

                                  min |g| = 0.0790412
// The minimum magnitude of the G-vectors.

----------- Double Check Mixing Parameters Begin ------------
mixing_type: pulay
// The type of charge mixing used in the SCF procedure, which is Pulay mixing.

mixing_beta: 0.4
// The Pulay mixing beta parameter, which controls the degree of mixing.

mixing_gg0: 0
mixing_gg0_min: 0.1
mixing_ndim: 8
// Parameters related to the Pulay mixing method.

----------- Double Check Mixing Parameters End ------------

SETUP THE ELECTRONS NUMBER
// The section where the total number of electrons is set up based on the atomic composition.

[...]
// The setup for the electron numbers for each atom type and the total number of electrons.

DONE : SETUP UNITCELL Time : 0.611609 (SEC)
// The time taken to complete the unit cell setup.
```

### Doing symmetry analysis
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Doing symmetry analysis:                                           |
| We calculate the norm of 3 vectors and the angles between them,    |
| the type of Bravais lattice is given. We can judge if the unticell |
| is a primitive cell. Finally we give the point group operation for |
| this unitcell. We use the point group operations to do symmetry |
| analysis on given k-point mesh and the charge density.             |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

LATTICE VECTORS: (CARTESIAN COORDINATE: IN UNIT OF A0)
// The lattice vectors of the unit cell in Cartesian coordinates, normalized to the Bohr radius (a0).

              +7.9535                  +0                  +0
                   +0             +7.9535                  +0
                   +0                  +0             +7.9535

                       right hand lattice = 1
// Indicates the handedness of the lattice.

                                   NORM_A = 7.9535
                                   NORM_B = 7.9535
                                   NORM_C = 7.9535
// The norms (lengths) of the lattice vectors.

                           ALPHA (DEGREE) = 90
                           BETA  (DEGREE) = 90
                           GAMMA (DEGREE) = 90
// The angles between the lattice vectors in degrees.

The lattice vectors have been changed (STRU_SIMPLE.cif)
// Indicates that the lattice vectors have been updated or defined in the STRU_SIMPLE.cif file.

(for optimal symmetric configuration:)
// The following parameters describe the lattice for the optimal symmetric configuration.

                             BRAVAIS TYPE = 1
                     BRAVAIS LATTICE NAME = 01. Cubic P (simple)
                                    ibrav = 1
                                    IBRAV = 1
                                  BRAVAIS = SIMPLE CUBIC
                       LATTICE CONSTANT A = 7.9535
// The lattice constant for the simple cubic lattice.

optimized lattice volume: 503.124
// The optimized volume of the unit cell.

optimized primitive cell volume: 125.781
// The optimized volume of the primitive cell.

Original cell was built up by 4 primitive cells.
// The original unit cell is composed of 4 primitive cells.

                        ROTATION MATRICES = 48
// The number of rotation matrices for the symmetry operations.

              PURE POINT GROUP OPERATIONS = 24
                   SPACE GROUP OPERATIONS = 48
// The number of operations in the point group and the space group.

                                       C2 = 3
                                       C3 = 8
                                       C4 = 0
                                       C6 = 0
                                       S1 = 6
                                       S3 = 0
                                       S4 = 6
                                       S6 = 0
// The counts of different symmetry operations (rotations and reflections).

                              POINT GROUP = T_d
               POINT GROUP IN SPACE GROUP = O_h
// The point group and the space group for the unit cell.

DONE : SYMMETRY Time : 0.781568 (SEC)
// The time taken to complete the symmetry analysis.

```

### Setup K-points
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Setup K-points                                                     |
| We setup the k-points according to input parameters.               |
| The reduced k-points are set according to symmetry operations.     |
| We treat the spin as another set of k-points.                      |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SETUP K-POINTS
// This section describes the setup of k-points, which are used for integrating over the Brillouin zone.

                                    nspin = 1
// The number of spin dimensions; 1 for non-spin-polarized calculations.

                   Input type of k points = Monkhorst-Pack(Gamma)
// The method used for generating k-points, which is the Monkhorst-Pack scheme centered at the Gamma point.

                                   nkstot = 1
// The total number of k-points generated.

                       right hand lattice = 1
// Indicates the handedness of the lattice used for k-point generation.

(for reciprocal lattice: )
// The following parameters describe the reciprocal lattice used for k-point generation.

                             BRAVAIS TYPE = 1
                     BRAVAIS LATTICE NAME = 01. Cubic P (simple)
                                    ibrav = 1
                       right hand lattice = 1
// The Bravais lattice type and its properties for the reciprocal lattice.

(for k-lattice: )
// Parameters for the k-point lattice.

                             BRAVAIS TYPE = 1
                     BRAVAIS LATTICE NAME = 01. Cubic P (simple)
                                    ibrav = 1
                       right hand lattice = 1

                        ROTATION MATRICES = 48
// The number of rotation matrices for the symmetry operations in the reciprocal space.

                               nkstot_ibz = 1
// The number of k-points in the irreducible Brillouin zone.

K-POINTS REDUCTION ACCORDING TO SYMMETRY
// The reduction of k-points according to the symmetry operations to find the irreducible Brillouin zone.

     IBZ    DIRECT_X    DIRECT_Y    DIRECT_Z  WEIGHT  ibz2bz
       1  0.00000000  0.00000000  0.00000000  1.0000       0

                               nkstot now = 1
// The current total number of k-points after symmetry reduction.

K-POINTS DIRECT COORDINATES
// The direct coordinates of the k-points in the Brillouin zone.

 KPOINTS    DIRECT_X    DIRECT_Y    DIRECT_Z  WEIGHT
       1  0.00000000  0.00000000  0.00000000  1.0000

           k-point number in this process = 1
       minimum distributed K point number = 1

K-POINTS CARTESIAN COORDINATES
// The Cartesian coordinates of the k-points.

 KPOINTS CARTESIAN_X CARTESIAN_Y CARTESIAN_Z  WEIGHT
       1  0.00000000  0.00000000  0.00000000  2.0000

K-POINTS DIRECT COORDINATES
// The direct coordinates of the k-points, repeated for clarity.

 KPOINTS    DIRECT_X    DIRECT_Y    DIRECT_Z  WEIGHT
       1  0.00000000  0.00000000  0.00000000  2.0000

 DONE : INIT K-POINTS Time : 1.01765 (SEC)
// The time taken to initialize the k-points.
```


### Setup plane waves of wave functions
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Setup plane waves of wave functions:                               |
| Use the energy cutoff and the lattice vectors to generate the      |
| dimensions of FFT grid. The number of FFT grid on each processor   |
| is 'nrxx'. The number of plane wave basis in reciprocal space is   |
| different for charge/potential and wave functions. We also set    |
| the 'sticks' for the parallel of FFT. The number of plane wave of  |
| each k-point is 'npwk[ik]' in each processor                       |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SETUP PLANE WAVES FOR WAVE FUNCTIONS
// This section describes the setup of plane waves for the wave functions in the calculation.

     energy cutoff for wavefunc (unit:Ry) = 150
// The energy cutoff used for the wave functions in Rydberg units.

              fft grid for wave functions = [ 120, 120, 120 ]
// The dimensions of the FFT grid used specifically for wave functions.

                    number of plane waves = 105591
// The total number of plane waves used for the wave functions.

                         number of sticks = 2709
// The number of 'sticks' used in the FFT parallelization for wave functions.

PARALLEL PW FOR WAVE FUNCTIONS
// The distribution of plane waves among processors for wave function calculations.

     PROC   COLUMNS(POT)             PW
        1            339          13201
        2            338          13198
        3            338          13198
        4            338          13198
        5            339          13199
        6            339          13199
        7            339          13199
        8            339          13199
-------------- sum -------------------
        8           2709         105591

// The sum of plane waves and the total number of processors used for wave function calculations.

                            occupied bands = 136
                                   NLOCAL = 888
                                   NBANDS = 200
                                   NBANDS = 200
// Parameters indicating the number of occupied bands, local bands, and total bands in the calculation.

SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS
// The setup for nonlocal pseudopotential projectors, which are important for systems with nonlocal pseudopotentials.

SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS
// Repetition of the setup for nonlocal pseudopotential projectors.

max number of nonlocal projetors among all species is 8
// The maximum number of nonlocal projectors for all species in the system.

Warning_Memory_Consuming allocated:  TwoCenterTable: Kinetic 11.3983 MB
// A warning about the memory consumption for the kinetic part of the two-center table.

Warning_Memory_Consuming allocated:  TwoCenterTable: Overlap 11.3983 MB
// A warning about the memory consumption for the overlap part of the two-center table.

Warning_Memory_Consuming allocated:  TwoCenterTable: Nonlocal 5.26246 MB
// A warning about the memory consumption for the nonlocal part of the two-center table.

SETUP THE DIVISION OF H/S MATRIX
// The division of the Hamiltonian and overlap matrices using 2D block algorithms.

 divide the H&S matrix using 2D block algorithms.
                                     nb2d = 1
                                     nloc = 100352

// Parameters indicating the number of 2D blocks and the number of local orbitals.

-------------------------------------------
SELF-CONSISTENT
-------------------------------------------
// The beginning of the self-consistent field (SCF) calculation section.
```
## Running scf processes
### Search adjacent atoms (init)
```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|                                                                    |
| Search adjacent atoms:                                             |
| Set the adjacent atoms for each atom and set the periodic boundary |
| condition for the atoms on real space FFT grid. For k-dependent    |
| algorithm, we also need to set the sparse H and S matrix element   |
| for each atom.                                                     |
|                                                                    |
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

SETUP SEARCHING RADIUS FOR PROGRAM TO SEARCH ADJACENT ATOMS
// The program sets up a searching radius to find adjacent atoms for each atom in the unit cell.

                  longest orb rcut (Bohr) = 7
// The longest orbital radius cutoff used in the search.

   longest nonlocal projector rcut (Bohr) = 2.46
// The longest nonlocal projector radius cutoff.

              searching radius is (Bohr)) = 18.9
// The overall searching radius for finding adjacent atoms.

         searching radius unit is (Bohr)) = 1.89

SETUP EXTENDED REAL SPACE GRID FOR GRID INTEGRATION
// The setup of an extended real space grid for numerical integration.

                          real space grid = [ 120, 120, 120 ]
// The dimensions of the real space grid.

                 big cell numbers in grid = [ 40, 40, 40 ]
// The number of 'big cells' within the grid.

             meshcell numbers in big cell = [ 3, 3, 3 ]
// The number of mesh cells within each big cell.

                        extended fft grid = [ 19, 19, 19 ]
// The dimensions of the extended FFT grid.

                dimension of extened grid = [ 79, 79, 79 ]
// The actual dimensions of the extended grid.

                            UnitCellTotal = 27
// The total number of unit cells.

              Atom number in sub-FFT-grid = 24
// The number of atoms within the sub-FFT grid.

    Local orbitals number in sub-FFT-grid = 888
// The number of local orbitals within the sub-FFT grid.

                                ParaV.nnr = 441130
// A parameter related to the parallelization of the calculation.

                                     nnrg = 1878264
// The number of G-vectors used in the calculation.

Warning_Memory_Consuming allocated:  Gint::hRGint 14.5 MB
// A warning about the memory consumption for the Hamiltonian real space grid interaction.

Warning_Memory_Consuming allocated:  Gint::DMRGint 14.5 MB
// A warning about the memory consumption for the overlap real space grid interaction.

Warning_Memory_Consuming allocated:  pvpR_reduced 14.3 MB
// A warning about the memory consumption for the reduced pseudopotential real space grid.

Warning_Memory_Consuming allocated:  LOC::DM_R 14.3 MB
// A warning about the memory consumption for the local orbitals real space grid.

                                 init_chg = atomic
// The initial charge density is set to be atomic.

DONE : INIT SCF Time : 4.10934 (SEC)
// The time taken to initialize the SCF calculation.


```

### Scf iteration
```
LC AO ALGORITHM --------------- ION= 1 ELEC= 1--------------------------------

Density error is 0.0441106498563
// The error in the electron density from the current iteration of the SCF calculation.

----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2970.5084378350     -40415.8407116355    
// The Kohn-Sham energy, which is a key quantity in DFT calculations.
 E_Harris       -2971.5066745628     -40429.4224190859    
// The Harris energy, an alternative energy expression sometimes used in SCF iterations.
 E_Fermi        0.9352224901         12.7243547631        
// The Fermi energy, which represents the highest occupied energy level in the system.
----------------------------------------------------------

LCAO ALGORITHM --------------- ION=   1  ELEC=   2--------------------------------

 Density error is 0.07324240119
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2970.8761527408     -40420.8437295927    
 E_Harris       -2985.6598714921     -40621.9865422408    
 E_Fermi        0.9942982614         13.5281218666        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   3--------------------------------

 Density error is 0.0256800696305
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.1839545654     -40425.0315882625    
 E_Harris       -2973.1043954412     -40451.1605268455    
 E_Fermi        0.9712074216         13.2139548737        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   4--------------------------------

 Density error is 0.0102215946188
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2225078250     -40425.5561322694    
 E_Harris       -2974.8325144469     -40474.6727921456    
 E_Fermi        0.9535174273         12.9732701541        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   5--------------------------------

 Density error is 0.00465926410973
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2191319329     -40425.5102009013    
 E_Harris       -2972.1662486317     -40438.3963846764    
 E_Fermi        0.9547114124         12.9895151541        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   6--------------------------------

 Density error is 0.00290323197151
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2213953394     -40425.5409961268    
 E_Harris       -2971.3541089783     -40427.3466578181    
 E_Fermi        0.9597988625         13.0587334638        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   7--------------------------------

 Density error is 0.000225469262322
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221507153     -40425.5512735432    
 E_Harris       -2971.1959132112     -40425.1942939861    
 E_Fermi        0.9575455903         13.0280761230        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   8--------------------------------

 Density error is 0.000105298036339
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221543072     -40425.5513224128    
 E_Harris       -2971.2345871224     -40425.7204795427    
 E_Fermi        0.9573864389         13.0259107570        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=   9--------------------------------

 Density error is 1.73342865404e-05
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545084     -40425.5513251506    
 E_Harris       -2971.2284300744     -40425.6367086066    
 E_Fermi        0.9573421242         13.0253078248        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  10--------------------------------

 Density error is 2.44145213497e-05
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545026     -40425.5513250720    
 E_Harris       -2971.2237593087     -40425.5731595783    
 E_Fermi        0.9573530336         13.0254562547        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  11--------------------------------

 Density error is 1.77043483202e-05
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545218     -40425.5513253324    
 E_Harris       -2971.2233998438     -40425.5682688075    
 E_Fermi        0.9573601195         13.0255526627        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  12--------------------------------

 Density error is 1.08567789394e-06
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545365     -40425.5513255324    
 E_Harris       -2971.2226767646     -40425.5584308108    
 E_Fermi        0.9573493003         13.0254054601        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  13--------------------------------

 Density error is 1.39646278617e-05
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545147     -40425.5513252362    
 E_Harris       -2971.2222832103     -40425.5530762294    
 E_Fermi        0.9573484591         13.0253940153        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  14--------------------------------

 Density error is 3.71116003478e-07
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545386     -40425.5513255615    
 E_Harris       -2971.2222838082     -40425.5530843652    
 E_Fermi        0.9573492888         13.0254053046        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  15--------------------------------

 Density error is 2.46641928838e-07
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545386     -40425.5513255614    
 E_Harris       -2971.2221840706     -40425.5517273647    
 E_Fermi        0.9573494485         13.0254074765        
----------------------------------------------------------


 LCAO ALGORITHM --------------- ION=   1  ELEC=  16--------------------------------

 Density error is 2.05528957543e-07
----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545386     -40425.5513255622    
 E_Harris       -2971.2221699838     -40425.5515357039    
 E_Fermi        0.9573494392         13.0254073506        
----------------------------------------------------------

LC AO ALGORITHM --------------- ION= 1 ELEC= 17--------------------------------

Density error is 6.53205775131e-08
// The error in the electron density has become very small, indicating near-convergence of the SCF cycle.
// Additional iterations of the SCF calculation, showing the gradual convergence of the density error and energies.
```
### Scf results
```

----------------------------------------------------------
     Energy           Rydberg                 eV          
----------------------------------------------------------
 E_KohnSham     -2971.2221545387     -40425.5513255626    
 E_KS(sigma->0) -2971.1513960918     -40424.5886075033    
// The Kohn-Sham energy with the sigma component of the density matrix going to zero, an indicator of the kinetic energy contribution.
 E_Harris       -2971.2221622072     -40425.5514298979    
 E_band         -562.7029473746      -7655.9663656893     
// The band energy, associated with the dispersion of the energy bands.
 E_one_elec     -1460.3005783359     -19868.4086580639    
// The one-electron contribution to the energy.
 E_Hartree      672.7815109284       9153.6620576760      
// The Hartree energy, representing the classical electrostatic interaction between electron densities.
 E_xc           -486.3000633355      -6616.4517991236     
// The exchange-correlation energy, accounting for quantum mechanical effects.
 E_Ewald        -1697.2615069019     -23092.4274899325    
// The Ewald energy, the contribution to the total energy from the Ewald summation for long-range interactions.
 E_entropy(-TS) -0.1415168938        -1.9254361186        
// The entropy contribution to the free energy (negative of the product of temperature and entropy).
 E_descf        0.0000000000         0.0000000000         
// The contribution to the energy from the density-cutoff term in the pseudopotential.
 E_exx          0.0000000000         0.0000000000         
// The exact exchange energy, if applicable (usually zero in standard DFT calculations).
 E_Fermi        0.9573493498         13.0254061340        
// The Fermi energy at this stage of the calculation.

charge density convergence is achieved
// A statement indicating that the charge density has converged to the desired accuracy.

final etot is -40425.551326 eV
// The final total energy of the system after the SCF calculation has fully converged.

EFERMI = 13.025406134 eV
// The final Fermi energy of the system.

STATE ENERGY(eV) AND OCCUPATIONS    NSPIN == 1
// The energies of the electronic states and their occupations for each k-point.

1/1 kpoint (Cartesian) = 0.0000 0.0000 0.0000 (13201 pws)
// The k-point at the center of the Brillouin zone and the number of plane waves associated with it.
       1       -90.2506        2.00000
       2       -90.2494        2.00000
       3       -90.2494        2.00000
       4       -90.2494        2.00000
       5       -90.2494        2.00000
       6       -90.2494        2.00000
       7       -90.2494        2.00000
       8       -90.2482        2.00000
       9       -90.2482        2.00000
      10       -90.2482        2.00000
      11       -90.2482        2.00000
      12       -90.2482        2.00000
      13       -90.2482        2.00000
      14       -90.2482        2.00000
      15       -90.2482        2.00000
      16       -90.2482        2.00000
      17       -51.8782        2.00000
      18       -51.8681        2.00000
      19       -51.8681        2.00000
      20       -51.8681        2.00000
      21       -51.8681        2.00000
      22       -51.8681        2.00000
      23       -51.8681        2.00000
      24       -51.8583        2.00000
      25       -51.8583        2.00000
      26       -51.8583        2.00000
      27       -51.8580        2.00000
      28       -51.8580        2.00000
      29       -51.8580        2.00000
      30       -51.8580        2.00000
      31       -51.8580        2.00000
      32       -51.8580        2.00000
      33       -51.7617        2.00000
      34       -51.7617        2.00000
      35       -51.7612        2.00000
      36       -51.7612        2.00000
      37       -51.7612        2.00000
      38       -51.7612        2.00000
      39       -51.7612        2.00000
      40       -51.7612        2.00000
      41       -51.7604        2.00000
      42       -51.7604        2.00000
      43       -51.7604        2.00000
      44       -51.7591        2.00000
      45       -51.7591        2.00000
      46       -51.7591        2.00000
      47       -51.7591        2.00000
      48       -51.7591        2.00000
      49       -51.7591        2.00000
      50       -51.7578        2.00000
      51       -51.7578        2.00000
      52       -51.7578        2.00000
      53       -51.7578        2.00000
      54       -51.7578        2.00000
      55       -51.7578        2.00000
      56       -51.7556        2.00000
      57       -51.7556        2.00000
      58       -51.7556        2.00000
      59       -51.7556        2.00000
      60       -51.7556        2.00000
      61       -51.7556        2.00000
      62       -51.7539        2.00000
      63       -51.7539        2.00000
      64       -51.7539        2.00000
      65       -21.8364        2.00000
      66       -21.7006        2.00000
      67       -21.7006        2.00000
      68       -21.7006        2.00000
      69       -21.7006        2.00000
      70       -21.7006        2.00000
      71       -21.7006        2.00000
      72       -21.5621        2.00000
      73       -5.23097        2.00000
      74       -5.23097        2.00000
      75       -5.23097        2.00000
      76       -5.23097        2.00000
      77       -5.23097        2.00000
      78       -5.23097        2.00000
      79       -5.02839        2.00000
      80       -5.02839        2.00000
      81       -5.02839        2.00000
      82       -4.76611        2.00000
      83       -4.76611        2.00000
      84       -4.76611        2.00000
      85       -4.76611        2.00000
      86       -4.76611        2.00000
      87       -4.76611        2.00000
      88       -4.42730        2.00000
      89       -4.42730        2.00000
      90       -4.42730        2.00000
      91       -4.21659        2.00000
      92       -4.21659        2.00000
      93       -4.21659        2.00000
      94       -4.21659        2.00000
      95       -4.21659        2.00000
      96       -4.21659        2.00000
      97        3.53302        2.00000
      98        5.39913        2.00000
      99        5.39913        2.00000
     100        5.39913        2.00000
     101        5.39913        2.00000
     102        5.39913        2.00000
     103        5.39913        2.00000
     104        8.74954        2.00000
     105        9.08287        2.00000
     106        9.08287        2.00000
     107        9.08287        2.00000
     108        9.08287        2.00000
     109        9.08287        2.00000
     110        9.08287        2.00000
     111        9.57721        2.00000
     112        9.57721        2.00000
     113        9.57721        2.00000
     114        9.77291        2.00000
     115        9.77291        2.00000
     116        9.77291        2.00000
     117        9.77291        2.00000
     118        9.77291        2.00000
     119        9.77291        2.00000
     120        10.1656        2.00000
     121        10.1656        2.00000
     122        10.4293        2.00000
     123        10.4293        2.00000
     124        10.4293        2.00000
     125        12.1630        1.99999
     126        12.1630        1.99999
     127        12.1630        1.99999
     128        12.1630        1.99999
     129        12.1630        1.99999
     130        12.1630        1.99999
     131        12.9733        1.21344
     132        12.9733        1.21344
     133        12.9733        1.21344
     134        12.9733        1.21344
     135        12.9733        1.21344
     136        12.9733        1.21344
     137        13.1019       0.690941
     138        13.1019       0.690941
     139        13.1019       0.690941
     140        13.1019       0.690941
     141        13.1019       0.690940
     142        13.1019       0.690940
     143        13.1733       0.442128
     144        13.4597      0.0240081
     145        13.4597      0.0240081
     146        13.4597      0.0240081
     147        13.5463     0.00678529
     148        13.5463     0.00678529
     149        13.5463     0.00678528
     150        13.5463     0.00678527
     151        13.5463     0.00678527
     152        13.5463     0.00678526
     153        13.6052     0.00258518
     154        13.6052     0.00258518
     155        13.6052     0.00258518
     156        13.6052     0.00258517
     157        13.6052     0.00258517
     158        13.6052     0.00258517
     159        13.7097    0.000376289
     160        13.7097    0.000376289
     161        13.7097    0.000376289
     162        13.7131    0.000351588
     163        13.7131    0.000351587
     164        13.7131    0.000351586
     165        13.7131    0.000351586
     166        13.7131    0.000351586
     167        13.7131    0.000351585
     168        13.8238    3.33179e-05
     169        13.8238    3.33179e-05
     170        13.8238    3.33179e-05
     171        13.8978    5.78882e-06
     172        13.8978    5.78881e-06
     173        13.8978    5.78881e-06
     174        13.8978    5.78879e-06
     175        13.8978    5.78879e-06
     176        13.8978    5.78879e-06
     177        13.9176    3.54072e-06
     178        13.9176    3.54072e-06
     179        13.9176    3.54072e-06
     180        13.9299    2.58893e-06
     181        14.1919    1.33980e-09
     182        14.1919    1.33980e-09
     183        14.1919    1.33980e-09
     184        14.1919    1.33980e-09
     185        14.1919    1.33980e-09
     186        14.1919    1.33980e-09
     187        14.3213    1.64246e-11
     188        14.3213    1.64246e-11
     189        14.3213    1.64246e-11
     190        14.3213    1.64246e-11
     191        14.3213    1.64246e-11
     192        14.3213    1.64246e-11
     193        14.5304    5.21805e-15
     194        14.5304    5.21805e-15
     195        14.5304    5.21805e-15
     196        14.5304    5.21805e-15
     197        14.5304    5.21805e-15
     198        14.5304    5.21805e-15
     199        14.5998    3.33067e-16
     200        14.7810        0.00000
// The list of state energies and their occupations for the converged calculation.

Warning_Memory_Consuming allocated:  Force::dS_K 10.0967 MB
// A warning indicating the memory allocated for the calculation of forces, specifically the kinetic contribution.

Warning_Memory_Consuming allocated:  Stress::dHr 10.0967 MB
// A warning about the memory consumption for the calculation of stress, specifically the Hellmann-Feynman contribution.

Warning_Memory_Consuming allocated:  Stress::dSR 20.1933 MB
// A warning about the memory consumption for the calculation of stress, specifically the nonlocal pseudopotential contribution.

Warning_Memory_Consuming allocated:  Force::dTVNL 10.0967 MB
// A warning about the memory consumption for the calculation of forces, specifically the nonlocal pseudopotential contribution.

correction force for each atom along direction 1 is 2.90404e-14
correction force for each atom along direction 2 is 2.49665e-14
correction force for each atom along direction 3 is -3.50742e-14
// The correction forces applied to each atom in the unit cell along the three principal directions.
```
## Results summary
### TOTAL-FORCE
```
------------------------------------------------------------------------------------------
 TOTAL-FORCE (eV/Angstrom)                                                                
------------------------------------------------------------------------------------------
// The total force experienced by each atom in the unit cell, expressed in electronvolts per angstrom.
                       Ce1        -0.0000024495        -0.0000024495         0.0000024495 
                       Ce2        -0.0000024495        -0.0000024495        -0.0000024495 
                       Ce3        -0.0000024495         0.0000024495        -0.0000024495 
                       Ce4        -0.0000024495         0.0000024495         0.0000024495 
                       Ce5         0.0000024495        -0.0000024495        -0.0000024495 
                       Ce6         0.0000024495        -0.0000024495         0.0000024495 
                       Ce7         0.0000024495         0.0000024495         0.0000024495 
                       Ce8         0.0000024495         0.0000024495        -0.0000024495 
                       Al1         0.0000000000         0.0000042432        -0.0000042432 
                       Al2         0.0000000000         0.0000000000         0.0000000000 
                       Al3        -0.0000026880         0.0000026880         0.0000015552 
                       Al4        -0.0000042432         0.0000000000        -0.0000042432 
                       Al5         0.0000000000        -0.0000042432         0.0000042432 
                       Al6         0.0000000000         0.0000000000         0.0000000000 
                       Al7        -0.0000042432        -0.0000042432         0.0000000000 
                       Al8        -0.0000026880         0.0000015552         0.0000026880 
                       Al9        -0.0000015552         0.0000026880         0.0000026880 
                      Al10         0.0000000000         0.0000000000         0.0000000000 
                      Al11         0.0000042432         0.0000042432         0.0000000000 
                      Al12         0.0000042432         0.0000000000         0.0000042432 
                      Al13         0.0000015552        -0.0000026880        -0.0000026880 
                      Al14         0.0000000000         0.0000000000         0.0000000000 
                      Al15         0.0000026880        -0.0000026880        -0.0000015552 
                      Al16         0.0000026880        -0.0000015552        -0.0000026880 
------------------------------------------------------------------------------------------
// The table listing the total forces on each atom, indicating whether the system is close to a force minimum (stable configuration).
```

### TOTAL-STRESS
```
----------------------------------------------------------------
 TOTAL-STRESS (KBAR)                                            
----------------------------------------------------------------
// The total stress experienced by the unit cell, a measure of the pressure in different directions.

        20.9280307419        -0.0000000000         0.0000000000 
        -0.0000000000        20.9280307419         0.0000000000 
         0.0000000000         0.0000000000        20.9280307419 

TOTAL-PRESSURE: 20.928031 KBAR
// The total pressure exerted on the system, calculated from the stress tensor.
```

### FINAL_ETOT_IS
```
--------------------------------------------
 !FINAL_ETOT_IS -40425.5513255625628517 eV
 --------------------------------------------
// The final total energy of the system after the SCF calculation, expressed in electron volts.
```



### TIME STATISTICS
```
TIME STATISTICS
--------------------------------------------------------------------------------
      CLASS_NAME                  NAME            TIME/s  CALLS   AVG/s  PER/%
--------------------------------------------------------------------------------
// A summary of the time spent in different parts of the calculation, providing insights into the efficiency and potential bottlenecks.
                       total                      500.80 11       45.53  100.00 
 Driver                reading                    0.09   1        0.09   0.02   
 Input                 Init                       0.08   1        0.08   0.02   
 Input_Conv            Convert                    0.00   1        0.00   0.00   
 Driver                driver_line                500.71 1        500.71 99.98  
 UnitCell              check_tau                  0.00   1        0.00   0.00   
 ESolver_KS_LCAO       before_all_runners         2.10   1        2.10   0.42   
 PW_Basis_Sup          setuptransform             0.03   1        0.03   0.01   
 PW_Basis_Sup          distributeg                0.02   1        0.02   0.00   
 mymath                heapsort                   0.02   2016     0.00   0.00   
 Symmetry              analy_sys                  0.17   1        0.17   0.03   
 PW_Basis_K            setuptransform             0.01   1        0.01   0.00   
 PW_Basis_K            distributeg                0.01   1        0.01   0.00   
 PW_Basis              setup_struc_factor         0.13   1        0.13   0.03   
 NOrbital_Lm           extra_uniform              0.18   26       0.01   0.04   
 Mathzone_Add1         SplineD2                   0.00   26       0.00   0.00   
 Mathzone_Add1         Cubic_Spline_Interpolation 0.01   26       0.00   0.00   
 Mathzone_Add1         Uni_Deriv_Phi              0.16   26       0.01   0.03   
 ppcell_vl             init_vloc                  0.12   1        0.12   0.02   
 Ions                  opt_ions                   498.27 1        498.27 99.50  
 ESolver_KS_LCAO       runner                     373.62 1        373.62 74.60  
 ESolver_KS_LCAO       before_scf                 1.64   1        1.64   0.33   
 ESolver_KS_LCAO       beforesolver               0.76   1        0.76   0.15   
 ESolver_KS_LCAO       set_matrix_grid            0.58   1        0.58   0.12   
 atom_arrange          search                     0.00   1        0.00   0.00   
 Grid_Technique        init                       0.52   1        0.52   0.10   
 Grid_BigCell          grid_expansion_index       0.03   2        0.01   0.01   
 Record_adj            for_2d                     0.02   1        0.02   0.00   
 Grid_Driver           Find_atom                  0.03   792      0.00   0.01   
 LCAO_domain           grid_prepare               0.00   1        0.00   0.00   
 Veff                  initialize_HR              0.00   1        0.00   0.00   
 OverlapNew            initialize_SR              0.00   1        0.00   0.00   
 EkineticNew           initialize_HR              0.00   1        0.00   0.00   
 NonlocalNew           initialize_HR              0.00   1        0.00   0.00   
 Charge                set_rho_core               0.14   1        0.14   0.03   
 PW_Basis_Sup          recip2real                 3.44   122      0.03   0.69   
 PW_Basis_Sup          gathers_scatterp           1.72   122      0.01   0.34   
 Charge                atomic_rho                 0.14   1        0.14   0.03   
 Potential             init_pot                   0.30   1        0.30   0.06   
 Potential             update_from_charge         11.23  18       0.62   2.24   
 Potential             cal_fixed_v                0.01   1        0.01   0.00   
 PotLocal              cal_fixed_v                0.01   1        0.01   0.00   
 Potential             cal_v_eff                  11.19  18       0.62   2.23   
 H_Hartree_pw          v_hartree                  1.32   18       0.07   0.26   
 PW_Basis_Sup          real2recip                 3.89   144      0.03   0.78   
 PW_Basis_Sup          gatherp_scatters           1.92   144      0.01   0.38   
 PotXC                 cal_v_eff                  9.81   18       0.55   1.96   
 XC_Functional         v_xc                       11.76  20       0.59   2.35   
 Potential             interpolate_vrs            0.03   18       0.00   0.01   
 Symmetry              rhog_symmetry              12.41  18       0.69   2.48   
 Symmetry              group fft grids            2.64   18       0.15   0.53   
 H_Ewald_pw            compute_ewald              0.01   1        0.01   0.00   
 Charge_Mixing         init_mixing                0.00   1        0.00   0.00   
 HSolverLCAO           solve                      343.26 17       20.19  68.54  
 HamiltLCAO            updateHk                   120.31 17       7.08   24.02  
 OperatorLCAO          init                       119.59 51       2.34   23.88  
 Veff                  contributeHR               118.31 17       6.96   23.62  
 Gint_interface        cal_gint                   320.96 35       9.17   64.09  
 Gint_interface        cal_gint_vlocal            114.33 17       6.73   22.83  
 Gint_Tools            cal_psir_ylm               82.79  272000   0.00   16.53  
 Gint_k                transfer_pvpR              3.98   17       0.23   0.79   
 OverlapNew            calculate_SR               0.67   1        0.67   0.13   
 OverlapNew            contributeHk               0.06   17       0.00   0.01   
 EkineticNew           contributeHR               0.67   17       0.04   0.13   
 EkineticNew           calculate_HR               0.67   1        0.67   0.13   
 NonlocalNew           contributeHR               0.49   17       0.03   0.10   
 NonlocalNew           calculate_HR               0.45   1        0.45   0.09   
 OperatorLCAO          contributeHk               0.09   17       0.01   0.02   
 HSolverLCAO           hamiltSolvePsiK            90.50  17       5.32   18.07  
 ElecStateLCAO         psiToRho                   132.45 17       7.79   26.45  
 elecstate             cal_dm                     2.02   18       0.11   0.40   
 psiMulPsiMpi          pdgemm                     2.01   18       0.11   0.40   
 DensityMatrix         cal_DMR                    0.13   18       0.01   0.03   
 Gint                  transfer_DMR               2.98   17       0.18   0.60   
 Gint_interface        cal_gint_rho               127.22 17       7.48   25.40  
 Charge_Mixing         get_drho                   0.07   17       0.00   0.01   
 Charge                mix_rho                    0.74   16       0.05   0.15   
 Charge                Pulay_mixing               0.45   16       0.03   0.09   
 ESolver_KS_LCAO       out_deepks_labels          0.00   1        0.00   0.00   
 LCAO_Deepks_Interface out_deepks_labels          0.00   1        0.00   0.00   
 ESolver_KS_LCAO       cal_force                  124.65 1        124.65 24.89  
 Force_Stress_LCAO     getForceStress             124.65 1        124.65 24.89  
 Forces                cal_force_loc              0.62   1        0.62   0.12   
 Forces                cal_force_ew               0.53   1        0.53   0.11   
 Forces                cal_force_cc               1.67   1        1.67   0.33   
 Forces                cal_force_scc              1.20   1        1.20   0.24   
 Stress_Func           stress_loc                 2.00   1        2.00   0.40   
 Stress_Func           stress_har                 0.11   1        0.11   0.02   
 Stress_Func           stress_ewa                 0.49   1        0.49   0.10   
 Stress_Func           stress_cc                  2.99   1        2.99   0.60   
 Stress_Func           stress_gga                 0.64   1        0.64   0.13   
 Force_LCAO            ftable                     114.38 1        114.38 22.84  
 Force_LCAO            allocate                   0.00   1        0.00   0.00   
 LCAO_domain           build_ST_new               12.21  2        6.10   2.44   
 LCAO_domain           vnl_mu_new                 6.91   1        6.91   1.38   
 Force_LCAO_k          allocate_k                 0.00   1        0.00   0.00   
 Force_LCAO            cal_fedm                   1.59   1        1.59   0.32   
 Force_LCAO_k          cal_edm_2d                 0.00   1        0.00   0.00   
 Force_LCAO            cal_ftvnl_dphi             0.07   1        0.07   0.01   
 Force_LCAO            cal_fvl_dphi               79.41  1        79.41  15.86  
 Gint_interface        cal_gint_force             79.41  1        79.41  15.86  
 Gint_Tools            cal_dpsir_ylm              40.76  8000     0.01   8.14   
 Gint_Tools            cal_dpsirr_ylm             10.64  8000     0.00   2.12   
 Force_LCAO            cal_fvnl_dbeta             11.73  1        11.73  2.34   
 ESolver_KS_LCAO       cal_stress                 0.00   1        0.00   0.00   
 ESolver_KS_LCAO       after_all_runners          0.02   1        0.02   0.00   
 ModuleIO              write_istate_info          0.02   1        0.02   0.00   
--------------------------------------------------------------------------------
// The breakdown of time statistics for various components of the ABACUS calculation.
```

### MEMORY STATISTICS
```

NAME-------------------------|MEMORY(MB)--------
// A summary of the memory consumption by different parts of the calculation.
                         total      1473.1712
                   Stress::dSR       164.4884
                  Gint::hRGint       116.0392
                 Gint::DMRGint       115.7267
                  pvpR_reduced       114.6401
                     LOC::DM_R       114.6401
       TwoCenterTable: Kinetic        91.1862
       TwoCenterTable: Overlap        91.1862
                   Force::dS_K        82.2442
                   Stress::dHr        82.2442
                  Force::dTVNL        82.2442
                     FFT::grid        54.4922
      TwoCenterTable: Nonlocal        42.0997
                HamiltLCAO::hR        28.1311
            DensityMatrix::DMR        28.1311
                  SF::strucFac        25.7154
               LOC::wfc_k_grid        21.6797
               GT::ind_bigcell        15.0464
         GT::in_this_processor        15.0464
              GT::index2normal        15.0464
               GT::index2ucell        15.0464
                HamiltLCAO::sR        14.9202
                      Chg::rho        13.1836
                 Chg::rho_save        13.1836
                 Chg::rho_core        13.1836
                 Pot::veff_fix        13.1836
                     Pot::veff        13.1836
              Pot::veff_smooth        13.1836
            DensityMatrix::DMK        12.0322
                     Chg::rhog         6.4288
                Chg::rhog_save         6.4288
                Chg::rhog_core         6.4288
                  meshball_pos         6.2642
      GT::bigcell_on_processor         3.7616
                RealGauntTable         3.6165
                GT::which_atom         3.1321
             GT::which_bigcell         3.1321
            GT::which_unitcell         3.1321
                  PW_B_K::gcar         2.4168
                  SF::eigts123         2.1097
                Grid::AtomLink         1.2822
                    index_ball         1.0440
 -------------   < 1.0 MB has been ignored ----------------
 ----------------------------------------------------------
// The list of memory allocations for various components, indicating the memory usage efficiency.
```

### Start and end times
```
Start  Time  : Thu Jul 18 11:34:56 2024
Finish Time  : Thu Jul 18 11:43:20 2024
Total  Time  : 0 h 8 mins 24 secs
// The start and end times of the calculation, along with the total duration.
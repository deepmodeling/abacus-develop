# Basic electronic structure calculation with PW basis set

[back to main page](../../README.md)

For this example, the input files are:
- INPUT:
    ```
    INPUT_PARAMETERS

    ntype 1
    nbands 4
    
    ecutwfc 60

    basis_type pw
    suffix Si2_diamond
    symmetry 1
    scf_nmax 60
    scf_thr 1.0e-9
    out_chg 1
    ```

    The meanings of the above parameters are:

    - ntype

        how many types of elements in the unit cell

    - nbands

        the number of bands to be calculated

    - ecutwfc
        
        the plane-wave energy cutoff for the wave function expansion (UNIT: Rydberg)

    - basis_type
    
        The type of basis. The default value is pw, meaning that the plane wave set is used.

        Attention: *This is a very important parameter in ABACUS. For more information, please refer to the [instruction](../input-main.md)*.

    - suffix

        Suffix of output directory. In this example the name of the output directory will be OUT.Si2_diamond. The default value is ABACUS.
      
    - symmetry

        Use symmetry(=1) or not(=0) in the calculation. The default value is 0.

    - scf_nmax
        The maximal iteration number for electronic-structure calculations.
    - scf_thr
        Tolerance of the difference of charge density, below which the self-consistent calculation is considered to be converged.
    - out_chg
        Print out the charge density(=1) or not(=0).

    A complete list of INPUT keyewords can be found in the [instruction](../input-main.md).

- STRU
    ```
    ATOMIC_SPECIES
    Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file
    NUMERICAL_ORBITAL
    Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
    LATTICE_CONSTANT
    10.2 // lattice scaling factor (Bohr)
    LATTICE_VECTORS
    0.5 0.5 0.0 // latvec1
    0.5 0.0 0.5 // latvec2
    0.0 0.5 0.5 // latvec3
    ATOMIC_POSITIONS
    Direct //Cartesian or Direct coordinate.
    Si // Element type
    0.0 // magnetism
    2 // number of atoms
    0.00 0.00 0.00 0 0 0
    0.25 0.25 0.25 1 1 1
    ```
    We provide an [instruction](../input-stru.md) on the speicifications of the STRU file.

- KPT
    ```
    K_POINTS
    0
    Gamma
    4 4 4 0.0 0.0 0.0
    ```
    More information on the KPT file can be found in this [short instruction](../input-kpt).


The following typical output information will be printed to the screen:
```
 *********************************************************
 *                                                       *
 *                  WELCOME TO ABACUS                    *
 *                                                       *
 *            'Atomic-orbital Based Ab-initio            *
 *                  Computation at UStc'                 *
 *                                                       *
 *          Website: http://abacus.ustc.edu.cn/          *
 *                                                       *
 *********************************************************
 Fri Jul  2 11:06:27 2021
 MAKE THE DIR         : OUT.Si2_diamond/
 DONE(0.0892439  SEC) : SETUP UNITCELL
 DONE(0.134056   SEC) : SYMMETRY
 DONE(0.136979   SEC) : INIT K-POINTS
 ---------------------------------------------------------
 Self-consistent calculations for electrons
 ---------------------------------------------------------
 SPIN    KPOINTS         PROCESSORS  
 1       8               1           
 ---------------------------------------------------------
 Use plane wave basis
 ---------------------------------------------------------
 ELEMENT NATOM       XC          
 Si      2           PZ-LDA
 ---------------------------------------------------------
 Initial plane wave basis and FFT box
 ---------------------------------------------------------
 DONE(0.420921   SEC) : INIT PLANEWAVE
 UNIFORM GRID DIM     : 36 * 36 * 36
 UNIFORM GRID DIM(BIG): 36 * 36 * 36
 MEMORY FOR PSI (MB)  : 1.03516
 DONE(0.434272   SEC) : LOCAL POTENTIAL
 DONE(0.45294    SEC) : NON-LOCAL POTENTIAL
 START POTENTIAL      : atomic
 DONE(0.463443   SEC) : INIT POTENTIAL
 DONE(0.63761    SEC) : INIT BASIS
 -------------------------------------------
 SELF-CONSISTENT : 
 -------------------------------------------
 ITER   ETOT(eV)       EDIFF(eV)      SCF_THR      CG_ITER    TIME(S)    
 CG1    -2.154524e+02  0.000000e+00   6.855e-02  1.000e+00  1.290e+00  
 CG2    -2.154992e+02  -4.673475e-02  2.378e-03  2.000e+00  8.400e-01  
 CG3    -2.155050e+02  -5.882715e-03  8.220e-05  2.594e+00  9.000e-01  
 CG4    -2.155054e+02  -3.840554e-04  7.415e-06  2.938e+00  1.090e+00  
 CG5    -2.155054e+02  -1.762864e-05  2.296e-07  2.875e+00  8.600e-01  
 CG6    -2.155054e+02  -7.438412e-07  9.511e-10  2.938e+00  9.200e-01  

  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
 B Run_pw              plane_wave_line     2.2255         1         2.2       99        %
 B PW_Basis            gen_pw              0.28373        1         0.28      13        %
 C Cell_PW             opt_cells_pw        1.8114         1         1.8       80        %
 X FFT                 FFT3D               0.89387        1876      0.00048   40        %
 C wavefunc            wfcinit             0.13829        1         0.14      6.1       %
 G Hamilt_PW           diagH_subspace      0.31748        48        0.0066    14        %
 H Hamilt_PW           h_psi               0.95722        635       0.0015    42        %
 H Hamilt_PW           vloc                0.89713        635       0.0014    40        %
 C Ions                opt_ions_pw         1.5966         1         1.6       71        %
 D Electrons           self_consistent     1.5966         1         1.6       71        %
 E Electrons           c_bands             0.65341        4         0.16      29        %
 F Hamilt              diagH_pw            0.89723        56        0.016     40        %
 G Diago_CG            diag                0.70296        56        0.013     31        %
 E electrons           c_bands             0.57657        4         0.14      26        %
 ----------------------------------------------------------------------------------------

 START  Time  : Fri Jul  2 11:06:27 2021
 FINISH Time  : Fri Jul  2 11:06:29 2021
 TOTAL  Time  : 2
 SEE INFORMATION IN : OUT.Si2_diamond/

```

[back to top](#Basic-electronic-structure-calculation-with-PW-basis-set)
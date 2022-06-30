# Stochastic DFT and mix stochastic-deterministic DFT

[back to main page](../../README.md)

For sto-scf example (examples/04a_Si_diamond_pw_stoscf), the input files are:
- INPUT:
    ```
    INPUT_PARAMETERS

    #Parameters	(General)
    calculation     sto-scf
    pseudo_dir		./	
    ntype			1	
    nbands			4
    nbands_sto      64
    nche_sto        100
    method_sto      1

    #Parameters (Accuracy)
    ecutwfc			50
    scf_nmax		20
    symmetry		1

    #Parameters (Smearing)
    smearing_method     fd #We can only use fd in "sto-scf/md" calculation
    smearing_sigma      0.6
    ```

    The meanings of the above parameters are:

    - nbands

        number of KS orbitals

    - nbands_sto

        number of stochastic orbitals

    - nche_sto

        Chebyshev expansion orders

    - method_sto

        method to do stochastic calculations
    

    A complete list of INPUT keyewords can be found in the [instruction](../input-main.md).

- STRU
    ```
    ATOMIC_SPECIES
    Si 1.000 Si.pz-vbc.UPF 	#Element, Mass, Pseudopotential
    
    LATTICE_CONSTANT
    10.2  			#Lattice constant
    
    LATTICE_VECTORS
    0.5 0.5 0.0 		#Lattice vector 1
    0.5 0.0 0.5 		#Lattice vector 2
    0.0 0.5 0.5 		#Lattice vector 3
    
    ATOMIC_POSITIONS
    Cartesian 		#Cartesian(Unit is LATTICE_CONSTANT)
    Si 			#Name of element	
    0.0			#Magnetic for this element.
    2			#Number of atoms
    0.00 0.00 0.00 1 1 1
    0.25 0.25 0.25 1 1 1
    ```
    We provide an [instruction](../input-stru.md) on the speicifications of the STRU file.

- KPT
    ```
    K_POINTS
    0
    Gamma
    2 2 2 0 0 0
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
 Wed Jun 29 17:07:44 2022
 MAKE THE DIR         : OUT.ABACUS/
 UNIFORM GRID DIM     : 36 * 36 * 36
 UNIFORM GRID DIM(BIG): 36 * 36 * 36
 DONE(0.041765   SEC) : SETUP UNITCELL
 DONE(0.0794668  SEC) : SYMMETRY
 DONE(0.0900106  SEC) : INIT K-POINTS
 ---------------------------------------------------------
 Self-consistent calculations for electrons
 ---------------------------------------------------------
 SPIN    KPOINTS         PROCESSORS  
 1       3               4           
 ---------------------------------------------------------
 Use plane wave basis
 ---------------------------------------------------------
 ELEMENT NATOM       XC          
 Si      2           
 ---------------------------------------------------------
 Initial plane wave basis and FFT box
 ---------------------------------------------------------
 DONE(0.0919447  SEC) : INIT PLANEWAVE
 DONE(0.0927064  SEC) : INIT CHARGE
 MEMORY FOR PSI (MB)  : 0.0736084
 DONE(0.0937741  SEC) : LOCAL POTENTIAL
 DONE(0.101361   SEC) : NON-LOCAL POTENTIAL
 START POTENTIAL      : atomic
 DONE(0.105865   SEC) : INIT POTENTIAL
 DONE(0.12326    SEC) : INIT BASIS
 -------------------------------------------
 SELF-CONSISTENT : 
 -------------------------------------------
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
 CG1    -2.988673e+02  0.000000e+00   8.370e-03  8.749e+00  
 CG2    -2.988733e+02  -5.935642e-03  5.393e-04  8.713e+00  
 CG3    -2.988722e+02  1.108353e-03   6.628e-06  8.304e+00  
 CG4    -2.988727e+02  -5.074167e-04  6.481e-08  8.684e+00  
 CG5    -2.988727e+02  7.521919e-06   2.904e-09  8.569e+00  
 CG6    -2.988727e+02  -7.194917e-06  4.086e-11  8.483e+00  

  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                       total               51.638         19        2.7       1e+02     %
   Run_pw              plane_wave_line     51.626         1         52        1e+02     %
   PW_Basis_K          recip2real          23.632         231192    0.0001    46        %
   PW_Basis_K          gathers_scatterp    8.462          231192    3.7e-05   16        %
   PW_Basis_K          real2recip          21.244         229968    9.2e-05   41        %
   PW_Basis_K          gatherp_scatters    6.5738         229968    2.9e-05   13        %
   Cell_PW             opt_cells_pw        51.508         1         52        1e+02     %
   Ions                opt_ions_pw         51.507         1         52        1e+02     %
   ESolver_SDFT_PW     Run                 51.507         1         52        1e+02     %
   HSolverPW_SDFT      solve               51.481         6         8.6       1e+02     %
   Stochastic_hchi     hchi_reciprocal     50.717         5148      0.0099    98        %
   Stochastic_hchi     vloc                48.819         5148      0.0095    95        %
   Stochastic_hchi     vnl                 1.1863         5148      0.00023   2.3       %
   Stochastic_Iter     calPn               25.41          18        1.4       49        %
   Stochastic_Iter     sum_stoband         25.467         6         4.2       49        %
 ----------------------------------------------------------------------------------------

 START  Time  : Wed Jun 29 17:07:44 2022
 FINISH Time  : Wed Jun 29 17:08:35 2022
 TOTAL  Time  : 51
 SEE INFORMATION IN : OUT.ABACUS/

```

For sto-md example (examples/03a_Sn64_lcao_md), the input files are:
- INPUT:
    ```
    INPUT_PARAMETERS

    #Parameters	(General)
    calculation     sto-md
    pseudo_dir		./	
    ntype			1	
    nbands			0
    nbands_sto      64
    nche_sto        20
    method_sto      2

    #Parameters (Accuracy)
    ecutwfc			50
    scf_nmax		20
    scf_thr         1e-6
    symmetry		1

    #Parameters (Smearing)
    smearing_method     fd
    smearing_sigma      7.34986072

    #Parameters (MD)
    md_tfirst      1160400
    md_dt          0.2
    md_nstep       10
    ```
    A complete list of INPUT keyewords can be found in the [instruction](../input-main.md).

- STRU
    ```
    ATOMIC_SPECIES
    Al 26.98 Al_ONCV_PBE_sr.upf	#Element, Mass, Pseudopotential

    LATTICE_CONSTANT
    12.14569  			#Lattice constant

    LATTICE_VECTORS
    1.0 0.0 0.0 		#Lattice vector 1
    0.0 1.0 0.0		    #Lattice vector 2
    0.0 0.0 1.0 		#Lattice vector 3

    ATOMIC_POSITIONS
    Cartesian 		#Cartesian(Unit is LATTICE_CONSTANT)
    Al 			#Name of element	
    0.0			#Magnetic for this element.
    16			#Number of atoms
    0.0598184   0.789049    0.599279    1   1   1
    0.0379145   0.996286    0.225113    1   1   1
    0.970468    0.216979    0.785658    1   1   1
    0.625149    0.0674556   0.582236    1   1   1
    0.744445    0.531421    0.825099    1   1   1
    0.913838    0.377547    0.417868    1   1   1
    0.733554    0.756325    0.367093    1   1   1
    0.348249    0.284465    0.805284    1   1   1
    0.442703    0.538407    0.0932358   1   1   1
    0.575306    0.363412    0.338145    1   1   1
    0.241851    0.860453    0.996106    1   1   1
    0.717794    0.905869    0.887799    1   1   1
    0.450016    0.564658    0.611414    1   1   1
    0.401454    0.9478      0.374336    1   1   1
    0.662405    0.250878    0.0369698   1   1   1
    0.0750363   0.548996    0.0543649   1   1   1
    ```
    We provide an [instruction](../input-stru.md) on the speicifications of the STRU file.

- KPT
    ```
    K_POINTS
    0
    Gamma
    1 1 1 0 0 0
    ```
    More information on the KPT file can be found in this [short instruction](../input-kpt).


[back to top](#stochastic-dft-and-mix-stochastic-deterministic-dft)
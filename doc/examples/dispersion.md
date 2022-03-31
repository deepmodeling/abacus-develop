# DFT + dispersion calculations

[back to main page](../../README.md)

LDA and GGAs cannot describe van der Waals (vdW) interactions in a physically correct way. In order to describe materials where vdW interactions are important, one needs to go beyond LDA and GGAs. To this end, one simple and popular approach is to add a Lennard-Jones type term in terms of atom-based pairwise C6=R6 summations to the existing GGA functionals. There are different ways to do this. Currently ABACUS provides three Grimme DFT-D methods, including [D2](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20495), [D3(0)](https://aip.scitation.org/doi/10.1063/1.3382344) and [D3(BJ)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21759), to describe van der Waals interactions. Among them, the D3 method has been implemented in ABACUS based on the
dftd3 [program](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3) written by Stefan Grimme, Stephan Ehrlich and Helge Krieg.

The user can set the vdW-corrected methods in INPUT file with the parameter, vdw_method:

   - none: the vdW corrections are not included
   - d2: DFT-D2 method
   - d3_0: DFT-D3(0) method
   - d3_bj: DFT-D3(BJ) method

For example, if you want to run a DFT-D2 calculation, you should add the following keyword to the INPUT file.
```
vdw_method d2
```
By this way, ABACUS will calculate total energies and forces with DFT-D2 correction. For more parameters please refer to the [list of input variables](../input-main.md#vdw-correction).

Below is an example using DFT + dispersion calculation:
- INPUT:

```
INPUT_PARAMETERS
ntype                         2
ecutwfc                       20
scf_thr_rho                           1e-06
scf_nmax                         400
basis_type                    lcao
ks_solver                     genelpa
smearing_method                      gaussian
smearing_sigma                         0.02
mixing_type                   pulay
mixing_beta                   0.4
vdw_method                    d2
calculation                   scf
cal_force                         1
cal_stress                        1
```

- STRU:

```
ATOMIC_SPECIES
Sn  118.710  Sn_ONCV_PBE-1.0.upf
Te  127.603  Te_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Sn_pbe_9.0au_100Ry_2s2p2d
Te_pbe_9.0au_100Ry_2s2p2d

LATTICE_CONSTANT
1.89035917

LATTICE_VECTORS
4.646 0 0 #latvec1
0 4.561 0 #latvec2
0 0 16.32 #latvec3

ATOMIC_POSITIONS
Direct

Sn #label
0 #magnetism
1 #number of atoms
0.5493  0.5000000000000000  0.2 1 1 1

Te #label
0 #magnetism
1 #number of atoms
0.0000000000000000  0.0000000000000000  0.2 1 1 1
```

- KPT:

```
K_POINTS
0
Monkhorst-Pack
2 2 2 0 0 0
```

- Pseudopotential and atomic basis
The files (`Sn_ONCV_PBE-1.0.upf`,`Te_ONCV_PBE-1.0.upf`,`Sn_pbe_9.0au_100Ry_2s2p2d`,`Te_pbe_9.0au_100Ry_2s2p2d`) can be found in the directory tests/tools/PP_ORB/

- Output
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
 Fri Jul  2 11:46:13 2021
 MAKE THE DIR         : OUT.ABACUS/
 DONE(0.127021   SEC) : SETUP UNITCELL
 DONE(0.128917   SEC) : INIT K-POINTS
 ---------------------------------------------------------
 Self-consistent calculations for electrons
 ---------------------------------------------------------
 SPIN    KPOINTS         PROCESSORS  NBASE       
 1       8               1           36          
 ---------------------------------------------------------
 Use Systematically Improvable Atomic bases
 ---------------------------------------------------------
 ELEMENT ORBITALS        NBASE       NATOM       XC          
 Sn      2s2p2d-9au      18          1           PBE
 Te      2s2p2d-9au      18          1           PBE
 ---------------------------------------------------------
 Initial plane wave basis and FFT box
 ---------------------------------------------------------
 SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS
 SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS
 DONE(4.49436    SEC) : INIT PLANEWAVE
 UNIFORM GRID DIM     : 30 * 30 * 90
 UNIFORM GRID DIM(BIG): 15 * 15 * 45
 DONE(4.49762    SEC) : INIT CHARGE
 DONE(4.4993     SEC) : INIT POTENTIAL
 START POTENTIAL      : atomic
 ---------------------------------------------------------
 SELF-CONSISTENT : 
 ---------------------------------------------------------
 ITER   ETOT(eV)       EDIFF(eV)      SCF_THR_RHO      TIME(s)    
 GE1    -4.263615e+03  0.000000e+00   1.030e-01  1.030e+00  
 GE2    -4.262783e+03  8.325232e-01   6.041e-02  1.010e+00  
 GE3    -4.262610e+03  1.722243e-01   2.263e-02  1.040e+00  
 GE4    -4.262486e+03  1.247259e-01   4.662e-03  1.000e+00  
 GE5    -4.262576e+03  -9.092446e-02  5.323e-04  9.700e-01  
 GE6    -4.262571e+03  5.648702e-03   1.910e-04  9.700e-01  
 GE7    -4.262569e+03  1.741079e-03   3.454e-05  9.600e-01  
 GE8    -4.262569e+03  1.839695e-04   7.188e-06  9.500e-01  
 GE9    -4.262569e+03  -1.606959e-05  1.728e-06  1.010e+00  
 GE10   -4.262569e+03  -2.716369e-05  4.516e-07  9.300e-01  
 ><><><><><><><><><><><><><><><><><><><><><><
 TOTAL-STRESS (KBAR):
 ><><><><><><><><><><><><><><><><><><><><><><
 -2.121e+02     -8.470e-07     -2.354e-10     
 -8.470e-07     1.142e+02      -3.656e-13     
 -5.953e-10     -2.662e-13     -8.366e+01     

  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
 A Run_lcao            lcao_line           15.77          1         16        1e+02     %
 B ORB_control         set_orb_tables      4.3448         1         4.3       28        %
 C LCAO_Orbitals       Read_Orbitals       2.327          1         2.3       15        %
 C ORB_gen_tables      gen_tables          2.0178         1         2         13        %
 D ORB_table_phi       init_Table          1.119          1         1.1       7.1       %
 D ORB_table_beta      init_Table_Beta     0.79881        1         0.8       5.1       %
 X FFT                 FFT3D               0.17197        175       0.00098   1.1       %
 E Potential           v_of_rho            0.50032        12        0.042     3.2       %
 B LOOP_ions           opt_ions            11.079         1         11        70        %
 C LOOP_elec           solve_elec_stru     4.5258         1         4.5       29        %
 D LOOP_elec           before_solver       0.53419        1         0.53      3.4       %
 E LCAO_Hamilt         set_lcao_matrices   0.53412        1         0.53      3.4       %
 G LCAO_gen_fixedH     build_Nonlocal_mu   1.3624         2         0.68      8.6       %
 X ORB_gen_tables      snap_psibeta        2.6413         209952    1.3e-05   17        %
 D LOOP_elec           solver              3.9777         1         4         25        %
 D ELEC_scf            scf                 3.9775         1         4         25        %
 E ELEC_cbands_k       cal_bands           1.8355         10        0.18      12        %
 K Gint_k              vlocal              1.635          10        0.16      10        %
 G Diago_LCAO_Matrix   genelpa             0.10715        40        0.0027    0.68      %
 G Diago_LCAO_Matrix   genelpa2            0.10148        80        0.0013    0.64      %
 E Local_Orbital_Cha   sum_bands           1.6351         10        0.16      10        %
 F Gint_k              charge              1.6077         10        0.16      10        %
 E Charge              mix_rho             0.10381        10        0.01      0.66      %
 F Force_LCAO_k        ftable_k            6.3736         1         6.4       40        %
 G Force_LCAO_k        cal_fvl_dphi_k      3.8463         1         3.8       24        %
 G Force_LCAO_k        cal_fvnl_dbeta_k    1.5456         1         1.5       9.8       %
 ----------------------------------------------------------------------------------------

 START  Time  : Fri Jul  2 11:46:13 2021
 FINISH Time  : Fri Jul  2 11:46:29 2021
 TOTAL  Time  : 16
 SEE INFORMATION IN : OUT.ABACUS/

```

[back to top](#DFT-+-dispersion-calculations)
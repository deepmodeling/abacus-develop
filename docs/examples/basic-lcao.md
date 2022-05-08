# Basic electronic structure calculation with LCAO basis set

[back to main page](../../README.md)

In this section we will describe how to do LCAO calculations using ABACUS. Again the crystal Si in the diamond structure will be taken as an example.

In the work directory, copy the `INPUT`, `KPT` and `STRU` file, the pseudopotential file(`Si.pz-vbc.UPF`), and in addition the numerical atomic orbital file(`Si_lda_8.0au_50Ry_2s2p1d`) from the `examples/02a_Si_diamond_lcao_scf/` directory.

Note the keyword
```
basis_type              lcao
```
in the INPUT file.

The information printed on the screen is different from that obtained using the plane-wave
basis:
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
 Fri Jul  2 11:23:16 2021
 MAKE THE DIR         : OUT.ABACUS/
 DONE(0.0615518  SEC) : SETUP UNITCELL
 DONE(0.0656908  SEC) : INIT K-POINTS
 ---------------------------------------------------------
 Self-consistent calculations for electrons
 ---------------------------------------------------------
 SPIN    KPOINTS         PROCESSORS  NBASE       
 1       64              1           26          
 ---------------------------------------------------------
 Use Systematically Improvable Atomic bases
 ---------------------------------------------------------
 ELEMENT ORBITALS        NBASE       NATOM       XC          
 Si      2s2p1d-8au      13          2           PZ-LDA
 ---------------------------------------------------------
 Initial plane wave basis and FFT box
 ---------------------------------------------------------
 SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS
 DONE(1.73298    SEC) : INIT PLANEWAVE
 UNIFORM GRID DIM     : 36 * 36 * 36
 UNIFORM GRID DIM(BIG): 18 * 18 * 18
 DONE(1.7346     SEC) : INIT CHARGE
 DONE(1.73558    SEC) : INIT POTENTIAL
 START POTENTIAL      : atomic
 ---------------------------------------------------------
 SELF-CONSISTENT : 
 ---------------------------------------------------------
 ITER   ETOT(eV)       EDIFF(eV)      SCF_THR      TIME(s)    
 GE1    -2.138667e+02  0.000000e+00   1.652e-01  4.690e+00  
 GE2    -2.139153e+02  -4.859216e-02  3.480e-02  4.970e+00  
 GE3    -2.139161e+02  -8.407097e-04  4.131e-03  4.760e+00  
 GE4    -2.139161e+02  -1.051285e-05  4.859e-05  4.630e+00  
 GE5    -2.139161e+02  2.951636e-08   6.059e-06  5.120e+00  
 GE6    -2.139161e+02  -3.891962e-09  1.608e-07  5.210e+00  
 GE7    -2.139161e+02  1.061725e-10   7.775e-10  4.930e+00  

  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
 A Run_lcao            lcao_line           25.505         1         26        1e+02     %
 B ORB_control         set_orb_tables      1.6457         1         1.6       6.4       %
 C LCAO_Orbitals       Read_Orbitals       1.0609         1         1.1       4.2       %
 C ORB_gen_tables      gen_tables          0.58477        1         0.58      2.3       %
 D ORB_table_phi       init_Table          0.33613        1         0.34      1.3       %
 B LOOP_ions           opt_ions            23.744         1         24        93        %
 C LOOP_elec           solve_elec_stru     23.744         1         24        93        %
 D LOOP_elec           set_matrix_grid     0.13358        1         0.13      0.52      %
 D Grid_Technique      init                0.10495        1         0.1       0.41      %
 D LOOP_elec           before_solver       5.2475         1         5.2       21        %
 E LCAO_Hamilt         set_lcao_matrices   5.2453         1         5.2       21        %
 G LCAO_gen_fixedH     build_Nonlocal_mu   5.1953         1         5.2       20        %
 X ORB_gen_tables      snap_psibeta        3.6577         1703858   2.1e-06   14        %
 D LOOP_elec           solver              18.363         1         18        72        %
 D ELEC_scf            scf                 18.363         1         18        72        %
 E ELEC_cbands_k       cal_bands           13.126         7         1.9       51        %
 K Gint_k              vlocal              8.6586         7         1.2       34        %
 F LCAO_Hamilt         calculate_Hk        3.927          448       0.0088    15        %
 G LCAO_nnr            folding_fixedH      3.7241         448       0.0083    15        %
 G Diago_LCAO_Matrix   genelpa             2.2571         224       0.01      8.8       %
 G Diago_LCAO_Matrix   genelpa1            0.1584         448       0.00035   0.62      %
 G Diago_LCAO_Matrix   genelpa2            0.29915        448       0.00067   1.2       %
 E Local_Orbital_Cha   sum_bands           5.1631         7         0.74      20        %
 F LCAO_Charge         cal_dk_k            0.20373        7         0.029     0.8       %
 F Gint_k              charge              4.9517         7         0.71      19        %
 ----------------------------------------------------------------------------------------

 START  Time  : Fri Jul  2 11:23:16 2021
 FINISH Time  : Fri Jul  2 11:23:42 2021
 TOTAL  Time  : 26
 SEE INFORMATION IN : OUT.ABACUS/

```

The string GEn in the first column means that the genelpa eigenvalue solver is used, and this is the n-th self-consistent KS iteration. In contrast, the output information from the PW calculation has the string CGn in its first column, indicating the conjugate gradients (CG) method is used to solve the Kohn-Sham equation.

In many cases, only one &Gamma; point (i.e., k=0) is needed in the calculations. In these cases, one can set in the INPUT file:

- gamma_only If set to 1, only &Gamma; point is used in the calculation. The gamma_only algorithm is expected to be faster than the standard algorithm.

[back to top](#basic-electronic-structure-calculation-with-lcao-basis-set)
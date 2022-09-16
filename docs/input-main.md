# INPUT file

## Table of contents

### [Structure of the file](#structure-of-the-file)

### [List of keywords](#list-of-keywords)

- [System variables](#system-variables)

    [suffix](#suffix) | [ntype](#ntype) | [calculation](#calculation) | [symmetry](#symmetry) | [kpar](#kpar) | [bndpar](#bndpar) | [latname](#latname) | [init_wfc](#init_wfc) | [init_chg](#init_chg) | [init_vel](#init_vel) | [nelec](#nelec) | [tot_magnetization](#tot_magnetization) | [dft_functional](#dft_functional) | [xc_temperature](#xc_temperature) | [pseudo_type](#pseudo_type) |  [pseudo_rcut](#pseudo_rcut) | [pseudo_mesh](#pseudo_mesh) | [mem_saver](#mem_saver) | [diago_proc](#diago_proc) | [nbspline](#nbspline) | [kspacing](#kspacing)  | [min_dist_coef](#min_dist_coef)

- [Variables related to input files](#variables-related-to-input-files)

    [stru_file](#stru_file) | [kpoint_file](#kpoint_file) | [pseudo_dir](#pseudo_dir) | [orbital_dir](#orbital_dir) | [read_file_dir](#read_file_dir)

- [Plane wave related variables](#plane-wave-related-variables)

    [ecutwfc](#ecutwfc) | [nx,ny,nz](#nx-ny-nz) | [pw_seed](#pw_seed) | [pw_diag_thr](#pw_diag_thr) | [pw_diag_nmax](#pw_diag_nmax) | [pw_diag_ndim](#pw_diag_ndim)

- [Numerical atomic orbitals related variables](#numerical-atomic-orbitals-related-variables)

    [nb2d](#nb2d) | [lmaxmax](#lmaxmax) | [lcao_ecut](#lcao_ecut) | [lcao_dk](#lcao_dk) | [lcao_dr](#lcao_dr) | [lcao_rmax](#lcao_rmax) | [search_radius](#search_radius) | [search_pbc](#search_pbc)

- [Electronic structure](#electronic-structure)

    [basis_type](#basis_type) | [ks_solver](#ks_solver) | [nbands](#nbands) | [nbands_istate](#nbands_istate) | [nspin](#nspin) | [occupations](#occupations) | [smearing_method](#smearing_method) | [smearing_sigma](#smearing_sigma) | [mixing_type](#mixing_type) | [mixing_beta](#mixing_beta) | [mixing_ndim](#mixing_ndim) | [mixing_gg0](#mixing_gg0) | [gamma_only](#gamma_only) | [printe](#printe) | [scf_nmax](#scf_nmax) | [scf_thr](#scf_thr) | [chg_extrap](#chg_extrap)

- [Electronic structure (SDFT)](#electronic-structure-sdft)

    [method_sto](#method_sto) | [nbands_sto](#nbands_sto) | [nche_sto](#nche_sto) | [emin_sto](#emin_sto) | [emax_sto](#emax_sto) | [seed_sto](#seed_sto) | [initsto_freq](#initsto_freq)

- [Geometry relaxation](#geometry-relaxation)

    [relax_nmax](#relax_nmax) | [relax_method](#relax_method) | [relax_cg_thr](#relax_cg_thr) | [relax_bfgs_w1](#relax_bfgs_w1) | [relax_bfgs_w2](#relax_bfgs_w2) | [relax_bfgs_rmax](#relax_bfgs_rmax) | [relax_bfgs_rmin](#relax_bfgs_rmin) | [relax_bfgs_init](#relax_bfgs_init) | [cal_force](#cal_force) | [force_thr](#force_thr) | [force_thr_ev](#force_thr_ev) | [cal_stress](#cal_stress) | [stress_thr](#stress_thr) | [press1, press2, press3](#press1-press2-press3) | [fixed_axes](#fixed_axes) | [cell_factor](#cell_factor)

- [Variables related to output information](#variables-related-to-output-information)

    [out_force](#out_force) | [out_mul](#out_mul) | [out_freq_elec](#out_freq_elec) | [out_freq_ion](#out_freq_ion) | [out_chg](#out_chg) | [out_pot](#out_pot) | [out_dm](#out_dm) | [out_wfc_pw](#out_wfc_pw) | [out_wfc_r](#out_wfc_r) | [out_wfc_lcao](#out_wfc_lcao) | [out_dos](#out_dos) | [out_band](#out_band) | [out_proj_band](#out_proj_band) | [out_stru](#out_stru) | [out_level](#out_level) | [out_alllog](#out_alllog) | [out_mat_hs](#out_mat_hs) | [out_mat_r](#out_mat_r) | [out_mat_hs2](#out_mat_hs2) | [out_element_info](#out_element_info) | [restart_save](#restart_save) | [restart_load](#restart_load)

- [Density of states](#density-of-states)

    [dos_edelta_ev](#dos_edelta_ev) | [dos_sigma](#dos_sigma) | [dos_scale](#dos_scale) | [dos_nche](#dos_nche)

- [Exact exchange](#exact-exchange) (Under tests)

    [exx_hybrid_type](#exx_hybrid_type) | [exx_hybrid_alpha](#exx_hybrid_alpha) | [exx_hse_omega](#exx_hse_omega) | [exx_separate_loop](#exx_separate_loop) | [exx_hybrid_step](#exx_hybrid_step) | [exx_lambda](#exx_lambda) | [exx_pca_threshold](#exx_pca_threshold) | [exx_c_threshold](#exx_c_threshold) | [exx_v_threshold](#exx_v_threshold) | [exx_dm_threshold](#exx_dm_threshold) | [exx_schwarz_threshold](#exx_schwarz_threshold) | [exx_cauchy_threshold](#exx_cauchy_threshold) | [exx_ccp_threshold](#exx_ccp_threshold) | [exx_ccp_rmesh_times](#exx_ccp_rmesh_times) | [exx_distribute_type](#exx_distribute_type) | [exx_opt_orb_lmax](#exx_opt_orb_lmax) | [exx_opt_orb_ecut](#exx_opt_orb_ecut) | [exx_opt_orb_tolerence](#exx_opt_orb_tolerence)

- [Molecular dynamics](#molecular-dynamics)

    [md_type](#md_type) | [md_nstep](#md_nstep) | [md_ensolver](#md_ensolver) | [md_restart](#md_restart) | [md_dt](#md_dt) | [md_t](#md_t) | [md_dumpfreq](#md_dumpfreq) | [md_restartfreq](#md_restartfreq) | [md_seed](#md_seed) | [md_tfreq](#md_tfreq) | [md_mnhc](#md_mnhc) | [lj_rcut](#lj_rcut) | [lj_epsilon](#lj_epsilon) | [lj_sigma](#lj_sigma) | [msst_direction](#msst_direction) | [msst_vel](#msst_vel) | [msst_vis](#msst_vis) | [msst_tscale](#msst_tscale) | [msst_qmass](#msst_qmass) | [md_damp](#md_damp)

- [vdW correction](#vdw-correction)

    [vdw_method](#vdw_method) | [vdw_s6](#vdw_s6) | [vdw_s8](#vdw_s8) | [vdw_a1](#vdw_a1) | [vdw_a2](#vdw_a2) | [vdw_d](#vdw_d) | [vdw_abc](#vdw_abc) | [vdw_C6_file](#vdw_C6_file) | [vdw_C6_unit](#vdw_C6_unit) | [vdw_R0_file](#vdw_R0_file) | [vdw_R0_unit](#vdw_R0_unit) | [vdw_model](#vdw_model) | [vdw_radius](#vdw_radius) | [vdw_radius_unit](#vdw_radius_unit) | [vdw_cn_radius](#vdw_cn_radius) | [vdw_cn_radius_unit](#vdw_cn_radius_unit) | [vdw_period](#vdw_period)

- [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)

    [berry_phase](#berry_phase) | [gdir](#gdir) | [towannier90](#towannier90) | [nnkpfile](#nnkpfile) | [wannier_spin](#wannier_spin)

- [TDDFT: time dependent density functional theory](#tddft-time-dependent-density-functional-theory) (Under tests)

    [tddft](#tddft) | [td_scf_thr](#td_scf_thr) | [td_dt](#td_dt) | [td_force_dt](#td_force_dt) | [td_vext](#td_vext) | [td_vext_dire](#td_vext_dire) | [td_timescale](#td_timescale) | [td_vexttype](#td_vexttype) | [td_vextout](#td_vextout) | [td_dipoleout](#td_dipoleout) | [ocp](#ocp) | [ocp_set](#ocp_set)

- [DFT+U correction](#DFTU-correction) (Under tests)

    [dft_plus_u](#dft_plus_u) | [orbital_corr](#orbital_corr) | [hubbard_u](#hubbard_u) | [hund_j](#hund_j) | [yukawa_potential](#yukawa_potential) | [omc](#omc)

- [Variables useful for debugging](#variables-useful-for-debugging)

    [nurse](#nurse) | [t_in_h](#t_in_h) | [vl_in_h](#vl_in_h) | [vnl_in_h](#vnl_in_h) | [test_force](#test_force) | [test_stress](#test_stress) | [colour](#colour)

- [DeePKS](#deepks)

    [deepks_out_labels](#deepks_out_labels) | [deepks_descriptor_lmax](#deepks_descriptor_lmax) | [deepks_scf](#deepks_scf) | [deepks_model](#deepks_model)

- [Electric field and dipole correction](#electric-field-and-dipole-correction)

    [efield_flag](#efield_flag) | [dip_cor_flag](#dip_cor_flag) | [efield_dir](#efield_dir) | [efield_pos_max](#efield_pos_max) | [efield_pos_dec](#efield_pos_dec) | [efield_amp ](#efield_amp )

- [Electronic conductivities](#electronic-conductivities)

    [cal_cond](#cal_cond) | [cond_nche](#cond_nche) | [cond_dw](#cond_dw) | [cond_wcut](#cond_wcut) | [cond_wenlarge](#cond_wenlarge) | [cond_fwhm ](#cond_fwhm )

- [Implicit solvation model](#implicit-solvation-model)

    [imp_sol](#imp_sol) | [eb_k](#eb_k) | [tau](#tau) | [sigma_k](#sigma_k) | [nc_k](#nc_k) 

[back to main page](../README.md)

## Structure of the file

Below is an example INPUT file with some of the most important parameters that need to be set:

```
INPUT_PARAMETERS
#Parameters (General)
ntype 1
nbands 4
#Parameters (Accuracy)
ecutwfc 60
```

Parameters list starts with key word `INPUT_PARAMETERS`. Any content before `INPUT_PARAMETERS` will be ignored.

Any line starting with `#` or `/` will also be ignored.

Each parameter value is provided by specifying the name of the input variable
and then putting the value after the name, separated by one or more blank characters(space
or tab). The following characters(≤ 150) in the same line will be neglected.

Depending on the input variable, the value may be an integer, a real number or a string. The parameters can be given in any order, but only one parameter should be given per line.

Furthermore, if a given parameter name appeared more than once in the input file, only the last value will be taken.

> Note: if a parameter name is not recognized by the program, the program will stop with an error message.

In the above example, the meanings of the parameters are:

- `ntype` : how many types of elements in the unit cell
- `nbands` : the number of bands to be calculated
- `ecutwfc` : the plane-wave energy cutoff for the wave function expansion (UNIT: Rydberg)

## List of keywords

### System variables

This part of variables are used to control general system parameters.

#### suffix

- **Type**: String
- **Description**: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the `suffix` is the name you can pick up for your convenience.
- **Default**: ABACUS

#### ntype

- **Type**: Integer
- **Description**: Number of different atom species in this calculations. This value must be set. If the number you set is smaller than the atom species in the STRU file, ABACUS only read the `wrong number` of atom information. If the number is larger than the atom species in the STRU file, ABACUS may stop and quit.
- **Default**: ***No default value***

#### calculation

- **Type**: String
- **Description**: Specify the type of calculation.
  - *scf*: do self-consistent electronic structure calculation
  - *relax*: do structure relaxation calculation, one can ues `relax_nmax` to decide how many ionic relaxations you want.
  - *cell-relax*: do cell relaxation calculation.
  - *nscf*: do the non self-consistent electronic structure calculations. For this option, you need a charge density file. For nscf calculations with planewave basis set, pw_diag_thr should be <= 1d-3.
  - *istate*: Please see the explanation for variable `nbands_istate`.
  - *ienvelope*: Please see the explanation for variable `nbands_istate`.
  - *md*: molecular dynamics
  - *sto-scf*: do self-consistent electronic structure calculation with [stochastic DFT](#electronic-structure-sdft)
  - *sto-md*: molecular dynamics with [stochastic DFT](#electronic-structure-sdft)
  - *test_memory* : checks memory required for the calculation. The number is not quite reliable, please use with care
  - *test_neighbour* : only performs neighbouring atom search
  - *gen_jle* : generates projectors for DeePKS; see also deepks_lmax_descriptor

  > Note: *istate* and *ienvelope* only work for LCAO basis set and are not working right now.
- **Default**: scf

#### symmetry

- **Type**: Integer
- **Description**: takes value 1, 0 and -1. 
  - if set to 1, symmetry analysis will be performed to determine the type of Bravais lattice and associated symmetry operations. (point groups only)
  - if set to 0, only time reversal symmetry would be considered in symmetry operations, which implied k point and -k point would be treated as one double weight k point.
  - if set to -1, any symmetry will not be considered.
- **Default**: 0

#### kpar

- **Type**: Integer
- **Description**: devide all processors into kpar groups, and k points will be distributed among each group. The value taken should be less than or equal to the number of k points as well as the number of MPI threads.
- **Default**: 1

#### bndpar

- **Type**: Integer
- **Description**: devide all processors into bndpar groups, and bands (only stochastic orbitals now) will be distributed among each group. It should be larger than 0.
- **Default**: 1

#### latname

- **Type**: String
- **Description**: Specifies the type of Bravias lattice. When set to `test`, the three lattice vectors are supplied explicitly in STRU file. When set to certain Bravais lattice type, there is no need to provide lattice vector, but a few lattice parameters might be required. For more information regarding this parameter, consult the [page on STRU file](input-stru.md).
    Available options are:
  - `test`: free strcture.
  - `sc`: simple cubie.
  - `fcc`: face-centered cubic.
  - `bcc`: body-centered cubic.
  - `hexagonal`: hexagonal.
  - `trigonal`: trigonal.
  - `st`: simple tetragonal.
  - `bct`: body-centered tetragonal.
  - `so`: orthorhombic.
  - `baco`: base-centered orthorhombic.
  - `fco`: face-centered orthorhombic.
  - `bco`: body-centered orthorhombic.
  - `sm`: simple monoclinic.
  - `bacm`: base-centered monoclinic.
  - `triclinic`: triclinic.
- **Default**: `test`

#### init_wfc

- **Type**: String
- **Description**: Only useful for plane wave basis only now. It is the name of the starting wave functions. In the future we should also make this         variable available for localized orbitals set.
    Available options are:
  - `atomic`: from atomic pseudo wave functions. If they are not enough, other wave functions are initialized with random numbers.
  - `atomic+random`: add small random numbers on atomic pseudo-wavefunctions
  - `file`: from file
  - `random`: random numbers
- **Default**:`atomic`

#### init_chg

- **Type**: String
- **Description**: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density. If set this to `atomic`, the density is starting from summation of atomic density of single atoms. If set this to `file`, the density will be read in from file. Besides, when you do `nspin=1` calculation, you only need the density file SPIN1_CHGCAR. However, if you do `nspin=2` calculation, you also need the density file SPIN2_CHGCAR. The density file should be output with these names if you set out_chg = 1 in INPUT file.
- **Default**: atomic

#### init_vel

- **Type**: Boolean
- **Description**: Read the atom velocity from the atom file (STRU) if set to true.
- **Default**: false

#### nelec

- **Type**: Real
- **Description**: If >0.0, this denotes total number of electrons in the system. Must be less than 2*nbands. If set to 0.0, the total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
- **Default**: 0.0

#### tot_magnetization

- **Type**: Real
- **Description**: Total magnetization of the system.
- **Default**: 0.0

#### dft_functional

- **Type**: String
- **Description**: type of exchange-correlation functional used in calculation. If dft_functional is not set, the program will adopt the functional used to generate pseudopotential files, provided all of them are generated using the same functional. For example, we present a few lines in Si’s GGA pseudopotential file Si_ONCV_PBE-1.0.upf:
        ```
        ...
        <PP_HEADER
        generated="Generated using ONCVPSP code by D. R. Hamann"
        author="Martin Schlipf and Francois Gygi"
        date="150105"
        comment=""
        element="Si"
        pseudo_type="NC"
        relativistic="scalar"
        is_ultrasoft="F"
        is_paw="F"
        is_coulomb="F"
        has_so="F"
        has_wfc="F"
        has_gipaw="F"
        core_correction="F"
        functional="PBE"
        z_valence=" 4.00"
        total_psenergy=" -3.74274958433E+00"
        rho_cutoff=" 6.01000000000E+00"
        ```
    According to the information above, this pseudopotential is generated using PBE functional.
    On the other hand, if dft_functional is specified, it will overwrite the functional from pseudopotentials and performs calculation with whichever functional the user prefers. We further offer two ways of supplying exchange-correlation functional. The first is using 'short-hand' names such as 'LDA', 'PBE', 'SCAN'. A complete list of 'short-hand' expressions can be found in [source code](../source/module_xc/xc_functional.cpp). The other way is only available when ***compiling with LIBXC***, and it allows for supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by plus sign, for example, 'dft_functional='LDA_X_1D_EXPONENTIAL+LDA_C_1D_CSC'. The list of LIBXC keywords can be found on its [website](https://www.tddft.org/programs/libxc/functionals/). In this way, **we support all the LDA,GGA and mGGA functionals provided by LIBXC**.

    Furthermore, the old INPUT parameter exx_hybrid_type for hybrid functionals has been absorbed into dft_functional. Options are `hf` (pure Hartree-Fock), `pbe0`(PBE0), `hse` (Note: in order to use HSE functional, LIBXC is required). Note also that HSE has been tested while PBE0 has NOT been fully tested yet, and the maximum parallel cpus for running exx is Nx(N+1)/2, with N being the number of atoms. And forces for hybrid functionals are not supported yet.

    If set to `opt_orb`, the program will not perform hybrid functional calculation. Instead, it is going to generate opt-ABFs as discussed in this [article](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c00481).
    
- **Default**: same as UPF file.

#### xc_temperature
- **Type**: Double
- **Description**: specifies temperature when using temperature-dependent XC functionals (KSDT and so on); unit in Rydberg
- **Default** : 0.0

#### pseudo_type

- **Type**: String
- **Description**: the format of pseudopotential files. Accepted value s are:
  - upf : .UPF format
  - vwr : .vwr format
  - upf201 : the new UPF format
  - blps : bulk local pseudopotential
- **Default** : upf

#### pseudo_rcut

- **Type**: Real
- **Description**: Cut-off of radial integration for pseudopotentials, in Bohr.
- **Default**: 15

#### pseudo_mesh

- **Type**: Integer
- **Description**: If set to 0, then use our own mesh for radial integration of pseudopotentials; if set to 1, then use the mesh that is consistent with quantum espresso.
- **Default**: 0

#### mem_saver

- **Type**: Boolean
- **Description**: Used only for nscf calculations. If set to 1, then a memory saving technique will be used for many k point calculations.
- **Default**: 0

#### diago_proc

- **Type**: Integer
- **Descrption**: If set to a positive number, then it specifies the number of threads used for carrying out diagonalization. Must be less than or equal to total number of MPI threads. Also, when cg diagonalization is used, diago_proc must be same as total number of MPI threads. If set to 0, then it will be set to the number of MPI threads. Normally, it is fine just leaving it to default value. Only used for pw base.
- **Default**: 0

#### nbspline

- **Type**: Integer
- **Descrption**: If set to a natural number, a Cardinal B-spline interpolation will be used to calculate Structure Factor. `nbspline` represents the order of B-spline basis and larger one can get more accurate results but cost more.
    It is turned off by default.
- **Default**: -1

#### kspacing

- **Type**: double
- **Descrption**: Set the smallest allowed spacing between k points, unit in 1/bohr. It should be larger than 0.0, and suggest smaller than 0.25. When you have set this value > 0.0, then the KPT file is unneccessary, and the number of K points nk_i = max(1,int(|b_i|/KSPACING)+1), where b_i is the reciprocal lattice vector. The default value 0.0 means that ABACUS will read the applied KPT file. Notice: if gamma_only is set to be true, kspacing is invalid.
- **Default**: 0.0

#### min_dist_coef

- **Type**: double
- **Descrption**: a factor related to the allowed minimum distance between two atoms. At the begining, ABACUS will check the structure, and if the distance of two atoms is shorter than min_dist_coef*(standard covalent bond length), we think this structure is unreasonable. If you want to calculate some structures in extreme condition like high pressure, you should set this parameter as a smaller value or even 0.
- **Default**: 0.2

### Variables related to input files

This part of variables are used to control input files related parameters.

#### stru_file

- **Type**: String
- **Description**: This parameter specifies the name of structure file which contains various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
- **Default**: STRU

#### kpoint_file

- **Type**: String
- **Description**: This parameter specifies the name of k-points file. Note that if you use atomic orbitals as basis, and you only use gamma point, you don't need to have k-point file in your directory, ABACUS will automatically generate `KPT` file. Otherwise, if you use more than one k-point, please do remember the algorithm in ABACUS is different for gamma only and various k-point dependent simulations. So first you should turn off the k-point algorithm by set `gamma_only = 0` in `INPUT` and then you should setup your own k-points file.
- **Default**: KPT

#### pseudo_dir

- **Type**: String
- **Description**: This parameter specifies pseudopotential directory.
- **Default**: ./

#### orbital_dir

- **Type**: String
- **Description**: This parameter specifies orbital file directory.
- **Default**: ./

#### read_file_dir

- **Type**: String
- **Description**: when the program needs to read files such as electron density(`SPIN1_CHG`) as a starting point, this variables tells the location of the files. For example, './' means the file is located in the working directory.
- **Default**: OUT.$suffix

### Plane wave related variables

This part of variables are used to control the plane wave related parameters.

#### ecutwfc

- **Type**: Real
- **Description**: Energy cutoff for plane wave functions, the unit is **Rydberg**. Note that even for localized orbitals basis, you still need to setup a energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set, etc. Also, because our orbitals are generated by matching localized orbitals to a chosen set of wave functions from certain energy cutoff, so this set of localize orbitals are most accurate under this same plane wave energy cutoff.
- **Default**: 50

#### nx, ny, nz

- **Type**: Integer
- **Description**: If set to a positive number, then the three variables specify the numbers of FFT grid points in x, y, z directions, respectively. If set to 0, the number will be calculated from ecutwfc.
- **Default**: 0

#### pw_seed

- **Type**: Integer
- **Description**: Only useful for plane wave basis only now. It is the random seed to initialize wave functions. Only positive integers are avilable.
- **Default**:0

#### pw_diag_thr

- **Type**: Real
- **Description**: Only used when you use `diago_type = cg` or `diago_type = david`. It indicates the threshold for the first electronic iteration, from the second iteration the pw_diag_thr will be updated automatically. **For nscf calculations with planewave basis set, pw_diag_thr should be <= 1d-3.**
- **Default**: 0.01

#### pw_diag_nmax

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = cg` or `ks_solver = dav`. It indicates the maximal iteration number for cg/david method.
- **Default**: 40

#### pw_diag_ndim

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = dav`. It indicates the maximal dimension for the Davidson method.
- **Default**: 10

### Numerical atomic orbitals related variables

This part of variables are used to control the numerical atomic orbitals related parameters.

#### nb2d

- **Type**: Integer
- **Description**: In LCAO calculations, we arrange the total number of processors in an 2D array, so that we can partition the wavefunction matrix (number of bands*total size of atomic orbital basis) and distribute them in this 2D array. When the system is large, we group processors into sizes of nb2d, so that multiple processors take care of one row block (a group of atomic orbitals) in the wavefunction matrix. If set to 0, nb2d will be automatically set in the program according to the size of atomic orbital basis:
  - if size <= 500 : nb2d = 1
  - if 500 < size <= 1000 : nb2d = 32
  - if size > 1000 : nb2d = 64;
- **Default**: 0

#### lmaxmax

- **Type**: Integer
- **Description**: If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.
- **Default**: 2.

#### lcao_ecut

- **Type**: Real
- **Description**: Energy cutoff when calculating LCAO two-center integrals. In Ry.
- **Default**: 50

#### lcao_dk

- **Type**: Real
- **Description**: Delta k for 1D integration in LCAO
- **Default**: 0.01

#### lcao_dr

- **Type**: Real
- **Description**: Delta r for 1D integration in LCAO
- **Default**: 0.01

#### lcao_rmax

- **Type**: Real
- **Description**: Max R for 1D two-center integration table
- **Default**: 30

#### search_radius

- **Type**: Real
- **Description**: Set the search radius for finding neighbouring atoms. If set to -1, then the radius will be set to maximum of projector and orbital cut-off.
- **Default**: -1

#### search_pbc

- **Type**: Boolean
- **Description**: In searching for neighbouring atoms, if set to 1, then periodic images will also be searched. If set to 0, then periodic images will not be searched.
- **Default**: 1

### Electronic structure

This part of variables are used to control the electronic structure and geometry relaxation
calculations.

#### basis_type

- **Type**: String
- **Description**: This is very important parameters to choose basis set in ABACUS.
  - *pw*: Using plane-wave basis set only.
  - *lcao_in_pw*: Expand the localized atomic set in plane-wave basis.
  - lcao: Using localized atomic orbital sets.
- **Default**: pw

#### ks_solver

- **Type**: String
- **Description**: It`s about choice of diagonalization methods for hamiltonian matrix expanded in a certain basis set.

  For plane-wave basis,
  - cg: cg method.
  - dav: the Davidson algorithm. (Currently not working with Intel MKL library).

  For atomic orbitals basis,
  - genelpa: This method should be used if you choose localized orbitals.
  - hpseps: old method, still used.
  - lapack: lapack can be used for localized orbitals, but is only used for single processor.
  - cusolver: this method needs building with the cusolver component for lcao and at least one gpu is available.

   If you set ks_solver=`hpseps` for basis_type=`pw`, the program will be stopped with an error message:

    ```txt
    hpseps can not be used with plane wave basis.
    ```

    Then the user has to correct the input file and restart the calculation.
- **Default**: `cg` (pw) or `genelpa` (lcao)

#### nbands

- **Type**: Integer
- **Description**: Number of Kohn-Sham orbitals to calculate. It is recommended you setup this value, especially when you use smearing techniques, more bands should be included.
- **Default**:
  - nspin=1: 1.2\*occupied_bands, occupied_bands+10)
  - nspin=2: max(1.2\*nelec, nelec+20)

#### nbands_istate

- **Type**: Integer
- **Description**: Only used when `calculation = ienvelope` or `calculation = istate`, this variable indicates how many bands around Fermi level you would like to calculate. `ienvelope` means to calculate the envelope functions of wave functions $\Psi_{i}=\Sigma_{\mu}C_{i\mu}\Phi_{\mu}$, where $\Psi_{i}$ is the ith wave function with the band index $i$ and $\Phi_{\mu}$ is the localized atomic orbital set. `istate` means to calculate the density of each wave function $|\Psi_{i}|^{2}$. Specifically, suppose we have highest occupied bands at 100th wave functions. And if you set this variable to 5, it will print five wave functions from 96th to 105th. But before all this can be carried out, the wave functions coefficients  should be first calculated and written into a file by setting the flag `out_wfc_lcao = 1`.
- **Default**: 5

#### nspin

- **Type**: Integer
- **Description**: Number of spin components of wave functions. There are only two choices now: 1 or 2, meaning non spin or collinear spin.
- **Default**: 1

#### occupations

- **Type**: String
- **Description**: Specifies how to calculate the occupations of bands. Available options are:
  - 'smearing' : gaussian smearing for metals; see also variables `smearing_method` and `smearing_sigma`.
  - 'tetrahedra' : Tetrahedron method, Bloechl's version: [P.E. Bloechl, PRB 49, 16223 (1994)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223). Requires a uniform grid of k-points that are automatically generated. Well suited for calculation of DOS, less so (because not variational) for force/optimization/dynamics calculations.
  - 'fixed' : for insulators with a gap
- **Default**: 'smearing'

#### smearing_method

- **Type**: String
- **Description**: It indicates which occupation and smearing method is used in the calculation.
  - fixed: use fixed occupations.
  - gauss or gaussian: use gaussian smearing method.
  - mp: use methfessel-paxton smearing method. The method recommends for metals.
- **Default**: fixed

#### smearing_sigma

- **Type**: Real
- **Description**: energy range for smearing, the unit is Rydberg.
- **Default**: 0.001

#### mixing_type

- **Type**: String
- **Description**: Charge mixing methods.
  - plain: Just simple mixing.
  - kerker: Use kerker method, which is the mixing method in G space.
  - pulay: Standard Pulay method.
  - pulay-kerker:
- **Default**: pulay

#### mixing_beta

- **Type**: Real
- **Description**: mixing parameter: 0 means no new charge
- **Default**: 0.7

#### mixing_ndim

- **Type**: Integer
- **Description**: It indicates the mixing dimensions in Pulay, Pulay method use the density from previous mixing_ndim steps and do a charge mixing based on these density.
- **Default**: 8

#### mixing_gg0

- **Type**: Real
- **Description**: used in pulay-kerker mixing method
- **Default**: 1.5

#### gamma_only

- **Type**: Integer
- **Description**: It is an important parameter **only to be used in localized orbitals set**.
    It you set gamma_only = 1, ABACUS use gamma only, the algorithm is fast and you don't need to specify the k-points file. If you set gamma_only = 0, more than one k-point is used and the ABACUS is slower compared to gamma only algorithm.
- **Default**: 0

#### printe

- **Type**: Integer
- **Description**: Print out energy for each band for every printe steps
- **Default**: 100

#### scf_nmax

- **Type**: Integer
- **Description**:This variable indicates the maximal iteration number for electronic iterations.
- **Default**: 40

#### scf_thr

- **Type**: Real
- **Description**: An important parameter in ABACUS. It`s the threshold for electronic iteration. It represents the charge density error between two sequential density from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
- **Default**:1e-06

#### chg_extrap

- **Type**: String
- **Description**: Methods to do extrapolation of density when ABACUS is doing geometry relaxations.
  - atomic: atomic extrapolation
  - first-order: first-order extrapolation
  - second-order: second-order extrapolation
- **Default**:atomic

 by time(NULL).
- **Default**:0

### Electronic structure (SDFT)

This part of variables are used to control the parameters of stochastic DFT (SDFT),  mix stochastic-deterministic DFT (MDFT), or complete-basis Chebyshev method (CT). To use it,  [calculation](#calculation) need to be set to "sto-scf" and "sto-md". We suggest using SDFT to calculate high-temperature systems and we only support [smearing_method](#smearing_method) "fd".

#### method_sto

- **Type**: Integer
- **Description**:
  - Different method to do SDFT.
  - 1: SDFT calculates $T_n(\hat{h})\ket{\chi}$ twice, where $T_n(x)$ is the n-th order Chebyshev polynomial and $\hat{h}=\frac{\hat{H}-\bar{E}}{\Delta E}$ owning eigen-value $\in(-1,1)$. This method cost less memory but slow.
  - 2: SDFT calculates $T_n(\hat{h})\ket{\chi}$ once but need much more memory. This method is much faster. Besides, it calculate $N_e$ with $\bra{\chi}\sqrt{\hat f}\sqrt{\hat f}\ket{\chi}$, which needs smaller [nche_sto](#nche_sto). However, when memory is not enough, only method 1 can be used.
  - other: use 2
- **Default**: 2

#### nbands_sto

- **Type**: Integer
- **Description**:
  - nbands_sto>0: Number of stochastic orbitals to calculate in SDFT and MDFT.  More bands obtain more precise results or smaller stochastic errors ($ \propto 1/\sqrt{N_{\chi}}$);
  - nbands_sto=0: Complete basis will be used to replace stochastic orbitals with the Chebyshev method (CT) and it will get the results the same as KSDFT without stochastic errors.
  - If you want to do MDFT. [nbands](#nbands) which represents the number of KS orbitals should be set.
- **Default**: 256

#### nche_sto

- **Type**: Integer
- **Description**: Chebyshev expansion orders for SDFT, MDFT, CT methods.
- **Default**:100

#### emin_sto

- **Type**: Real
- **Description**: Trial energy to guess the lower bound of eigen energies of the Hamitonian Operator $\hat{H}$. The unit is Ry.
- **Default**:0.0

#### emax_sto

- **Type**: Real
- **Description**: Trial energy to guess the upper bound of eigen energies of the Hamitonian Operator $\hat{H}$. The unit is Ry.
- **Default**:0.0

#### seed_sto

- **Type**: Integer
- **Description**: The random seed to generate stochastic orbitals.
  - seed_sto>=0: Stochastic orbitals have the form of $\exp(i2\pi\theta(G))$, where $\theta$ is a uniform distribution in $(0,1)$. If seed_sto = 0, the seed is decided by time(NULL).
  - seed_sto<=-1: Stochastic orbitals have the form of $\pm1$ with the equal probability. If seed_sto = -1, the seed is decided by time(NULL).
- **Default**:0

#### initsto_freq

- **Type**: Integer
- **Description**: Frequency (once each initsto_freq steps) to generate new stochastic orbitals when running md.
- **Default**:1000

#### npart_sto

- **Type**: Integer
- **Description**: Make memory cost to 1/npart_sto times of previous one when running post process of SDFT like DOS with method_sto = 2.
- **Default**:1

### Geometry relaxation

This part of variables are used to control the geometry relaxation.

#### relax_nmax

- **Type**: Integer
- **Description**: The maximal number of ionic iteration steps, the minimal value is 1.
- **Default**: 1

#### cal_force

- **Description**: If set to 1, calculate the force at the end of the electronic iteration. 0 means the force calculation is turned off.
- **Default**: 0

#### force_thr

- **Type**: Real
- **Description**: The threshold of the force convergence, it indicates the largest force among all the atoms, the unit is Ry=Bohr
- **Default**: 0.000388935 Ry/Bohr = 0.01 eV/Angstrom

#### force_thr_ev

- **Type**: Real
- **Description**: The threshold of the force convergence, has the same function as force_thr, just the unit is different, it is eV=Angstrom, you can choose either one as you like. The recommendation value for using atomic orbitals is 0:04 eV/Angstrom.
- **Default**: 0.01 eV/Angstrom

#### relax_bfgs_w1

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.01

#### relax_bfgs_w2

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.5

#### relax_bfgs_rmax

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the maximal movement of all the atoms. The sum of the movements from all atoms can be increased during the optimization steps. However, it will not be larger than relax_bfgs_rmax Bohr.
- **Default**: 0.8

#### relax_bfgs_rmin

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the minimal movement of all the atoms. When the movement of all the atoms is smaller than relax_bfgs_rmin Bohr , and the force convergence is still not achieved, the calculation will break down.
- **Default**: 1e-5

#### relax_bfgs_init

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the initial movement of all the atoms. The sum of the movements from all atoms is relax_bfgs_init Bohr.
- **Default**: 0.5

#### cal_stress

- **Type**: Integer
- **Description**: If set to 1, calculate the stress at the end of the electronic iteration. 0 means the stress calculation is turned off.
- **Default**: 0

#### stress_thr

- **Type**: Real
- **Description**: The threshold of the stress convergence, it indicates the largest stress among all the directions, the unit is KBar,
- **Default**: 0.01

#### press1, press2, press3

- **Type**: Real
- **Description**: the external pressures along three axes,the compressive stress is taken to be positive, the unit is KBar.
- **Default**: 0

#### fixed_axes

- **Type**: String
- **Description**:which axes are fixed when do cell relaxation. Possible choices are:
  - None : default; all can relax
  - volume : relaxation with fixed volume
  - a : fix a axis during relaxation
  - b : fix b axis during relaxation
  - c : fix c axis during relaxation
  - ab : fix both a and b axes during relaxation
  - ac : fix both a and c axes during relaxation
  - bc : fix both b and c axes during relaxation
  - abc : fix all three axes during relaxation
- **Default**: None

#### relax_method

- **Type**: String
- **Description**: The method to do geometry optimizations. If set to bfgs, using BFGS algorithm. If set to cg, using cg algorithm. If set to sd, using steepest-descent lgorithm.
- **Default**: cg

#### relax_cg_thr

- **Type**: Real
- **Description**: When move-method is set to 'cg-bfgs', a mixed cg-bfgs algorithm is used. The ions first move according to cg method, then switched to bfgs when maximum of force on atoms is reduced below cg-threshold. Unit is eV/Angstrom.
- **Default**: 0.5

#### cell_factor

- **Type**: Real
- **Description**: Used in the construction of the pseudopotential tables. It should exceed the maximum linear contraction of the cell during a simulation.
- **Default**: 1.2

### Variables related to output information

This part of variables are used to control the output of properties.

#### out_force

- **Type**: Integer
- **Description**: Determines whether to output the out_force into a file named `Force.dat` or not. If 1, then force will be written; if 0, then the force will not be written.
- **Default**: 0

#### out_mul

- **Type**: Integer
- **Description**: If set to 1, ABACUS will output the Mulliken population analysis result. The name of the output file is mulliken.txt
- **Default**: 0

#### out_freq_elec

- **Type**: Integer
- **Description**: If set to >1, it represents the frequency of electronic iters to output charge density (if [out_chg](#out_chg) is turned on) and wavefunction (if [out_wfc_pw](#out_wf_pw) or [out_wfc_r](#out_wfc_r) is turned on). If set to 0, ABACUS will output them only when converged in SCF. Used for the restart of SCF.
- **Default**: 0

#### out_freq_ion

- **Type**: Integer
- **Description**: If set to >1, it represents the frequency of ionic steps to output charge density (if [out_chg](#out_chg) is turned on) and wavefunction (if [out_wfc_pw](#out_wf_pw) or [out_wfc_r](#out_wfc_r) is turned on). If set to 0, ABACUS will output them only when ionic steps reach its maximum step. Used for the restart of MD or Relax.
- **Default**: 0

#### out_chg

- **Type**: Integer
- **Description**: If set to 1, ABACUS will output the charge density on real space grid. The name of the density file is SPIN1_CHGCAR and SPIN2_CHGCAR (if nspin = 2). Suppose each density on grid has coordinate (x; y; z). The circle order of the density on real space grid is: z is the outer loop, then y and finally x (x is moving fastest).
- **Default**: 0

#### out_pot

- **Type**: Integer
- **Description**: If set to 1, ABACUS will output the local potential on real space grid. The name of the file is SPIN1_POT and SPIN2_POT (if nspin = 2). If set to 2, ABACUS will output the electrostatic potential on real space grid. The name of the file is ElecStaticPot and ElecStaticP ot_AV E (along the z-axis).
- **Default**: 0

#### out_dm

- **Type**: Integer
- **Description**: If set to 1, ABACUS will output the density matrix of localized orbitals, only useful for localized orbitals set. The name of the output file is SPIN1_DM and SPIN2_DM in the output directory.
- **Default**: 0

#### out_wfc_pw

- **Type**: Integer
- **Description**: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs the coefficients of wave functions into text files. The file names are WAVEFUNC$K.txt, where $K is the index of k point. When set this variable to 2, results are stored in binary files. The file names are WAVEFUNC$K.dat.
- **Default**: 0

#### out_wfc_r

- **Type**: Integer
- **Description**: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs real-space wave functions into  `OUT.suffix/wfc_realspace/`. The file names are wfc_realspace$K$B, where $K is the index of k point, $B is the index of band.
- **Default**: 0

#### out_wfc_lcao

- **Type**: Integer
- **Description**: **Only used in localized orbitals set**. If set to 1, ABACUS will output the wave functions coefficients.
- **Default**: 0

#### out_dos

- **Type**: Integer
- **Description**: Controls whether to output the density of state (DOS). For more information, refer to the [worked example](examples/dos.md).
- **Default**: 0

#### out_band

- **Type**: Integer
- **Description**: Controls whether to output the band structure. For mroe information, refer to the [worked example](examples/band-struc.md)
- **Default**: 0

#### out_proj_band

- **Type**: Integer
- **Description**: Controls whether to output the projected band structure. For mroe information, refer to the [worked example](examples/band-struc.md)
- **Default**: 0

#### out_stru

- **Type**: Boolean
- **Description**: If set to 1, then tje structure files will be written after each ion step
- **Default**: 0

#### out_level

- **Type**: String
- **Description**: Controls the level of output. `ie` means write output at electron level; `i` means write additional output at ions level.
- **Default**: ie

#### out_alllog

- **Type**: Integer
- **Description**: determines whether to write log from all ranks in an MPI run. If set to be 1, then each rank will write detained running information to a file named running_${calculation}\_(${rank}+1).log. If set to 0, log will only be written from rank 0 into a file named running_${calculation}.log.
- **Default**: 0

#### out_mat_hs

- **Type**: Boolean
- **Description**: Only for LCAO calculations. When set to 1, ABACUS will generate two lists of files `data-$k-H` and `data-$k-S` that store the Hamiltonian and S matrix for each k point in k space, respectively.
- **Default**: 0

#### out_mat_r

- **Type**: Boolean
- **Description**: Only for LCAO and not gamma_only calculations. When set to 1, ABACUS will generate a file with name staring with `data-rR-tr` which stores overlap matrix as a function of R, in units of lattice vectors.
- **Default**: 0

#### out_mat_hs2

- **Type**: Boolean
- **Description**: Only for LCAO and not gamma_only calculations. When set to 1, ABACUS will generate two files starting with `data-HR-sparse` and `data-SR-sparse` that store the Hamiltonian and S matrix in real space, respectively, as functions of R, in units of lattice vectors.
- **Default**: 0

#### out_element_info

- **Type**: Boolean
- **Description**: When set to 1, ABACUS will generate a new directory under OUT.suffix path named as element name such as 'Si', which contained files "Si-d1-orbital-dru.dat  Si-p2-orbital-k.dat    Si-s2-orbital-dru.dat
    Si-d1-orbital-k.dat    Si-p2-orbital-r.dat    Si-s2-orbital-k.dat
    Si-d1-orbital-r.dat    Si-p2-orbital-ru.dat   Si-s2-orbital-r.dat
    Si-d1-orbital-ru.dat   Si-p-proj-k.dat        Si-s2-orbital-ru.dat
    Si.NONLOCAL            Si-p-proj-r.dat        Si-s-proj-k.dat
    Si-p1-orbital-dru.dat  Si-p-proj-ru.dat       Si-s-proj-r.dat
    Si-p1-orbital-k.dat    Si-s1-orbital-dru.dat  Si-s-proj-ru.dat
    Si-p1-orbital-r.dat    Si-s1-orbital-k.dat    v_loc_g.dat
    Si-p1-orbital-ru.dat   Si-s1-orbital-r.dat
Si-p2-orbital-dru.dat  Si-s1-orbital-ru.dat" for example.

- **Default**: 0

#### restart_save

- **Type**: Boolean
- **Description**: Only for LCAO, store charge density file and H matrix file every scf step for restart.
- **Default**: 0

#### restart_load

- **Type**: Boolean
- **Description**: Only for LCAO, used for restart, only if that:
  - set restart_save as true and do scf calculation before.
  - please ensure suffix is same with calculation before and density file and H matrix file is exist.
    restart from stored density file and H matrix file.
- **Default**: 0

### Density of states

This part of variables are used to control the calculation of DOS.

#### dos_edelta_ev

- **Type**: Real
- **Description**: controls the step size in writing DOS (in eV).
- **Default**: 0.1

#### dos_sigma

- **Type**: Real
- **Description**: controls the width of Gaussian factor when obtaining smeared DOS (in eV).
- **Default**: 0.07

#### dos_scale

- **Type**: Real
- **Description**: the energy range of dos output is given by (emax-emin)*(1+dos_scale), centered at (emax+emin)/2. This parameter will be used when dos_emin and dos_emax are not set.
- **Default**: 0.01

#### dos_emin_ev

- **Type**: Real
- **Description**: minimal range for dos (in eV). If we set it, "dos_scale" will be ignored.
- **Default**: minimal eigenenergy of $\hat{H}$

#### dos_emax_ev

- **Type**: Real
- **Description**: maximal range for dos (in eV). If we set it, "dos_scale" will be ignored.
- **Default**: maximal eigenenergy of $\hat{H}$

#### dos_nche

- **Type**: Integer
- **Description**: orders of Chebyshev expansions when using SDFT to calculate DOS
- **Default**: 100

### DeePKS

This part of variables are used to control the usage of DeePKS method (a comprehensive data-driven approach to improve accuracy of DFT).
Warning: this function is not robust enough for the current version. Please try the following variables at your own risk:

#### deepks_out_labels

- **Type**: Boolean
- **Description**: when set to 1, ABACUS will calculate and output descriptor for DeePKS training. In `LCAO` calculation, a path of *.orb file is needed to be specified under `NUMERICAL_DESCRIPTOR`in `STRU`file. For example:

    ```txt
    NUMERICAL_ORBITAL
    H_gga_8au_60Ry_2s1p.orb
    O_gga_7au_60Ry_2s2p1d.orb

    NUMERICAL_DESCRIPTOR
    jle.orb
    ```

- **Default**: 0

#### deepks_descriptor_lmax

- **Type**: Integer
- **Description**: when generating projectors, this variable controls the max angular momentum of descriptor basis.
- **Default**: 2

#### deepks_scf

- **Type**: Boolean
- **Description**: only when deepks is enabled in `LCAO` calculation can this variable set to 1. Then, a trained, traced model file is needed for self-consistant field iteration in DeePKS method.
- **Default**: 0

#### deepks_model

- **Type**: String
- **Description**: the path of the trained, traced NN model file (generated by deepks-kit). used when deepks_scf is set to 1.
- **Default**: None

### Electric field and dipole correction

This part of variables are relevant to electric field and dipole correction

#### efield_flag

- **Type**: Boolean
- **Description**: If set to true, a saw-like potential simulating an electric field
is added to the bare ionic potential.
- **Default**: false

#### dip_cor_flag

- **Type**: Boolean
- **Description**: If dip_cor_flag == true and efield_flag == true,  a dipole correction is also
added to the bare ionic potential. If you want no electric field, parameter efield_amp  should be zero. Must be used ONLY in a slab geometry for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
- **Default**: false

#### efield_dir

- **Type**: Integer
- **Description**: The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir = 0, 1 or 2. Used only if efield_flag == true.
- **Default**: 2

#### efield_pos_max

- **Type**: Real
- **Description**: Position of the maximum of the saw-like potential along crystal axis efield_dir, within the  unit cell, 0 < efield_pos_max < 1. Used only if efield_flag == true.
- **Default**: 0.5

#### efield_pos_dec

- **Type**: Real
- **Description**: Zone in the unit cell where the saw-like potential decreases, 0 < efield_pos_dec < 1. Used only if efield_flag == true.
- **Default**: 0.1

#### efield_amp

- **Type**: Real
- **Description**: Amplitude of the electric field, in ***Hartree*** a.u.; 1 a.u. = 51.4220632*10^10 V/m. Used only if efield_flag == true. The saw-like potential increases with slope efield_amp  in the region from (efield_pos_max+efield_pos_dec-1) to (efield_pos_max), then decreases until (efield_pos_max+efield_pos_dec), in units of the crystal vector efield_dir. Important: the change of slope of this potential must be located in the empty region, or else unphysical forces will result.
- **Default**: 0.0

### Exact Exchange

This part of variables are relevant when using hybrid functionals

#### exx_hybrid_alpha

- **Type**: Real
- **Description**: fraction of Fock exchange in hybrid functionals, so that $E_{X}=\alpha F_{X}+(1-\alpha)E_{X,LDA/GGA}$
- **Default**: 0.25

#### exx_hse_omega

- **Type**: Real
- **Description**: range-separation parameter in HSE functional, such that $1/r=erfc(\omega r)/r+erf(\omega r)/r$.
- **Default**: 0.11

#### exx_separate_loop

- **Type**: Boolean
- **Description**: There are two types of iterative approach provided by ABACUS to evaluate Fock exchange. If this parameter is set to 0, it will start with a GGA-Loop, and then Hybrid-Loop, in which EXX Hamiltonian $H_{exx}$ is updated with electronic iterations. If this parameter is set to 1, a two-step method is employed, i.e. in the inner iterations, density matrix is updated, while in the outer iterations, $H_{exx}$ is calculated based on density matrix that converges in the inner iteration. 
- **Default**: 1

#### exx_hybrid_step

- **Type**: Integer
- **Description**: This variable indicates the maximal electronic iteration number in the evaluation of Fock exchange.
- **Default**: 100

#### exx_lambda

- **Type**: Real
- **Description**: It is used to compensate for divergence points at G=0 in the evaluation of Fock exchange using *lcao_in_pw* method.
- **Default**: 0.3

#### exx_pca_threshold

- **Type**: Real
- **Description**: To accelerate the evaluation of four-center integrals ($ik|jl$), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$. The size of the ABF (i.e. number of $P_{k}$) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculations becomes faster. However this comes at the cost of computational accuracy. A relatively safe choice of the value is 1d-4.
- **Default**: 0

#### exx_c_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). Smaller components (less than exx_c_threshold) of the $C^{k}_{ij}$ matrix is neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-4.
- **Default**: 0

#### exx_v_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). With the approximation $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$, the four-center integral in Fock exchange is expressed as $(ik|jl)=\Sigma_{a,b}C^{a}_{ij}V_{ab}C^{b}_{kl}$, where $V_{ab}=(P_{a}|P_{b})$ is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 0

#### exx_dm_threshold

- **Type**: Real
- **Description**: The Fock exchange can be expressed as $\Sigma_{k,l}(ik|jl)D_{kl}$ where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-4.
- **Default**: 0

#### exx_schwarz_threshold

- **Type**: Real
- **Description**: In practice the four-center integrals are sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each integral before carrying out explicit evaluations. Those that are smaller than exx_schwarz_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-5.
- **Default**: 0

#### exx_cauchy_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-7.
- **Default**: 0

#### exx_ccp_threshold

- **Type**: Real
- **Description**: It is related to the cutoff of on-site Coulomb potentials, currently not used.
- **Default**: 1e-8

#### exx_ccp_rmesh_times

- **Type**: Real
- **Description**: This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. For HSE1, setting it to 1 is enough. But for PBE0, a much larger number must be used.
- **Default**: 10

#### exx_distribute_type

- **Type**: String
- **Description**: When running in parallel, the evaluation of Fock exchange is done by distributing atom pairs on different threads, then gather the results. exx_distribute_type governs the mechanism of distribution. Available options are `htime`, `order`, `kmean1` and `kmeans2`. `order` is where atom pairs are simply distributed by their orders. `hmeans` is a distribution where the balance in time is achieved on each processor, hence if the memory is sufficient, this is the recommended method. `kmeans1` and `kmeans2` are two methods where the k-means clustering method is used to reduce memory requirement. They might be necessary for very large systems.
- **Default**: `htime`

#### exx_opt_orb_lmax

- **Type**: Integer
- **Description**: See also the entry [exx_hybrid_type](#exx_hybrid_type). This parameter is only relevant when exx_hybrid_type=`opt_orb`. The radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. exx_opt_orb_lmax gives the maximum l of the spherical Bessel functions. A reasonable choice is 2.
- **Default**: 0

#### exx_opt_orb_ecut

- **Type**: Real
- **Description**: See also the entry [exx_hybrid_type](#exx_hybrid_type). This parameter is only relevant when exx_hybrid_type=`opt_orb`. A plane wave basis is used to optimize the radial ABFs. This parameter thus gives the cut-off of plane wave expansion, in Ry. A reasonable choice is 60.
- **Default**: 0

#### exx_opt_orb_tolerence

- **Type**: Real
- **Description**: See also the entry [exx_hybrid_type](#exx_hybrid_type). This parameter is only relevant when exx_hybrid_type=`opt_orb`. exx_opt_orb_tolerence determines the threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.
- **Default**: 0

### Molecular dynamics

This part of variables are used to control the molecular dynamics calculations.

#### md_type

- **Type**: Integer
- **Description**: control the ensemble to run md.
  - -1: FIRE method to relax;
  - 0: NVE ensemble;
  - 1: NVT ensemble with Nose Hoover Chain;
  - 2: NVT ensemble with Langevin method;
  - 3: NVT ensemble with Anderson thermostat;
  - 4: MSST method;

  ***Note: when md_type is set to 1, md_tfreq is required to stablize temperature. It is an empirical parameter whose value is system-dependent, ranging from 1/(40\*md_dt) to 1/(100\*md_dt). An improper choice of its value might lead to failure of job.***
- **Default**: 1

#### md_nstep

- **Type**: Integer
- **Description**: the total number of md steps.
- **Default**: 10

#### md_ensolver

- **Type**: String
- **Description**: choose the energy solver for MD.
  - FP: First-Principles MD;
  - LJ: Leonard Jones potential;
  - DP: DeeP potential;
- **Default**: FP

#### md_restart

- **Type**: Boolean
- **Description**: to control whether restart md.
  - 0: When set to 0, ABACUS will calculate md normolly.
  - 1: When set to 1, ABACUS will calculate md from last step in your test before.
- **Default**: 0

#### md_dt

- **Type**: Double
- **Description**: This is the time step(fs) used in md simulation .
- **Default**: 1.0

#### md_tfirst & md_tlast

- **Type**: Double
- **Description**: This is the temperature (K) used in md simulation, md_tlast`s default value is md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.
- **Default**: No default

#### md_dumpfreq

- **Type**: Integer
- **Description**:This is the frequence to dump md information.
- **Default**: 1

#### md_restartfreq

- **Type**: Integer
- **Description**:This is the frequence to output restart information.
- **Default**: 5

#### md_seed

- **Type**: Integer
- **Description**:
  - md_seed < 0: No srand() in MD initialization.
  - md_seed >= 0: srand(md_seed) in MD initialization.
- **Default**: -1

#### md_tfreq

- **Type**: Real
- **Description**:
  - When md_type = 1, md_tfreq controls the frequency of the temperature oscillations during the simulation. If it is too large, the
temperature will fluctuate violently; if it is too small, the temperature will take a very long time to equilibrate with the atomic system.
  - When md_type = 3, md_tfreq*md_dt is the collision probability in Anderson method.
  - If md_tfreq is not set in INPUT, md_tfreq will be autoset to be 1/40/md_dt.
- **Default**: 1/40/md_dt

#### md_mnhc

- **Type**: Integer
- **Description**: Number of Nose-Hoover chains.
- **Default**: 4

#### lj_rcut

- **Type**: Real
- **Description**: Cut-off radius for Leonard Jones potential (angstrom).
- **Default**: 8.5 (for He)

#### lj_epsilon

- **Type**: Real
- **Description**: The value of epsilon for Leonard Jones potential (eV).
- **Default**: 0.01032 (for He)

#### lj_sigma

- **Type**: Real
- **Description**: The value of sigma for Leonard Jones potential (angstrom).
- **Default**: 3.405 (for He)

#### msst_direction

- **Type**: Integer
- **Description**: the direction of shock wave for MSST.
- **Default**: 2 (z direction)

#### msst_vel

- **Type**: Real
- **Description**: the velocity of shock wave ($\AA$/fs) for MSST.
- **Default**: 0.0

#### msst_vis

- **Type**: Real
- **Description**: artificial viscosity (mass/length/time) for MSST.
- **Default**: 0.0

#### msst_tscale

- **Type**: Real
- **Description**: reduction in initial temperature (0~1) used to compress volume in MSST.
- **Default**: 0.01

#### msst_qmass

- **Type**: Real
- **Description**: Inertia of extended system variable. Used only when md_type is 4, you should set a number which is larger than 0. Note that Qmass of NHC is set by md_tfreq.
- **Default**: No default

#### md_damp

- **Type**: Real
- **Description**: damping parameter (fs) used to add force in Langevin method.
- **Default**: 1.0

### DFT+U correction

This part of variables are used to control DFT+U correlated parameters

#### dft_plus_u

- **Type**: Boolean
- **Description**: If set to 1, ABCUS will calculate plus U correction, which is especially important for correlated electron.
- **Default**: 0

#### orbital_corr

- **Type**: Int
- **Description**: $l_1,l_2,l_3,\ldots$ for atom type 1,2,3 respectively.(usually 2 for d electrons and 3 for f electrons) .Specify which orbits need plus U correction for each atom. If set to -1, the correction would not be calculate for this atom.
- **Default**: None

#### hubbard_u

- **Type**: Real
- **Description**: Hubbard Coulomb interaction parameter U(ev) in plus U correction,which should be specified for each atom unless Yukawa potential is use. ABACUS use a simplified scheme which only need U and J for each atom.
- **Default**: 0.0

#### hund_j

- **Type**: Real
- **Description**: Hund exchange parameter J(ev) in plus U correction ,which should be specified for each atom unless Yukawa potential is use. ABACUS use a simplified scheme which only need U and J for each atom.
- **Default**: 0.0

#### yukawa_potential

- **Type**: Boolean
- **Description**: whether use the local screen Coulomb potential method to calculate the value of U and J. If this is set to 1, hubbard_u and hund_j do not need to be specified.
- **Default**: 0

#### omc

- **Type**: Boolean
- **Description**: whether turn on occupation matrix control method or not
- **Default**: 0

### vdW correction

This part of variables are used to control vdW-corrected related parameters.

#### vdw_method

- **Type**: String
- **Description**: If set to d2 ,d3_0 or d3_bj, ABACUS will calculate corresponding vdW correction, which is DFT-D2, DFT-D3(0) or DFTD3(BJ) method. And this correction includes energy and forces. `none` means that no vdW-corrected method has been used.
- **Default**: none

#### vdw_s6

- **Type**: Real
- **Description**: This scale factor is to optimize the interaction energy deviations. For DFT-D2, it is found to be 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). For DFT-D3, recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
- **Default**: 0.75 if vdw_method is chosen to be d2; 1.0 if vdw_method is chosen to be d3_0 or d3_bj

#### vdw_s8

- **Type**: Real
- **Description**: This scale factor is only the parameter of DFTD3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
- **Default**: 0.722 if vdw_method is chosen to be d3_0; 0.7875 if vdw_method is chosen to be d3_bj

#### vdw_a1

- **Type**: Real
- **Description**: This damping function parameter is only the parameter of DFT-D3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
- **Default**: 1.217 if vdw_method is chosen to be d3_0; 0.4289 if vdw_method is chosen to be d3_bj

#### vdw_a2

- **Type**: Real
- **Description**: This damping function arameter is only the parameter of DFT-D3(BJ) approach. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
- **Default**: 1.0 if vdw_method is chosen to be d3_0; 4.4407 if vdw_method is chosen to be d3_bj

#### vdw_d

- **Type**: Real
- **Description**: The variable is to control the dumping speed of dumping function of DFT-D2.
- **Default**: 20

#### vdw_abc

- **Type**: Integer
- **Description**: The variable is to control the calculation of three-body term of DFT-D3 approachs, including D3(0) and D3(BJ). If set to 1, ABACUS will calculate three-body term, otherwise, the three-body term is not included.
- **Default**: 0

#### vdw_C6_file

- **Type**: String
- **Description**: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent's $C_6$ Parameters file. If you dont setup this, ABACUS will use the default $C_6$ Parameters stored in the programme already. We've stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new $C_6$ Parameters, you should provide a file contained all the $C_6$ Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
- **Default**: default

#### vdw_C6_unit

- **Type**: String
- **Description**: This variable which is useful only when set vdw_method to d2 specifies unit of $C_6$ Parameters. Two kinds of unit is available: `J·nm6/mol` (means J·nm^{6}/mol) and eVA(means eV·Angstrom)
- **Default**: Jnm6/mol

#### vdw_R0_file

- **Type**: String
- **Description**: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent's $R_0$ Parameters file. If you don't setup this, ABACUS will use the default $R_0$ Parameters stored in the programme already. We've stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new $R_0$ Parameters, you should provide a file contained all the $R_0$ Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
- **Default**: default

#### vdw_R0_unit

- **Type**: String
- **Description**: This variable which is useful only when set vdw_method to d2 specifies unit of $R_0$ Parameters. Two kinds of unit is available: A(means Angstrom) and Bohr.
- **Default**: A

#### vdw_model

- **Type**: String
- **Description**: To calculate the periodic structure, you can assign the number of lattice cells calculated. This variable specifies the kind of model to assign. If set to period, ABACUS will calculate a cubic periodic structure as assigned in vdw_period. If set to radius, ABACUS will calculate a cubic periodic structure containing a sphere, whose radius is vdw_radius and centre is origin point.
- **Default**: radius

#### vdw_radius

- **Type**: Real
- **Description**: If vdw_model is set to radius, this variable specifies the radius of the calculated sphere. For DFT-D2, the default value is 56.6918 ,while it is 95 for DFT-D3. This variable will be set to default, if no vdW-corrected method has been used.
- **Default**: 56.6918 if vdw_method is chosen to be d2; 95 if vdw_method is chosen to be d3_0 or d3_bj

#### vdw_radius_unit

- **Type**: String
- **Description**: If vdw_model is set to radius, this variable specifies the unit of vdw_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
- **Default**: Bohr

#### vdw_cn_radius

- **Type**: Real
- **Description**: This cutoff is chosen for the calculation of the coordination number (CN) in DFT-D3 approachs, including D3(0) and D3(BJ).
- **Default**: 40.0

#### vdw_cn_radius_unit

- **Type**: String
- **Description**: This variable specifies the unit of vdw_cn_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
- **Default**: Bohr

#### vdw_period

- **Type**: Int Int Int
- **Description**: If vdw_model is set to period, these variables specify the number of x, y and z periodic.
- **Default**: 3 3 3

### Berry phase and wannier90 interface

This part of variables are used to control berry phase and wannier90 interfacae parameters.

#### berry_phase

- **Type**: Integer
- **Description**: 1, calculate berry phase; 0, no calculate berry phase.
- **Default**: 0

#### gdir

- **Type**: Integer
- **Description**:
  - 1: calculate the polarization in the direction of the lattice vector a_1 that is defined in STRU file.
  - 2: calculate the polarization in the direction of the lattice vector a_2 that is defined in STRU file.
  - 3: calculate the polarization in the direction of the lattice vector a_3 that is defined in STRU file.
- **Default**: 3

#### towannier90

- **Type**: Integer
- **Description**: 1, generate files for wannier90 code; 0, no generate.
- **Default**: 0

#### nnkpfile

- **Type**: String
- **Description**: the file name when you run “wannier90 -pp ...”.
- **Default**: seedname.nnkp

#### wannier_spin

- **Type**: String
- **Description**: If nspin is set to 2,
  - up: calculate spin up for wannier function.
  - down: calculate spin down for wannier function.
- **Default**: up

### TDDFT: time dependent density functional theory

#### tddft

- **Type**: Integer
- **Description**:
  - 1: calculate the real time time dependent density functional theory (TDDFT).
  - 0: do not calculate TDDFT.
- **Default**: 0

#### td_scf_thr

- **Type**: Double
- **Description**: Accuracy of electron convergence when doing time-dependent evolution.
- **Default**: 1e-9

#### td_dt

- **Type**: Double
- **Description**: Time-dependent evolution time step. (fs)
- **Default**: 0.02

#### td_force_dt

- **Type**: Double
- **Description**: Time-dependent evolution force changes time step. (fs)
- **Default**: 0.02

#### td_vext

- **Type**: Integer
- **Description**:
  - 1: add a laser material interaction (extern laser field).
  - 0: no extern laser field.
- **Default**: 0

#### td_vext_dire

- **Type**: Integer
- **Description**:
  - 1: the direction of external light field is along x axis.
  - 2: the direction of external light field is along y axis.
  - 3: the direction of external light field is along z axis.
- **Default**: 1

#### td_timescale

- **Type**: Double
- **Description**: Time range of external electric field application. (fs)
- **Default**: 0.5

#### td_vexttype

- **Type**: Integer
- **Description**:
  - 1: Gaussian-type light field.
  - 2: Delta function form light field.
  - 3: Trigonometric function form light field.
- **Default**: 1

#### td_vextout

- **Type**: Integer
- **Description**:
  - 1: Output external electric field.
  - 0: do not Output external electric field.
- **Default**: 0

#### td_dipoleout

- **Type**: Integer
- **Description**:
  - 1: Output dipole.
  - 0: do not Output dipole.
- **Default**: 0

#### ocp

- **Type**: Boolean
- **Description**: option for choose whether calcualting constrained DFT or not.
    Only used for TDDFT.
- **Default**:0

#### ocp_set

- **Type**: string
- **Description**: If ocp is true, the ocp_set is a string to set the number of occupancy, like 1 10 * 1 0 1 representing the 13 band occupancy, 12th band occupancy 0 and the rest 1, the code is parsing this string into an array through a regular expression.
- **Default**:none

### Variables useful for debugging

#### nurse

- **Type**: Boolean
- **Description**: If set to 1, the Hamiltonian matrix and S matrix in each iteration will be written in output.
- **Default**: 0

#### t_in_h

- **Type**: Boolean
- **Description**: If set to 0, then kinetic term will not be included in obtaining the Hamiltonian.
- **Default**: 1

#### vl_in_h

- **Type**: Boolean
- **Description**: If set to 0, then local pseudopotential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

#### vnl_in_h

- **Type**: Boolean
- **Description**:  If set to 0, then non-local pseudopotential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

#### test_force

- **Type**: Boolean
- **Description**: If set to 1, then detailed components in forces will be written to output.
- **Default**: 0

#### test_stress

- **Type**: Boolean
- **Description**: If set to 1, then detailed components in stress will be written to output.
- **Default**: 0

#### colour

- **Type**: Boolean
- **Description**: If set to 1, output to terminal will have some color.
- **Default**: 0

### Electronic conductivities

Frequency-dependent electronic conductivities can be calculated with Kubo-Greenwood formula[Phys. Rev. B 83, 235120 (2011)].
Onsager coefficiencies:
$$L_{mn}(\omega)=(-1)^{m+n}\frac{2\pi e^2\hbar^2}{3m_e^2\omega\Omega}\\
\times\sum_{ij\alpha\mathbf{k}}W(\mathbf{k})\left(\frac{\epsilon_{i\mathbf{k}}+\epsilon_{j\mathbf{k}}}{2}-\mu\right)^{m+n-2}|
\langle\Psi_{i\mathbf{k}}|\nabla_\alpha|\Psi_{j\mathbf{k}}\rangle|^2\\
\times[f(\epsilon_{i\mathbf{k}})-f(\epsilon_{j\mathbf{k}})]\delta(\epsilon_{j\mathbf{k}}-\epsilon_{i\mathbf{k}}-\hbar\omega).$$
They can also computed by $j$-$j$ correlation function.
$$L_{mn}=\frac{2e^{m+n-2}}{3\Omega\hbar\omega}\Im[\tilde{C}_{mn}(\omega)]\\
\tilde{C}_{mn}=\int_0^\infty C_{mn}(t)e^{-i\omega t}e^{-\frac{1}{2}(\Delta E)^2t^2}dt\\
C_{mn}(t)=-2\theta(t)\Im\left\{Tr\left[\sqrt{\hat f}\hat{j}_m(1-\hat{f})e^{i\frac{\hat{H}}{\hbar}t}\hat{j}_ne^{-i\frac{\hat{H}}{\hbar}t}\sqrt{\hat f}\right]\right\},$$
where $j_1$ is electric flux and $j_2$ is thermal flux.
Frequency-dependent electric conductivities:    $\sigma(\omega)=L_{11}(\omega).$
Frequency-dependent thermal conductivities: $\kappa(\omega)=\frac{1}{e^2T}\left(L_{22}-\frac{L_{12}^2}{L_{11}}\right).$
DC electric conductivities: $\sigma = \lim_{\omega\to 0}\sigma(\omega)$
Thermal conductivities: $\kappa = \lim_{\omega\to 0}\kappa(\omega)$

#### cal_cond

- **Type**: Boolean
- **Description**: If set to 1, electronic conductivities will be calculated. Only supported in calculations of SDFT and KSDFT_PW. 
- **Default**: 0

#### cond_nche

- **Type**: Integer
- **Description**: Chebyshev expansion orders for stochastic Kubo Greenwood. Only used when the calculation is SDFT.
- **Default**: 20

#### cond_dw

- **Type**: Real
- **Description**: Frequency interval ($d\omega$) for frequency-dependent conductivities. The unit is eV.
- **Default**: 0.1

#### cond_wcut

- **Type**: Real
- **Description**: Cutoff frequency for frequency-dependent conductivities. The unit is eV.
- **Default**: 10.0

#### cond_wenlarge

- **Type**: Integer
- **Description**: Control the t interval: dt = $\frac{\pi}{\omega_{cut}\times\omega enlarge}$
- **Default**: 10

#### cond_fwhm

- **Type**: Integer
- **Description**: We use gaussian functions to approxiamte $\delta(E)\approx \frac{1}{\sqrt{2\pi}\Delta E}e^{-\frac{E^2}{2{\Delta E}^2}}$. FWHM for conductivities, $FWHM=2*\sqrt{2\ln2}\cdot \Delta E$. The unit is eV.
- **Default**: 0.3

### Implicit solvation model

This part of variables are used to control the usage of implicit solvation model. This approach treats the solvent as a continuous medium instead of individual “explicit” solvent molecules, which means that the solute embedded in an implicit solvent and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath.

#### imp_sol

- **Type**: Boolean
- **Description**: If set to 1, an implicit solvation correction is considered.
- **Default**: 0

#### eb_k

- **Type**: Real
- **Description**: The relative permittivity of the bulk solvent, 80 for water. Used only if `imp_sol` == true.
- **Default**: 80

#### tau

- **Type**: Real
- **Description**: The effective surface tension parameter, which describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent that are not captured by the electrostatic terms. The unit is $Ry/Bohr^{2}$.
- **Default**: 1.0798e-05

#### sigma_k

- **Type**: Real
- **Description**: We assume a diffuse cavity that is implicitly determined by the electronic structure of the solute.
`sigma_k` is the parameter that describes the width of the diffuse cavity.
- **Default**: 0.6

#### nc_k

- **Type**: Real
- **Description**: It determines at what value of the electron density the dielectric cavity forms. 
The unit is $Bohr^{-3}$.
- **Default**: 0.00037
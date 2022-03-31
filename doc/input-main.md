# INPUT file
- [Structure of the file](#structure-of-the-file)
- [List of keywords](#list-of-keywords)
    - [System variables](#system-variables)

        [suffix](#suffix) | [ntype](#ntype) | [calculation](#calculation) | [dft_functional](#dft-functional) | [pseudo_type](#pseudo-type) | [kpar](#kpar) | [symmetry](#symmetry) | [pseudo_rcut](#pseudo-rcut) | [pseudo_mesh](#pseudo_mesh) | [nelec](#nelec) | [tot_magnetization](#tot-magnetization) | [mem_saver](#mem-saver) | [latname](#latname) | [init_wfc](#init_wfc) | [pw_seed](#pw_seed) | [init_chg](#init_chg) | [init_vel](#init_vel) | [diago_proc](#diago_proc) | [nbspline](#nbspline)

    - [Variables related to input files](#variables-related-to-input-files)

        [stru_file](#stru_file) | [kpoint_file](#kpoint-file) | [pseudo_dir](#pseudo-dir) | [orbital_dir](#orbital-dir) | [read_file_dir](#read-file-dir)

    - [Plane wave related variables](#plane-wave-related-variables)
    
        [ecutwfc](#ecutwfc) | [nx,ny,nz](#nx) | [diag_thr_e](#diag_thr_e) | [diag_cg_nmax](#diago-cg-maxiter) | [diago_david_ndim](#diago-david-ndim)

    - [Numerical atomic orbitals related variables](#numerical-atomic-orbitals-related-variables)

        [nb2d](#nb2d) | [lmaxmax](#lmaxmax) | [lcao_ecut](lcao-ecut) | [lcao_dk](#lcao-dk) | [lcao_dr](#lcao-dr) | [lcao_rmax](#lcao-rmax) | [search_radius](#search-radius) | [search_pbc](#search-pbc)

    - [Electronic structure](#electronic-structure)
    
        [basis_type](#basis-type) | [ks_solver](#ks-solver) | [nbands](#nbands) | [nbands_istate](#nbands-istate) | [nspin](#nspin) | [occupations](#occupations) | [smearing_method](#smearing_method) | [smearing_sigma](#smearing_sigma) | [mixing_type](#mixing-type) | [mixing_beta](#mixing-beta) | [mixing_ndim](#mixing-ndim) | [mixing_gg0](#mixing-gg0) | [gamma_only](#gamma-only) | [printe](#printe) | [scf_nmax](#scf_nmax) | [scf_thr_rho](#scf_thr_rho) | [chg_extrap](#chg_extrap) | [ocp](#ocp) | [ocp_set](#ocp_set)

    - [Geometry relaxation](#geometry-relaxation)
    
        [relax_nmax](#relax_nmax) | [cal_force](#cal_force) | [force_thr](#force-thr) | [force_thr_ev](#force-thr-ev) | [out_force](#out_force) | [relax_bfgs_w1](#bfgs-w1) | [relax_bfgs_w2](#bfgs-w2) | [relax_bfgs_rmax](#relax_bfgs_rmax) | [relax_bfgs_rmin](#relax_bfgs_rmin) | [relax_bfgs_init](#relax_bfgs_init) | [cal_stress](#cal_stress) | [stress_thr](#stress-thr) | [press1, press2, press3](#press) | [fixed_axes](#fixed-axes) | [relax_method](#relax_method) | [relax_cg_thr_e](#relax_cg_thr_e) | [cell_factor](#cell-factor)

    - [Variables related to program output](#variables-related-to-program-output)

        [out_mul](#out_mul) | [out_chg](#out_chg) | [out_pot](#out_pot) | [out_dm](#out-dm) | [out_wfc_pw](#out_wfc_pw) | [out_wfc_r](#out_wfc_r) | [out_wfc_lcao](#out_wfc_lcao) | [out_dos](#out-dos) | [out_band](#out-band) | [out_stru](#out-stru) | [out_level](#out_level) | [out_alllog](#out-alllog) | [out_mat_hs](#out_mat_hs) | [out_mat_r](#out_mat_r) | [out_mat_hs2](#out_mat_hs2) | [out_element_info](#out-element-info) | [restart_save](#restart_save) | [restart_load](#restart_load)

    - [Density of states](#density-of-states)

        [dos_edelta_ev](#dos-edelta-ev) | [dos_sigma](#dos-sigma) | [dos_scale](#dos-scale)
    
    - [Exact exchange](#exact-exchange) (under tests)
    
        [exx_hybrid_type](#exx-hybrid-type) | [exx_hybrid_alpha](#exx-hybrid-alpha) | [exx_hse_omega](#exx-hse-omega) | [exx_separate_loop](#exx-separate-loop) | [exx_hybrid_step](#exx-hybrid-step) | [exx_lambda](#exx-lambda) | [exx_pca_threshold](#exx-pca-threshold) | [exx_c_threshold](#exx-c-threshold) | [exx_v_threshold](#exx-v-threshold) | [exx_dm_threshold](#exx-dm-threshold) | [exx_schwarz_threshold](#exx-schwarz-threshold) | [exx_cauchy_threshold](#exx-cauchy-threshold) | [exx_ccp_threshold](#exx-ccp-threshold) | [exx_ccp_rmesh_times](#exx-ccp-rmesh-times) | [exx_distribute_type](#exx-distribute-type) | [exx_opt_orb_lmax](#exx-opt-orb-lmax) | [exx_opt_orb_ecut](#exx-opt-orb-ecut) | [exx_opt_orb_tolerence](#exx-opt-orb-tolerence)

    - [Molecular dynamics](#molecular-dynamics)

        [md_type](#md-type) | [md_nstep](#md_nstep) | [md_ensolver](#md-ensolver) | [md_restart](#md-restart) | [md_dt](#md-dt) | [md_t](#md-t) | [md_dumpfreq](#md-dumpfreq) | [md_restartfreq](#md-restartfreq) | [md_tfreq](#md-tfreq) | [md_mnhc](#md-mnhc) | [lj_rcut](#lj-rcut) | [lj_epsilon](#lj-epsilon) | [lj_sigma](#lj-sigma) | [msst_direction](#msst-direction) | [msst_vel](#msst-vel) | [msst_vis](#msst-vis) | [msst_tscale](#msst-tscale) | [msst_qmass](#msst-qmass) | [md_damp](#md-damp)

    - [DFT+U correction](#DFT_U-correction)

        [dft_plus_u](#dft_plus_u) | [orbital_corr](#orbital_corr) | [hubbard_u](#hubbard_u) | [hund_j](#hund_j) | [yukawa_potential](#yukawa_potential) | [omc](#omc)

    - [VdW correction](#vdw-correction)

        [vdw_method](#vdw-method) | [vdw_s6](#vdw-s6) | [vdw_s8](#vdw-s8) | [vdw_a1](#vdw-a1) | [vdw_a2](#vdw-a2) | [vdw_d](#vdw-d) | [vdw_abc](#vdw-abc) | [vdw_C6_file](#vdw-C6-file) | [vdw_C6_unit](#vdw-C6-unit) | [vdw_R0_file](#vdw-R0-file) | [vdw_R0_unit](#vdw-R0-unit) | [vdw_model](#vdw-model) | [vdw_radius](#vdw-radius) | [vdw_radius_unit](#vdw-radius-unit) | [vdw_cn_radius](#vdw-cn-radius) | [vdw_cn_radius_unit](#vdw-cn-radius-unit) | [vdw_period](#vdw-period)
        
    - [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)
    
        [berry_phase](#berry-phase) | [gdir](#gdir) | [towannier90](#towannier90) | [nnkpfile](#nnkpfile) | [wannier_spin](#wannier-spin)

    - [TDDFT: time dependent density functional theory](#TDDFT-doc)
    
        [tddft](#tddft) | [td_scf_thr_rho](#td_scf_thr_rho) | [td_dt](#td_dt) | [td_force_dt](#td_force_dt) | [td_vext](#td_vext) | [td_vext_dire](#td_vext_dire) | [td_timescale](#td_timescale) | [td_vexttype](#td_vexttype) | [td_vextout](#td_vextout) | [td_dipoleout](#td_dipoleout)

    - [Variables useful for debugging](#variables-useful-for-debugging)

        [nurse](#nurse) | [t_in_h](#t-in-h) | [vl_in_h](#vl-in-h) | [vnl_in_h](#vnl-in-h) | [test_force](#test-force) | [test_stress](#test-stress) | [colour](#colour) | [test_just_neighbor](#test-just-neighbor)
    - [DeePKS](#deepks)
    
        [deepks_out_labels](#out-descriptor) | [deepks_descriptor_lmax](#lmax-descriptor) | [deepks_scf](#deepks-scf) | [deepks_model](#model-file)

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

*Note: if a parameter name is not recognized by the program, the program will stop with an
error message.*

In the above example, the meanings of the parameters are:
- ntype : how many types of elements in the unit cell
- nbands : the number of bands to be calculated
- ecutwfc : the plane-wave energy cutoff for the wave function expansion (UNIT: Rydberg)

[back to top](#input-file)

## List of keywords
    
### System variables
This part of variables are used to control general system parameters.
- suffix<a id="suffix"></a>
    - *Type*: String
    - *Description*: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the ‘suffix’ is the name you can pick up for your convenience.
    - *Default*: ABACUS

    [back to top](#input-file)
- ntype<a id="ntype"></a>
    - *Type*: Integer
    - *Description*: Number of different atom species in this calculations. This value must be set. If the number you set is smaller than the atom species in the STRU file, ABACUS only read the ‘wrong number’ of atom information. If the number is larger than the atom species in the STRU file, ABACUS may stop and quit.
    - *Default*: ***No default value***

    [back to top](#input-file)

- calculation<a id="calculation"></a>
    - *Type*: String
    - *Description*: Specify the type of calculation.
        - *scf*: do self-consistent electronic structure calculation
        - *relax*: do structure relaxation calculation, one can ues ‘relax_nmax’ to decide how many ionic relaxations you want.
        - *cell-relax*: do cell relaxation calculation.
        - *nscf*: do the non self-consistent electronic structure calculations. For this option, you need a charge density file. For nscf calculations with planewave basis set, diag_thr_e should be <= 1d-3.
        - *istate*: Please see the explanation for variable `nbands_istate`.
        - *ienvelope*: Please see the explanation for variable `nbands_istate`.
        - *md*: molecular dynamics
    
        <mark>Note: *istate* and *ienvelope* only work for LCAO basis set and are not working right now.</mark> 
    - *Default*: scf
    
    [back to top](#input-file)

- dft_functional<a id="dft-functional"></a>
    - *Type*: String
    - *Description*: In our package, the XC functional can either be set explicitly using the dft_functional keyword as explained below, or set implicitly according to the XC functional information read from pseudopotential file. The user should ensure that the XC functional set in the INPUT file and the pseudopotential file are consistent. If more than one element is present in the system, make sure all of pseudopotentials have the same XC functional. **Currently only LDA and GGA are supported.**

        To be specific, we briefly explain the format of the pseudopotential file and the key information it contains. There are a few lines in Si’s GGA pseudopotential file Si_ONCV_PBE-1.0.upf:
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

        Possible values of this variable are:
        - none: the functional is specified implicity by the input pseudopotential file
        - lda: Perdew-Zunger local density approximation
        - pbe: Perdew-Burke-Ernzerhof general gradient approximation
    
        If the functional specified by the user is not consistent with the pseudopotential file, the program will stop with an error message.
    - *Default*: none

    [back to top](#input-file)

- pseudo_type<a id="pseudo-type"></a>
    - *Type*: String
    - *Description*: the format of pseudopotential files. Accepted value s are:
        - upf : .UPF format
        - vwr : .vwr format
        - upf201 : the new UPF format
        - blps : bulk local pseudopotential
    - *Default* : upf

    [back to top](#input-file)

- kpar<a id="kpar"></a>
    - *Type*: Integer
    - *Description*: devide all processors into kpar groups, and k points will be distributed among each group. The value taken should be less than or equal to the number of k points as well as the number of MPI threads.
    - *Default*: 1

    [back to top](#input-file)

- symmetry<a id="symmetry"></a>
    - *Type*: Integer
    - *Description*: takes value 0 and 1, if set to 1, symmetry analysis will be performed to determine the type of Bravais lattice and associated symmetry operations.
    - *Default*: 0

    [back to top](#input-file)

- pseudo_rcut<a id="pseudo-rcut"></a>
    - *Type*: Real
    - *Description*: Cut-off of radial integration for pseudopotentials, in Bohr.
    - *Default*: 15

    [back to top](#input-file)

- pseudo_mesh<a id="pseudo_mesh"></a>
    - *Type*: Integer
    - *Description*: If set to 0, then use our own mesh for radial integration of pseudopotentials; if set to 1, then use the mesh that is consistent with quantum espresso.
    - *Default*: 0

    [back to top](#input-file)

- nelec<a id="nelec"></a>
    - *Type*: Real
    - *Description*: If >0.0, this denotes total number of electrons in the system. Must be less than 2*nbands. If set to 0.0, the total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
    - *Default*: 0.0

    [back to top](#input-file)

- tot_magnetization<a id="tot-magnetization"></a>
    - *Type*: Real
    - *Description*: Total magnetization of the system.
    - *Default*: 0.0

    [back to top](#input-file)

- mem_saver<a id="mem-saver"></a>
    - *Type*: Boolean
    - *Description*: Used only for nscf calculations. If set to 1, then a memory saving technique will be used for many k point calculations.
    - *Default*: 0

    [back to top](#input-file)

- latname<a id="latname"></a>
    - *Type*: String
    - *Description*: Specifies the type of Bravias lattice. When set to "test", the three lattice vectors are supplied explicitly in STRU file. When set to certain Bravais lattice type, there is no need to provide lattice vector, but a few lattice parameters might be required. For more information regarding this parameter, consult the [page on STRU file](input-stru.md).
    Available options are:
        - "test": free strcture.
        - "sc": simple cubie.
        - "fcc": face-centered cubic.
        - "bcc": body-centered cubic.
        - "hexagonal": hexagonal.
        - "trigonal": trigonal.
        - "st": simple tetragonal.
        - "bct": body-centered tetragonal.
        - "so": orthorhombic.
        - "baco": base-centered orthorhombic.
        - "fco": face-centered orthorhombic.
        - "bco": body-centered orthorhombic.
        - "sm": simple monoclinic.
        - "bacm": base-centered monoclinic.
        - "triclinic": triclinic.
    - *Default*: "test"

    [back to top](#input-file)

- init_wfc<a id="init_wfc"></a>
    - *Type*: String
    - *Description*: Only useful for plane wave basis only now. It is the name of the starting wave functions. In the future we should also make this         variable available for localized orbitals set. 
    Available options are:
        - "atomic": from atomic pseudo wave functions. If they are not enough, other wave functions are initialized with random numbers.
        - "atomic+random": add small random numbers on atomic pseudo-wavefunctions
        - "file": from file
        - "random": random numbers
    - *Default*:"atomic"

    [back to top](#input-file)

- pw_seed<a id="pw_seed"></a>
    - *Type*: Integer
    - *Description*: Only useful for plane wave basis only now. It is the random seed to initialize wave functions. Only positive integers are avilable.
    - *Default*:0

    [back to top](#input-file)

- init_chg<a id="init_chg"></a>
    - *Type*: String
    - *Description*: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density. If set this to ‘atomic’, the density is starting from summation of atomic density of single atoms. If set this to ‘file’, the density will be read in from file. Besides, when you do ‘nspin=1’ calculation, you only need the density file SPIN1_CHGCAR. However, if you do ‘nspin=2’ calculation, you also need the density file SPIN2_CHGCAR. The density file should be output with these names if you set out_chg = 1 in INPUT file.
    - *Default*:atomic

    [back to top](#input-file)

- init_vel<a id="init_vel"></a>
    - *Type*: Boolean
    - *Description*: Read the atom velocity from the atom file (STRU) if set to true.
    - *Default*: false

    [back to top](#input-file)

- diago_proc<a id="diago_proc"></a>
    - *Type*: Integer
    - *Descrption*: If set to a positive number, then it specifies the number of threads used for carrying out diagonalization. Must be less than or equal to total number of MPI threads. Also, when cg diagonalization is used, diago_proc must be same as total number of MPI threads. If set to 0, then it will be set to the number of MPI threads. Normally, it is fine just leaving it to default value. Only used for pw base.
    - *Default*: 0

    [back to top](#input-file)

- nbspline<a id="nbspline"></a>
    - *Type*: Integer
    - *Descrption*: If set to a natural number, a Cardinal B-spline interpolation will be used to calculate Structure Factor. "nbspline" represents the order of B-spline basis and larger one can get more accurate results but cost more.
    It is turned off by default.
    - *Default*: -1

    [back to top](#input-file)

### Variables related to input files
This part of variables are used to control input files related parameters.

- stru_file<a id="stru_file"></a>
    - *Type*: String
    - *Description*: This parameter specifies the name of structure file which contains various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
    - *Default*: STRU

    [back to top](#input-file)

- kpoint_file<a id="kpoint-file"></a>
    - *Type*: String
    - *Description*: This parameter specifies the name of k-points file. Note that if you use atomic orbitals as basis, and you only use gamma point, you don’t need to have k-point file in your directory, ABACUS will automatically generate ‘KPT’ file. Otherwise, if you use more than one k-point, please do remember the algorithm in ABACUS is different for gamma only and various k-point dependent simulations. So first you should turn off the k-point algorithm by set `gamma_only = 0` in `INPUT` and then you should setup your own k-points file.
    - *Default*: KPT

    [back to top](#input-file)

- pseudo_dir<a id="pseudo-dir"></a>
    - *Type*: String
    - *Description*: This parameter specifies pseudopotential directory.
    - *Default*: ./

    [back to top](#input-file)

- orbital_dir<a id="orbital-dir"></a>
    - *Type*: String
    - *Description*: This parameter specifies orbital file directory.
    - *Default*: ./

    [back to top](#input-file)

- read_file_dir<a id="read-file-dir"></a>
    - *Type*: String
    - *Description*: when the program needs to read files such as electron density(`SPIN1_CHG`) as a starting point, this variables tells the location of the files. For example, './' means the file is located in the working directory.
    - *Default*: OUT.$suffix

    [back to top](#input-file)

### Plane wave related variables
This part of variables are used to control the plane wave related parameters.

- ecutwfc<a id="ecutwfc"></a>
    - *Type*: Real
    - *Description*: Energy cutoff for plane wave functions, the unit is **Rydberg**. Note that even for localized orbitals basis, you still need to setup a energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set, etc. Also, because our orbitals are generated by matching localized orbitals to a chosen set of wave functions from certain energy cutoff, so this set of localize orbitals are most accurate under this same plane wave energy cutoff.
    - *Default*: 50

    [back to top](#input-file)

- nx, ny, nz<a id="nx"></a>
    - *Type*: Integer
    - *Description*: If set to a positive number, then the three variables specify the numbers of FFT grid points in x, y, z directions, respectively. If set to 0, the number will be calculated from ecutwfc.
    - *Default*: 0

    [back to top](#input-file)

- diag_thr_e<a id="diag_thr_e"></a>
    - *Type*: Real
    - *Description*: Only used when you use diago_type = cg or diago_type = david. It indicates the threshold for the first electronic iteration, from the second iteration the diag_thr_e will be updated automatically. **For nscf calculations with planewave basis set, diag_thr_e should be <= 1d-3.**
    - *Default*: 0.01

    [back to top](#input-file)

- diag_cg_nmax<a id="diago-cg-maxiter"></a>
    - *Type*: Integer
    - *Description*: Only useful when you use ks_solver = cg or ks_solver = dav. It indicates the maximal iteration number for cg/david method.
    - *Default*: 40

    [back to top](#input-file)

- diago_david_ndim<a id="diago-david-ndim"></a>
    - *Type*: Integer
    - *Description*: Only useful when you use ks_solver = dav. It indicates the maximal dimension for the Davidson method.
    - *Default*: 10

    [back to top](#input-file)

### Numerical atomic orbitals related variables
This part of variables are used to control the numerical atomic orbitals related parameters.

- nb2d<a id="nb2d"></a>
    - *Type*: Integer
    - *Description*: In LCAO calculations, we arrange the total number of processors in an 2D array, so that we can partition the wavefunction matrix (number of bands*total size of atomic orbital basis) and distribute them in this 2D array. When the system is large, we group processors into sizes of nb2d, so that multiple processors take care of one row block (a group of atomic orbitals) in the wavefunction matrix. If set to 0, nb2d will be automatically set in the program according to the size of atomic orbital basis: 
        - if size <= 500 : nb2d = 1
        - if 500 < size <= 1000 : nb2d = 32
        - if size > 1000 : nb2d = 64;
    - *Default*: 0

    [back to top](#input-file)

- lmaxmax<a id="lmaxmax"></a>
    - *Type*: Integer
    - *Description*: If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.
    - *Default*: 2.

    [back to top](#input-file)

- lcao_ecut<a id="lcao-ecut"></a>

    - *Type*: Real
    - *Description*: Energy cutoff when calculating LCAO two-center integrals. In Ry.
    - *Default*: 50

    [back to top](#input-file)
    
- lcao_dk<a id="lcao-dk"></a>

    - *Type*: Real
    - *Description*: Delta k for 1D integration in LCAO
    - *Default*: 0.01

    [back to top](#input-file)

- lcao_dr<a id="lcao-dr"></a>

    - *Type*: Real
    - *Description*: Delta r for 1D integration in LCAO
    - *Default*: 0.01

    [back to top](#input-file)

- lcao_rmax<a id="lcao-rmax"></a>

    - *Type*: Real
    - *Description*: Max R for 1D two-center integration table
    - *Default*: 30

    [back to top](#input-file)

- search_radius<a id="search-radius"></a>

    - *Type*: Real
    - *Description*: Set the search radius for finding neighbouring atoms. If set to -1, then the radius will be set to maximum of projector and orbital cut-off.
    - *Default*: -1

    [back to top](#input-file)

- search_pbc<a id="search-pbc"></a>

    - *Type*: Boolean
    - *Description*: In searching for neighbouring atoms, if set to 1, then periodic images will also be searched. If set to 0, then periodic images will not be searched.
    - *Default*: 1

    [back to top](#input-file)

### Electronic structure
This part of variables are used to control the electronic structure and geometry relaxation
calculations.
- basis_type<a id="basis-type"></a>
    - *Type*: String
    - *Description*: This is very important parameters to choose basis set in ABACUS.
        - *pw*: Using plane-wave basis set only.
        - *lcao_in_pw*: Expand the localized atomic set in plane-wave basis.
        - lcao: Using localized atomic orbital sets.
    - *Default*: pw

    [back to top](#input-file)
- ks_solver<a id="ks-solver"></a>
    - *Type*: String
    - *Description*: It’s about choice of diagonalization methods for hamiltonian matrix expanded in a certain basis set.

        For plane-wave basis,
        - cg: cg method.
        - dav: the Davidson algorithm. (Currently not working with Intel MKL library).

        For atomic orbitals basis,
        - genelpa: This method should be used if you choose localized orbitals.
        - hpseps: old method, still used.
        - lapack: lapack can be used for localized orbitals, but is only used for single processor.

        If you set ks_solver=‘hpseps’ for basis_type=‘pw’, the program will be stopped with an error message:

        ```
        hpseps can not be used with plane wave basis.
        ```

        Then the user has to correct the input file and restart the calculation.
    - *Default*: cg (pw) or genelpa (lcao)

    [back to top](#input-file)

- nbands<a id="nbands"></a>
    - *Type*: Integer
    - *Description*: Number of bands to calculate. It is recommended you setup this value, especially when you use smearing techniques, more bands should be included.
    - *Default (nspin=1)*: *1.2\*occupied_bands, occupied_bands+10)*
    - *Default (nspin=2)*: *max(1.2\*nelec, nelec+20)*

    [back to top](#input-file)
    
- nbands_istate<a id="nbands-istate"></a>
    - *Type*: Integer
    - *Description*: Only used when `calculation = ienvelope` or `calculation = istate`, this variable indicates how many bands around Fermi level you would like to calculate. `ienvelope` means to calculate the envelope functions of wave functions <em>&Psi;<sub>i</sub>=&Sigma;<sub>&mu;</sub>C<sub>i&mu;</sub>&Phi;<sub>&mu;</sub></em>, where <em>&Psi;<sub>i</sub></em> is the ith wave function with the band index <em>i</em> and <em></sub>&Phi;<sub>&mu;</sub></em> is the localized atomic orbital set. `istate` means to calculate the density of each wave function <em>|&Psi;<sub>i</sub>|<sup>2</sup></em>. Specifically, suppose we have highest occupied bands at 100th wave functions. And if you set this variable to 5, it will print five wave functions from 96th to 105th. But before all this can be carried out, the wave functions coefficients  should be first calculated and written into a file by setting the flag `out_wfc_lcao = 1`.
    -   *Default*: 5

    [back to top](#input-file)

- nspin<a id="nspin"></a>
    - *Type*: Integer
    - *Description*: Number of spin components of wave functions. There are only two choices now: 1 or 2, meaning non spin or collinear spin.
    - *Default*: 1

    [back to top](#input-file)

- occupations<a id="occupations"></a>
    - *Type*: String
    - *Description*: Specifies how to calculate the occupations of bands. Available options are:
        - 'smearing' : gaussian smearing for metals; see also variables `smearing_method` and `smearing_sigma`.
        - 'tetrahedra' : Tetrahedron method, Bloechl's version: [P.E. Bloechl, PRB 49, 16223 (1994)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223). Requires a uniform grid of k-points that are automatically generated. Well suited for calculation of DOS, less so (because not variational) for force/optimization/dynamics calculations.
        - 'fixed' : for insulators with a gap
    - *Default*: 'smearing'

    [back to top](#input-file)

- smearing_method<a id="smearing_method"></a>
    - *Type*: String
    - *Description*: It indicates which occupation and smearing method is used in the calculation.
        - fixed: use fixed occupations.
        - gauss or gaussian: use gaussian smearing method.
        - mp: use methfessel-paxton smearing method. The method recommends for metals. 
    - *Default*: fixed

    [back to top](#input-file)

- smearing_sigma<a id="smearing_sigma"></a>
    - *Type*: Real
    - *Description*: energy range for smearing, the unit is Rydberg.
    - *Default*: 0.001

    [back to top](#input-file)
    
- mixing_type<a id="mixing-type"></a>
    - *Type*: String
    - *Description*: Charge mixing methods.
        - plain: Just simple mixing.
        - kerker: Use kerker method, which is the mixing method in G space.
        - pulay: Standard Pulay method.
        - pulay-kerker:
    - *Default*: pulay

    [back to top](#input-file)

- mixing_beta<a id="mixing-beta"></a>
    - *Type*: Real
    - *Description*: mixing parameter: 0 means no new charge
    - *Default*: 0.7

    [back to top](#input-file)

- mixing_ndim<a id="mixing-ndim"></a>
    - *Type*: Integer
    - *Description*: It indicates the mixing dimensions in Pulay, Pulay method use the density from previous mixing_ndim steps and do a charge mixing based on these density.
    - *Default*: 8

    [back to top](#input-file)

- mixing_gg0<a id="mixing-gg0"></a>
    - *Type*: Real
    - *Description*: used in pulay-kerker mixing method
    - *Default*: 1.5

    [back to top](#input-file)

- gamma_only<a id="gamma-only"></a>
    - *Type*: Integer
    - *Description*: It is an important parameter **only to be used in localized orbitals set**.
    It you set gamma_only = 1, ABACUS use gamma only, the algorithm is fast and you don’t need to specify the k-points file. If you set gamma_only = 0, more than one k-point is used and the ABACUS is slower compared to gamma only algorithm.
    - *Default*: 0

    [back to top](#input-file)

- printe<a id="printe"></a>
    - *Type*: Integer
    - *Description*: Print out energy for each band for every printe steps
    - *Default*: 100

    [back to top](#input-file)

- scf_nmax<a id="scf_nmax"></a>
    - *Type*: Integer
    - *Description*:This variable indicates the maximal iteration number for electronic iterations.
    - *Default*: 40

    [back to top](#input-file)

- scf_thr_rho<a id="scf_thr_rho"></a>
    - *Type*: Real
    - *Description*: An important parameter in ABACUS. It’s the threshold for electronic iteration. It represents the charge density error between two sequential density from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
    - *Default*:1e-06

    [back to top](#input-file)

- chg_extrap<a id="chg_extrap"></a>
    - *Type*: String
    - *Description*: Methods to do extrapolation of density when ABACUS is doing geometry relaxations.
        - atomic: atomic extrapolation
        - first-order: first-order extrapolation
        - second-order: second-order extrapolation
    - *Default*:atomic

    [back to top](#input-file)


- ocp<a id="ocp"></a>
    - *Type*: Boolean
    - *Description*: option for choose whether calcualting constrained DFT or not.
    Only used for TDDFT.
    - *Default*:0

    [back to top](#input-file)

- ocp_set<a id="ocp_set"></a>
    - *Type*: string
    - *Description*: If ocp is true, the ocp_set is a string to set the number of occupancy, like 1 10 * 1 0 1 representing the 13 band occupancy, 12th band occupancy 0 and the rest 1, the code is parsing this string into an array through a regular expression.
    - *Default*:none

    [back to top](#input-file)

### Geometry relaxation
This part of variables are used to control the geometry relaxation.

- relax_nmax<a id="relax_nmax"></a>
    - *Type*: Integer
    - *Description*: The maximal number of ionic iteration steps, the minimal value is 1.
    - *Default*: 1

    [back to top](#input-file)

- cal_force<a id="cal_force"></a>
    - *Description*: If set to 1, calculate the force at the end of the electronic iteration. 0 means the force calculation is turned off.
    - *Default*: 0

    [back to top](#input-file)

- force_thr<a id="force-thr"></a>
    - *Type*: Real
    - *Description*: The threshold of the force convergence, it indicates the largest force among all the atoms, the unit is Ry=Bohr
    - *Default*: 0.000388935 Ry/Bohr = 0.01 eV/Angstrom

    [back to top](#input-file)

- force_thr_ev<a id="force-thr-ev"></a>
    - *Type*: Real
    - *Description*: The threshold of the force convergence, has the same function as force_thr, just the unit is different, it is eV=Angstrom, you can choose either one as you like. The recommendation value for using atomic orbitals is 0:04 eV/Angstrom.
    - *Default*: 0.01 eV/Angstrom

    [back to top](#input-file)

- out_force<a id="out_force"></a>
    - *Type*: Integer
    - *Description*: Determines whether to output the out_force into a file named `Force.dat` or not. If 1, then force will be written; if 0, then the force will not be written.
    - *Default*: 0

    [back to top](#input-file)

- relax_bfgs_w1<a id="bfgs-w1"></a>
    - *Type*: Real
    - *Description*: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
    - *Default*: 0.01

    [back to top](#input-file)

- relax_bfgs_w2<a id="bfgs-w2"></a>
    - *Type*: Real
    - *Description*: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
    - *Default*: 0.5

    [back to top](#input-file)

- relax_bfgs_rmax<a id="relax_bfgs_rmax"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the maximal movement of all the atoms. The sum of the movements from all atoms can be increased during the optimization steps. However, it will not be larger than relax_bfgs_rmax Bohr.
    - *Default*: 0.8

    [back to top](#input-file)

- relax_bfgs_rmin<a id="relax_bfgs_rmin"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the minimal movement of all the atoms. When the movement of all the atoms is smaller than relax_bfgs_rmin Bohr , and the force convergence is still not achieved, the calculation will break down.
    - *Default*: 1e-5

    [back to top](#input-file)

- relax_bfgs_init<a id="relax_bfgs_init"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the initial movement of all the atoms. The sum of the movements from all atoms is relax_bfgs_init Bohr.
    - *Default*: 0.5

    [back to top](#input-file)

- cal_stress<a id="cal_stress"></a>
    - *Type*: Integer
    - *Description*: If set to 1, calculate the stress at the end of the electronic iteration. 0 means the stress calculation is turned off.
    - *Default*: 0

    [back to top](#input-file)

- stress_thr<a id="stress-thr"></a>
    - *Type*: Real
    - *Description*: The threshold of the stress convergence, it indicates the largest stress among all the directions, the unit is KBar,
    - *Default*: 0.01

    [back to top](#input-file)

- press1, press2, press3<a id="press"></a>
    - *Type*: Real
    - *Description*: the external pressures along three axes,the compressive stress is taken to be positive, the unit is KBar.
    - *Default*: 0

    [back to top](#input-file)

- fixed_axes<a id="fixed-axes"></a>
    - *Type*: String
    - *Description*:which axes are fixed when do cell relaxation. Possible choices are:
        - None : default; all can relax
        - volume : relaxation with fixed volume
        - a : fix a axis during relaxation
        - b : fix b axis during relaxation
        - c : fix c axis during relaxation
        - ab : fix both a and b axes during relaxation
        - ac : fix both a and c axes during relaxation
        - bc : fix both b and c axes during relaxation
        - abc : fix all three axes during relaxation
    - *Default*: None

    [back to top](#input-file)

- relax_method<a id="relax_method"></a>
    - *Type*: String
    - *Description*: The method to do geometry optimizations. If set to bfgs, using BFGS algorithm. If set to cg, using cg algorithm. If set to sd, using steepest-descent lgorithm.
    - *Default*: cg

    [back to top](#input-file)

- relax_cg_thr_e<a id="relax_cg_thr_e"></a>
    - *Type*: Real
    - *Description*: When move-method is set to 'cg-bfgs', a mixed cg-bfgs algorithm is used. The ions first move according to cg method, then switched to bfgs when maximum of force on atoms is reduced below cg-threshold. Unit is eV/Angstrom.
    - *Default*: 0.5

    [back to top](#input-file)

- cell_factor<a id="cell-factor"></a>
    - *Type*: Real
    - *Description*: Used in the construction of the pseudopotential tables. It should exceed the maximum linear contraction of the cell during a simulation.
    - *Default*: 1.2

    [back to top](#input-file)

### Variables related to program output
This part of variables are used to control the output of properties.

- out_mul<a id="out_mul"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the Mulliken population analysis result. The name of the output file is mulliken.txt
    - *Default*: 0

    [back to top](#input-file)

- out_chg<a id="out_chg"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the charge density on real space grid. The name of the density file is SPIN1_CHGCAR and SPIN2_CHGCAR (if nspin = 2). Suppose each density on grid has coordinate (x; y; z). The circle order of the density on real space grid is: z is the outer loop, then y and finally x (x is moving fastest).
    - *Default*: 0

    [back to top](#input-file)
- out_pot<a id="out_pot"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the local potential on real space grid. The name of the file is SPIN1_POT and SPIN2_POT (if nspin = 2). If set to 2, ABACUS will output the electrostatic potential on real space grid. The name of the file is ElecStaticPot and ElecStaticP ot_AV E (along the z-axis).
    - *Default*: 0

    [back to top](#input-file)
- out_dm<a id="out-dm"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the density matrix of localized orbitals, only useful for localized orbitals set. The name of the output file is SPIN1_DM and SPIN2_DM in the output directory.
    - *Default*: 0

    [back to top](#input-file)
- out_wfc_pw<a id="out_wfc_pw"></a>
    - *Type*: Integer
    - *Description*: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs the coefficients of wave functions into text files. The file names are WAVEFUNC$K.txt, where $K is the index of k point. When set this variable to 2, results are stored in binary files. The file names are WAVEFUNC$K.dat.
    - *Default*: 0

    [back to top](#input-file)
- out_wfc_r<a id="out_wfc_r"></a>
    - *Type*: Integer
    - *Description*: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs real-space wave functions into  `OUT.suffix/wfc_realspace/`. The file names are wfc_realspace$K$B, where $K is the index of k point, $B is the index of band. 
    - *Default*: 0

    [back to top](#input-file)
- out_wfc_lcao<a id="out_wfc_lcao"></a>
    - *Type*: Integer
    - *Description*: **Only used in localized orbitals set**. If set to 1, ABACUS will output the wave functions coefficients.
    - *Default*: 0

    [back to top](#input-file)
- out_dos<a id="out-dos"></a>
    - *Type*: Integer
    - *Description*: Controls whether to output the density of state (DOS). For more information, refer to the [worked example](examples/dos.md).
    - *Default*: 0

    [back to top](#input-file)

- out_band<a id="out-band"></a>
    - *Type*: Integer
    - *Description*: Controls whether to output the band structure. For mroe information, refer to the [worked example](examples/band-struc.md)
    - *Default*: 0
    
    [back to top](#input-file)

- out_stru<a id="out-stru"></a>
    - *Type*: Boolean
    - *Description*: If set to 1, then tje structure files will be written after each ion step
    - *Default*: 0

    [back to top](#input-file)

- out_level<a id="out-level"></a>
    - *Type*: String
    - *Description*: Controls the level of output. "ie" means write output at electron level; "i" means write additional output at ions level.
    - *Default*: ie

    [back to top](#input-file)

- out_alllog<a id="out-alllog"></a>
    - *Type*: Integer
    - *Description*: determines whether to write log from all ranks in an MPI run. If set to be 1, then each rank will write detained running information to a file named running_${calculation}\_(${rank}+1).log. If set to 0, log will only be written from rank 0 into a file named running_${calculation}.log.
    - *Default*: 0

    [back to top](#input-file)

- out_mat_hs<a id="out_mat_hs"></a>
    - *Type*: Boolean
    - *Description*: Only for LCAO calculations. When set to 1, ABACUS will generate two files `data-H` and `data-S` that store the Hamiltonian and S matrix in k space, respectively.
    - *Default*: 0

    [back to top](#input-file)

- out_mat_r<a id="out_mat_r"></a>
    - *Type*: Boolean
    - *Description*: Only for LCAO and not gamma_only calculations. When set to 1, ABACUS will generate a file with name staring with `data-rR-tr` which stores overlap matrix as a function of R, in units of lattice vectors.
    - *Default*: 0

    [back to top](#input-file)

- out_mat_hs2<a id="out_mat_hs2"></a>
    - *Type*: Boolean
    - *Description*: Only for LCAO and not gamma_only calculations. When set to 1, ABACUS will generate two files starting with `data-HR-sparse` and `data-SR-sparse` that store the Hamiltonian and S matrix in real space, respectively, as functions of R, in units of lattice vectors.
    - *Default*: 0

    [back to top](#input-file)

- out_element_info<a id="out-element-info"></a>
    - *Type*: Boolean
    - *Description*: When set to 1, ABACUS will generate a new directory under OUT.suffix path named as element name such as 'Si', which contained files "Si-d1-orbital-dru.dat  Si-p2-orbital-k.dat    Si-s2-orbital-dru.dat
    Si-d1-orbital-k.dat    Si-p2-orbital-r.dat    Si-s2-orbital-k.dat
    Si-d1-orbital-r.dat    Si-p2-orbital-ru.dat   Si-s2-orbital-r.dat
    Si-d1-orbital-ru.dat   Si-p-proj-k.dat        Si-s2-orbital-ru.dat
    Si.NONLOCAL            Si-p-proj-r.dat        Si-s-proj-k.dat
    Si-p1-orbital-dru.dat  Si-p-proj-ru.dat       Si-s-proj-r.dat
    Si-p1-orbital-k.dat    Si-s1-orbital-dru.dat  Si-s-proj-ru.dat
    Si-p1-orbital-r.dat    Si-s1-orbital-k.dat    v_loc_g.dat
    Si-p1-orbital-ru.dat   Si-s1-orbital-r.dat
    Si-p2-orbital-dru.dat  Si-s1-orbital-ru.dat " for example.

    - *Default*: 0

    [back to top](#input-file)

- restart_save<a id="restart_save"></a>
    - *Type*: Boolean
    - *Description*: Only for LCAO, store charge density file and H matrix file every scf step for restart.
    - *Default*: 0

    [back to top](#input-file)

- restart_load<a id="restart_load"></a>
    - *Type*: Boolean
    - *Description*: Only for LCAO, used for restart, only if that:
        * set restart_save as true and do scf calculation before.
        * please ensure suffix is same with calculation before and density file and H matrix file is exist. 

      restart from stored density file and H matrix file.
    - *Default*: 0
    
    [back to top](#input-file)

### Density of states
This part of variables are used to control the calculation of DOS.

- dos_edelta_ev<a id="dos-edelta-ev"></a>
    - *Type*: Real
    - *Description*: controls the step size in writing DOS (in eV).
    - *Default*: 0.1

    [back to top](#input-file)

- dos_sigma<a id="dos-sigma"></a>
    - *Type*: Real
    - *Description*: controls the width of Gaussian factor when obtaining smeared DOS (in eV).
    - *Default*: 0.07

    [back to top](#input-file)

- dos_scale<a id="dos-scale"></a>
    - *Type*: Real
    - *Description*: the energy range of dos output is given by (emax-emin)*(1+dos_scale), centered at (emax+emin)/2.
    - *Default*: 0.01

    [back to top](#input-file)


### DeePKS
This part of variables are used to control the usage of DeePKS method (a comprehensive data-driven approach to improve accuracy of DFT).
Warning: this function is not robust enough for version 2.2.0. Please try these variables in https://github.com/deepmodeling/abacus-develop/tree/deepks .

- deepks_out_labels<a id="out-descriptor"></a>
    - *Type*: Boolean
    - *Description*: when set to 1, ABACUS will calculate and output descriptor for DeePKS training. In `LCAO` calculation, a path of *.orb file is needed to be specified under `NUMERICAL_DESCRIPTOR`in `STRU`file. For example: 
    ```
    NUMERICAL_ORBITAL
    H_gga_8au_60Ry_2s1p.orb
    O_gga_7au_60Ry_2s2p1d.orb

    NUMERICAL_DESCRIPTOR
    jle.orb
    ```
    - *Default*: 0

    [back to top](#input-file)
- deepks_descriptor_lmax<a id="lmax-descriptor"></a>
    - *Type*: Integer
    - *Description*: control the max angular momentum of descriptor basis. 
    - *Default*: 0

    [back to top](#input-file)
- deepks_scf<a id="deepks-scf"></a>
    - *Type*: Boolean
    - *Description*: only when deepks is enabled in `LCAO` calculation can this variable set to 1. Then, a trained, traced model file is needed for self-consistant field iteration in DeePKS method.
    - *Default*: 0

    [back to top](#input-file)
- deepks_model<a id="model-file"></a>
    - *Type*: String
    - *Description*: the path of the trained, traced NN model file (generated by deepks-kit). used when deepks_scf is set to 1.
    - *Default*: None

    [back to top](#input-file)

### Exact Exchange
This part of variables are relevant when using hybrid functionals

- exx_hybrid_type<a id="exx-hybrid-type"></a>
    - *Type*: String
    - *Description*: Type of hybrid functional used. Options are "hf" (pure Hartree-Fock), "pbe0"(PBE0), "hse" (Note: in order to use HSE functional, LIBXC is required). Note also that HSE has been tested while PBE0 has NOT been fully tested yet, and the maxmum parallel cpus for running exx is Nx(N+1)/2, with N being the number of atoms.

    
        If set to "no", then no hybrid functional is used (i.e.,Fock exchange is not included.)

        If set to "opt_orb", the program will not perform hybrid functional calculation. Instead, it is going to generate opt-ABFs as discussed in this [article](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c00481).
    - *Default*: "no"

    [back to top](#input-file)

- exx_hybrid_alpha<a id="exx-hybrid-alpha"></a>
    - *Type*: Real
    - *Description*: fraction of Fock exchange in hybrid functionals, so that E<sub>X</sub>=&alpha;F<sub>X</sub>+(1-&alpha;)E_<sub>X,LDA/GGA</sub>
    - *Default*: 0.25

    [back to top](#input-file)

- exx_hse_omega<a id="exx-hse-omega"></a>
    - *Type*: 
    - *Description*: range-separation parameter in HSE functional, such that 1/r=erfc(&omega;r)/r+erf(&omega;r)/r.
    - *Default*: 0.11

    [back to top](#input-file)
adial integration for pseudopotentials, in Bohr.
@@ -214,6 +279,13 @@ This part of variables are used to control general system para
    - *Default*: 0.3

    [back to top](#input-file)

- exx_pca_threshold<a id="exx-pca-threshold"></a>
    - *Type*: Real
    - *Description*: To accelerate the evaluation of four-center integrals (ik|jl), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): &phi;<sub>i</sub>&phi;<sub>j</sub>~C<sup>k</sup><sub>ij</sub>P<sub>k</sub>. The size of the ABF (i.e. number of P<sub>k</sub>) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculations becomes faster. However this comes at the cost of computational accuracy. A relatively safe choice of the value is 1d-4.
    - *Default*: 0

    [back to top](#input-file)

- exx_c_threshold<a id="exx-c-threshold"></a>
    - *Type*: Real
    - *Description*: See also the entry [exx_pca_threshold](#exx-pca-threshold). Smaller components (less than exx_c_threshold) of the C<sup>k</sup><sub>ij</sub> matrix is neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-4.
    - *Default*: 0

    [back to top](#input-file)

- exx_v_threshold<a id="exx-v-threshold"></a>
    - *Type*: Real
    - *Description*: See also the entry [exx_pca_threshold](#exx-pca-threshold). With the approximation &phi;<sub>i</sub>&phi;<sub>j</sub>~C<sup>k</sup><sub>ij</sub>P<sub>k</sub>, the four-center integral in Fock exchange is expressed as (ik|jl)=&Sigma;<sub>a,b</sub>C<sup>a</sup><sub>ij</sub>V<sub>ab</sub>C<sup>b</sup><sub>kl</sub>, where V<sub>ab</sub>=(P<sub>a</sub>|P<sub>b</sub>) is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
    - *Default*: 0

    [back to top](#input-file)

- exx_dm_threshold<a id="exx-dm-threshold"></a>
    - *Type*: Real
    - *Description*: The Fock exchange can be expressed as &Sigma;<sub>k,l</sub>(ik|jl)D<sub>kl</sub> where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-4.
    - *Default*: 0

    [back to top](#input-file)

- exx_schwarz_threshold<a id="exx-schwarz-threshold"></a>
    - *Type*: Real
    - *Description*: In practice the four-center integrals are sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each integral before carrying out explicit evaluations. Those that are smaller than exx_schwarz_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-5.
    - *Default*: 0

    [back to top](#input-file)

- exx_cauchy_threshold<a id="exx-cauchy-threshold"></a>
    - *Type*: Real
    - *Description*: In practice the Fock exchange matrix is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1d-7.
    - *Default*: 0

    [back to top](#input-file)

- exx_ccp_threshold<a id="exx-ccp-threshold"></a>
    - *Type*: Real
    - *Description*: It is related to the cutoff of on-site Coulomb potentials, currently not used.
    - *Default*: 1e-8

    [back to top](#input-file)

- exx_ccp_rmesh_times<a id="exx-ccp-rmesh-times"></a>
    - *Type*: Real
    - *Description*: This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. For HSE1, setting it to 1 is enough. But for PBE0, a much larger number must be used.
    - *Default*: 10

    [back to top](#input-file)

- exx_distribute_type<a id="exx-distribute-type"></a>
    - *Type*: String
    - *Description*: When running in parallel, the evaluation of Fock exchange is done by distributing atom pairs on different threads, then gather the results. exx_distribute_type governs the mechanism of distribution. Available options are "htime", "order", "kmean1" and "kmeans2". "order" is where atom pairs are simply distributed by their orders. "hmeans" is a distribution where the balance in time is achieved on each processor, hence if the memory is sufficient, this is the recommended method. "kmeans1" and "kmeans2" are two methods where the k-means clustering method is used to reduce memory requirement. They might be necessary for very large systems.
    - *Default*: "htime"

    [back to top](#input-file)

- exx_opt_orb_lmax<a id="exx-opt-orb-lmax"></a>
    - *Type*: Integer
    - *Description*: See also the entry [exx_hybrid_type](#exx-hybrid-type). This parameter is only relevant when exx_hybrid_type="opt_orb". The radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. exx_opt_orb_lmax gives the maximum l of the spherical Bessel functions. A reasonable choice is 2.
    - *Default*: 0

    [back to top](#input-file)

- exx_opt_orb_ecut<a id="exx-opt-orb-ecut"></a>
    - *Type*: Real
    - *Description*: See also the entry [exx_hybrid_type](#exx-hybrid-type). This parameter is only relevant when exx_hybrid_type="opt_orb". A plane wave basis is used to optimize the radial ABFs. This parameter thus gives the cut-off of plane wave expansion, in Ry. A reasonable choice is 60.
    - *Default*: 0

    [back to top](#input-file)

- exx_opt_orb_tolerence<a id="exx-opt-orb-tolerence"></a>
    - *Type*: Real
    - *Description*: See also the entry [exx_hybrid_type](#exx-hybrid-type). This parameter is only relevant when exx_hybrid_type="opt_orb". exx_opt_orb_tolerence determines the threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.
    - *Default*: 0

    [back to top](#input-file)

### Molecular dynamics
This part of variables are used to control the molecular dynamics calculations.

- md_type<a id="md-type"></a>
    - *Type*: Integer
    - *Description*: control the ensemble to run md.
        - -1: FIRE method to relax;
        - 0: NVE ensemble;
        - 1: NVT ensemble with Nose Hoover Chain;
        - 2: NVT ensemble with Langevin method;
        - 3: NVT ensemble with Anderson thermostat;
        - 4: MSST method; 
    - *Default*: 1

    [back to top](#input-file)

- md_nstep<a id="md_nstep"></a>
    - *Type*: Integer
    - *Description*: the total number of md steps.
    - *Default*: 10

    [back to top](#input-file)

- md_ensolver<a id="md-ensolver"></a>
    - *Type*: String
    - *Description*: choose the energy solver for MD.
        - FP: First-Principles MD;
        - LJ: Leonard Jones potential;
        - DP: DeeP potential;
    - *Default*: FP

    [back to top](#input-file)

- md_restart<a id="md-restart"></a>
    - *Type*: Boolean
    - *Description*: to control whether restart md.
        - 0: When set to 0, ABACUS will calculate md normolly.
        - 1: When set to 1, ABACUS will calculate md from last step in your test before.
    - *Default*: 0

    [back to top](#input-file)
- md_dt<a id="md-dt"></a>
    - *Type*: Double
    - *Description*: This is the time step(fs) used in md simulation .
    - *Default*: 1

    [back to top](#input-file)
- md_tfirst & md_tlast<a id="md-t"></a>
    - *Type*: Double
    - *Description*: This is the temperature (K) used in md simulation, md_tlast’s default value is md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.
    - *Default*: No default

    [back to top](#input-file)

- md_dumpfreq<a id="md-dumpfreq"></a>
    - *Type*: Integer
    - *Description*:This is the frequence to dump md information.
    - *Default*: 1

    [back to top](#input-file)

- md_restartfreq<a id="md-restartfreq"></a>
    - *Type*: Integer
    - *Description*:This is the frequence to output restart information.
    - *Default*: 5

    [back to top](#input-file)

- md_tfreq<a id="md-tfreq"></a>
    - *Type*: Real
    - *Description*: 
        - Oscillation frequency, used to determine Qmass of NHC; 
        - 1/(md_tfreq*md_dt) is collision probability in Anderson method.
    - *Default*: 1.0

    [back to top](#input-file)

- md_mnhc<a id="md-mnhc"></a>
    - *Type*: Integer
    - *Description*: Number of Nose-Hoover chains.
    - *Default*: 4

    [back to top](#input-file)

- lj_rcut<a id="lj-rcut"></a>
    - *Type*: Real
    - *Description*: Cut-off radius for Leonard Jones potential (angstrom).
    - *Default*: 8.5 (for He)

    [back to top](#input-file)

- lj_epsilon<a id="lj-epsilon"></a>
    - *Type*: Real
    - *Description*: The value of epsilon for Leonard Jones potential (eV).
    - *Default*: 0.01032 (for He)

    [back to top](#input-file)

- lj_sigma<a id="lj-sigma"></a>
    - *Type*: Real
    - *Description*: The value of sigma for Leonard Jones potential (angstrom).
    - *Default*: 3.405 (for He)

    [back to top](#input-file)

- msst_direction<a id="msst-direction"></a>
    - *Type*: Integer
    - *Description*: the direction of shock wave for MSST.
    - *Default*: 2 (z direction)

    [back to top](#input-file)

- msst_vel<a id="msst-vel"></a>
    - *Type*: Real
    - *Description*: the velocity of shock wave (\AA/fs) for MSST.
    - *Default*: 0

    [back to top](#input-file)

- msst_vis<a id="msst-vis"></a>
    - *Type*: Real
    - *Description*: artificial viscosity (mass/length/time) for MSST.
    - *Default*: 0

    [back to top](#input-file)

- msst_tscale<a id="msst-tscale"></a>
    - *Type*: Real
    - *Description*: reduction in initial temperature (0~1) used to compress volume in MSST.
    - *Default*: 0

    [back to top](#input-file)

- msst_qmass<a id="msst-qmass"></a>
    - *Type*: Double
    - *Description*: Inertia of extended system variable. Used only when md_type is 4, you should set a number which is larger than 0. Note that Qmass of NHC is set by md_tfreq.
    - *Default*: No default

    [back to top](#input-file)

- md_damp<a id="md-damp"></a>
    - *Type*: Real
    - *Description*: damping parameter (fs) used to add force in Langevin method.
    - *Default*: 1

    [back to top](#input-file)

### DFT+U correction
This part of variables are used to control DFT+U correlated parameters
- dft_plus_u<a id="dft_plus_u"></a>
    - *Type*: Boolean
    - *Description*: If set to 1, ABCUS will calculate plus U correction, which is especially important for correlated electron.
    - *Default*: 0

    [back to top](#input-file)

- orbital_corr<a id="orbital_corr"></a>
    - *Type*: Int
    - *Description*: $l_1,l_2,l_3,\ldots$ for atom type 1,2,3 respectively.(usually 2 for d electrons and 3 for f electrons) .Specify which orbits need plus U correction for each atom. If set to -1, the correction would not be calculate for this atom.
    - *Default*: None

    [back to top](#input-file)

- hubbard_u<a id="hubbard_u"></a>
    - *Type*: Real
    - *Description*: Hubbard Coulomb interaction parameter U(ev) in plus U correction,which should be specified for each atom unless Yukawa potential is use. ABACUS use a simplified scheme which only need U and J for each atom.
    - *Default*: 0.0 

    [back to top](#input-file)

- hund_j<a id="hund_j"></a>
    - *Type*: Real
    - *Description*: Hund exchange parameter J(ev) in plus U correction ,which should be specified for each atom unless Yukawa potential is use. ABACUS use a simplified scheme which only need U and J for each atom.
    - *Default*: 0.0 

    [back to top](#input-file)

- yukawa_potential<a id="yukawa_potential"></a>
    - *Type*: Boolean
    - *Description*: whether use the local screen Coulomb potential method to calculate the value of U and J. If this is set to 1, hubbard_u and hund_j do not need to be specified.
    - *Default*: 0

    [back to top](#input-file)

- omc<a id="omc"></a> 
    - *Type*: Boolean
    - *Description*: whether turn on occupation matrix control method or not
    - *Default*: 0

    [back to top](#input-file)


### VdW correction
This part of variables are used to control vdW-corrected related parameters.

- vdw_method<a id="vdw-method"></a>
    - *Type*: String
    - *Description*: If set to d2 ,d3_0 or d3_bj, ABACUS will calculate corresponding vdW correction, which is DFT-D2, DFT-D3(0) or DFTD3(BJ) method. And this correction includes energy and forces. "none" means that no vdW-corrected method has been used.
    - *Default*: none

    [back to top](#input-file)
- vdw_s6<a id="vdw-s6"></a>
    - *Type*: Real
    - *Description*: This scale factor is to optimize the interaction energy deviations. For DFT-D2, it is found to be 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). For DFT-D3, recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: 0.75 if vdw_method is chosen to be d2; 1.0 if vdw_method is chosen to be d3_0 or d3_bj

    [back to top](#input-file)
- vdw_s8<a id="vdw-s8"></a>
    - *Type*: Real
    - *Description*: This scale factor is only the parameter of DFTD3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: 0.722 if vdw_method is chosen to be d3_0; 0.7875 if vdw_method is chosen to be d3_bj

    [back to top](#input-file)
-   vdw_a1<a id="vdw-a1"></a>
    - *Type*: Real
    - *Description*: This damping function parameter is only the parameter of DFT-D3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: 1.217 if vdw_method is chosen to be d3_0; 0.4289 if vdw_method is chosen to be d3_bj

    [back to top](#input-file)
- vdw_a2<a id="vdw-a2"></a>
    - *Type*: Real
    - *Description*: This damping function arameter is only the parameter of DFT-D3(BJ) approach. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: 1.0 if vdw_method is chosen to be d3_0; 4.4407 if vdw_method is chosen to be d3_bj

    [back to top](#input-file)
- vdw_d<a id="vdw-d"></a>
    - *Type*: Real
    - *Description*: The variable is to control the dumping speed of dumping function of DFT-D2.
    - *Default*: 20

    [back to top](#input-file)
- vdw_abc<a id="vdw-abc"></a>
    - *Type*: Integer
    - *Description*: The variable is to control the calculation of three-body term of DFT-D3 approachs, including D3(0) and D3(BJ). If set to 1, ABACUS will calculate three-body term, otherwise, the three-body term is not included.
    - *Default*: 0

    [back to top](#input-file)
- vdw_C6_file<a id="vdw-C6-file"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent’s C<sub>6</sub> Parameters file. If you don’t setup this, ABACUS will use the default C<sub>6</sub> Parameters stored in the programme already. We’ve stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new C<sub>6</sub> Parameters, you should provide a file contained all the C<sub>6</sub> Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
    - *Default*: default

    [back to top](#input-file)
- vdw_C6_unit<a id="vdw-C6-unit"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies unit of C<sub>6</sub> Parameters. Two kinds of unit is available: `J·nm6/mol` (means J·nm<sup>6</sup>/mol) and eVA(means eV·Angstrom)
    - *Default*: Jnm6/mol

    [back to top](#input-file)
- vdw_R0_file<a id="vdw-R0-file"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent’s R<sub>0</sub> Parameters file. If you don’t setup this, ABACUS will use the default R<sub>0</sub> Parameters stored in the programme already. We’ve stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new R<sub>0</sub> Parameters, you should provide a file contained all the R<sub>0</sub> Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
    - *Default*: default

    [back to top](#input-file)
- vdw_R0_unit<a id="vdw-R0-unit"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies unit of R<sub>0</sub> Parameters. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: A

    [back to top](#input-file)
- vdw_model<a id="vdw-model"></a>
    - *Type*: String
    - *Description*: To calculate the periodic structure, you can assign the number of lattice cells calculated. This variable specifies the kind of model to assign. If set to period, ABACUS will calculate a cubic periodic structure as assigned in vdw_period. If set to radius, ABACUS will calculate a cubic periodic structure containing a sphere, whose radius is vdw_radius and centre is origin point.
    - *Default*: radius

    [back to top](#input-file)
- vdw_radius<a id="vdw-radius"></a>
    - *Type*: Real
    - *Description*: If vdw_model is set to radius, this variable specifies the radius of the calculated sphere. For DFT-D2, the default value is 56.6918 ,while it is 95 for DFT-D3. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: 56.6918 if vdw_method is chosen to be d2; 95 if vdw_method is chosen to be d3_0 or d3_bj

    [back to top](#input-file)
- vdw_radius_unit<a id="vdw-radius-unit"></a>
    - *Type*: String
    - *Description*: If vdw_model is set to radius, this variable specifies the unit of vdw_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: Bohr

    [back to top](#input-file)
- vdw_cn_radius<a id="vdw-cn-radius"></a>
    - *Type*: Real
    - *Description*: This cutoff is chosen for the calculation of the coordination number (CN) in DFT-D3 approachs, including D3(0) and D3(BJ).
    - *Default*: 40.0

    [back to top](#input-file)
- vdw_cn_radius_unit<a id="vdw-cn-radius-unit"></a>
    - *Type*: String
    - *Description*: This variable specifies the unit of vdw_cn_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: Bohr

    [back to top](#input-file)
- vdw_period<a id="vdw-period"></a>
    - *Type*: Int Int Int
    - *Description*: If vdw_model is set to period, these variables specify the number of x, y and z periodic.
    - *Default*: 3 3 3

    [back to top](#input-file)

### Berry phase and wannier90 interface
This part of variables are used to control berry phase and wannier90 interfacae parameters.

- berry_phase<a id="berry-phase"></a>
    - *Type*: Integer
    - *Description*: 1, calculate berry phase; 0, no calculate berry phase.
    - *Default*: 0

    [back to top](#input-file)
- gdir<a id="gdir"></a>
    - *Type*: Integer
    - *Description*:
        - 1: calculate the polarization in the direction of the lattice vector a_1 that is defined in STRU file.
        - 2: calculate the polarization in the direction of the lattice vector a_2 that is defined in STRU file.
        - 3: calculate the polarization in the direction of the lattice vector a_3 that is defined in STRU file.
    - *Default*: 3

    [back to top](#input-file)
- towannier90<a id="towannier90"></a>
    - *Type*: Integer
    - *Description*: 1, generate files for wannier90 code; 0, no generate.
    - *Default*: 0

    [back to top](#input-file)
- nnkpfile<a id="nnkpfile"></a>
    - *Type*: String
    - *Description*: the file name when you run “wannier90 -pp ...”.
    - *Default*: seedname.nnkp

    [back to top](#input-file)
- wannier_spin<a id="wannier-spin"></a>
    - *Type*: String
    - *Description*: If nspin is set to 2,
        - up: calculate spin up for wannier function.
        - down: calculate spin down for wannier function.
    - *Default*: up

    [back to top](#input-file)

### TDDFT: time dependent density functional theory
- tddft<a id="tddft"></a>
    - *Type*: Integer
    - *Description*:
        - 1: calculate the real time time dependent density functional theory (TDDFT).
        - 0: do not calculate TDDFT.
    - *Default*: 0

    [back to top](#input-file)
- td_scf_thr_rho<a id="td_scf_thr_rho"></a>
    - *Type*: Double
    - *Description*: Accuracy of electron convergence when doing time-dependent evolution.
    - *Default*: 1e-9

    [back to top](#input-file)
- td_dt<a id="td_dt"></a>
    - *Type*: Double
    - *Description*: Time-dependent evolution time step. (fs)
    - *Default*: 0.02

    [back to top](#input-file)
- td_force_dt<a id="td_force_dt"></a>
    - *Type*: Double
    - *Description*: Time-dependent evolution force changes time step. (fs)
    - *Default*: 0.02

    [back to top](#input-file)
- td_vext<a id="td_vext"></a>
    - *Type*: Integer
    - *Description*:
        - 1: add a laser material interaction (extern laser field).
        - 0: no extern laser field.
    - *Default*: 0

    [back to top](#input-file)
- td_vext_dire<a id="td_vext_dire"></a>
    - *Type*: Integer
    - *Description*:
        - 1: the direction of external light field is along x axis.
        - 2: the direction of external light field is along y axis.
        - 3: the direction of external light field is along z axis.
    - *Default*: 1

    [back to top](#input-file)
- td_timescale<a id="td_timescale"></a>
    - *Type*: Double
    - *Description*: Time range of external electric field application. (fs)
    - *Default*: 0.5

    [back to top](#input-file)
- td_vexttype<a id="td_vexttype"></a>
    - *Type*: Integer
    - *Description*:
        - 1: Gaussian-type light field.
        - 2: Delta function form light field.
        - 3: Trigonometric function form light field.
    - *Default*: 1

    [back to top](#input-file)
- td_vextout<a id="td_vextout"></a>
    - *Type*: Integer
    - *Description*:
        - 1: Output external electric field.
        - 0: do not Output external electric field.
    - *Default*: 0

    [back to top](#input-file)
- td_dipoleout<a id="td_dipoleout"></a>
    - *Type*: Integer
    - *Description*:
        - 1: Output dipole.
        - 0: do not Output dipole.
    - *Default*: 0

    [back to top](#input-file)

### Variables useful for debugging

- nurse(#nurse)

    - *Type*: Boolean
    - *Description*: If set to 1, the Hamiltonian matrix and S matrix in each iteration will be written in output.
    - *Default*: 0

    [back to top](#input-file)

- t_in_h<a id="t-in-h"></a>

    - *Type*: Boolean
    - *Description*: If set to 0, then kinetic term will not be included in obtaining the Hamiltonian.
    - *Default*: 1

    [back to top](#input-file)

- vl_in_h<a id="vl-in-h"></a>

    - *Type*: Boolean
    - *Description*: If set to 0, then local pseudopotential term will not be included in obtaining the Hamiltonian.
    - *Default*: 1

    [back to top](#input-file)

- vnl_in_h<a id="vnl-in-h"></a>

    - *Type*: Boolean
    - *Description*:  If set to 0, then non-local pseudopotential term will not be included in obtaining the Hamiltonian.
    - *Default*: 1

    [back to top](#input-file)

- test_force<a id="test-force"></a>

    - *Type*: Boolean
    - *Description*: If set to 1, then detailed components in forces will be written to output.
    - *Default*: 0

    [back to top](#input-file)

- test_stress<a id="test-stress"></a>

    - *Type*: Boolean
    - *Description*: If set to 1, then detailed components in stress will be written to output.
    - *Default*: 0

    [back to top](#input-file)

- colour<a id="colour"></a>

    - *Type*: Boolean
    - *Description*: If set to 1, output to terminal will have some color.
    - *Default*: 0

    [back to top](#input-file)

- test_just_neighbor<a id="test-just-neighbor"></a>
    - *Type*: Boolean
    - *Description*: If set to 1, then only perform the neighboring atoms search.
    - *Default*: 0

    [back to top](#input-file)

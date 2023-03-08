# Full List of INPUT Keywords

- [System variables](#system-variables)

  [suffix](#suffix) | [ntype](#ntype) | [calculation](#calculation) | [esolver_type](#esolver_type) | [symmetry](#symmetry) | [kpar](#kpar) | [bndpar](#bndpar) | [latname](#latname) | [init_wfc](#init_wfc) | [init_chg](#init_chg) | [init_vel](#init_vel) | [nelec](#nelec) | [nupdown](#nupdown) | [dft_functional](#dft_functional) | [xc_temperature](#xc_temperature) | [pseudo_rcut](#pseudo_rcut) | [pseudo_mesh](#pseudo_mesh) | [mem_saver](#mem_saver) | [diago_proc](#diago_proc) | [nbspline](#nbspline) | [kspacing](#kspacing)  | [min_dist_coef](#min_dist_coef) | [symmetry_prec](#symmetry_prec) | [device](#device)
- [Variables related to input files](#variables-related-to-input-files)

  [stru_file](#stru_file) | [kpoint_file](#kpoint_file) | [pseudo_dir](#pseudo_dir) | [orbital_dir](#orbital_dir) | [read_file_dir](#read_file_dir) | [wannier_card](#wannier_card)
- [Plane wave related variables](#plane-wave-related-variables)

  [ecutwfc](#ecutwfc) | [nx,ny,nz](#nx-ny-nz) | [pw_seed](#pw_seed) | [pw_diag_thr](#pw_diag_thr) | [pw_diag_nmax](#pw_diag_nmax) | [pw_diag_ndim](#pw_diag_ndim)
- [Numerical atomic orbitals related variables](#numerical-atomic-orbitals-related-variables)

  [nb2d](#nb2d) | [lmaxmax](#lmaxmax) | [lcao_ecut](#lcao_ecut) | [lcao_dk](#lcao_dk) | [lcao_dr](#lcao_dr) | [lcao_rmax](#lcao_rmax) | [search_radius](#search_radius) | [search_pbc](#search_pbc) | [bx,by,bz](#bx-by-bz)
- [Electronic structure](#electronic-structure)

  [basis_type](#basis_type) | [ks_solver](#ks_solver) | [nbands](#nbands) | [nbands_istate](#nbands_istate) | [nspin](#nspin) | [smearing_method](#smearing_method) | [smearing_sigma](#smearing_sigma) | [smearing_sigma_temp](#smearing_sigma_temp) | [mixing_type](#mixing_type) | [mixing_beta](#mixing_beta) | [mixing_ndim](#mixing_ndim) | [mixing_gg0](#mixing_gg0) | [mixing_tau](#mixing_tau) | [mixing_dftu](#mixing_dftu) | [gamma_only](#gamma_only) | [printe](#printe) | [scf_nmax](#scf_nmax) | [scf_thr](#scf_thr) | [chg_extrap](#chg_extrap) | [lspinorb](#lspinorb) | [noncolin](#noncolin) | [soc_lambda](#soc_lambda)
- [Electronic structure (SDFT)](#electronic-structure-sdft)

  [method_sto](#method_sto) | [nbands_sto](#nbands_sto) | [nche_sto](#nche_sto) | [emin_sto](#emin_sto) | [emax_sto](#emax_sto) | [seed_sto](#seed_sto) | [initsto_freq](#initsto_freq) | [npart_sto](#npart_sto)
- [Geometry relaxation](#geometry-relaxation)

  [relax_nmax](#relax_nmax) | [relax_method](#relax_method) | [relax_cg_thr](#relax_cg_thr) | [relax_bfgs_w1](#relax_bfgs_w1) | [relax_bfgs_w2](#relax_bfgs_w2) | [relax_bfgs_rmax](#relax_bfgs_rmax) | [relax_bfgs_rmin](#relax_bfgs_rmin) | [relax_bfgs_init](#relax_bfgs_init) | [cal_force](#cal_force) | [force_thr](#force_thr) | [force_thr_ev](#force_thr_ev) | [force_thr_ev2](#force_thr_ev2) | [cal_stress](#cal_stress) | [stress_thr](#stress_thr) | [press1, press2, press3](#press1-press2-press3) | [fixed_axes](#fixed_axes) | [cell_factor](#cell_factor) | [fixed_ibrav](#fixed_ibrav) | [relax_new](#relax_new) | [relax_scale_force](#relax_scale_force)
- [Variables related to output information](#variables-related-to-output-information)

  [out_mul](#out_mul) | [out_freq_elec](#out_freq_elec) | [out_freq_ion](#out_freq_ion) | [out_chg](#out_chg) | [out_pot](#out_pot) | [out_dm](#out_dm) | [out_wfc_pw](#out_wfc_pw) | [out_wfc_r](#out_wfc_r) | [out_wfc_lcao](#out_wfc_lcao) | [out_dos](#out_dos) | [out_band](#out_band) | [out_proj_band](#out_proj_band) | [out_stru](#out_stru) | [out_bandgap](#out_bandgap) | [out_level](#out_level) | [out_alllog](#out_alllog) | [out_mat_hs](#out_mat_hs) | [out_mat_r](#out_mat_r) | [out_mat_hs2](#out_mat_hs2) | [out_mat_t](#out_mat_t) | [out_mat_dh](#out_mat_dh) | [out_hs2_interval](#out_hs2_interval) | [out_element_info](#out_element_info) | [restart_save](#restart_save) | [restart_load](#restart_load) | [dft_plus_dmft](#dft_plus_dmft) | [rpa](#rpa)
- [Density of states](#density-of-states)

  [dos_edelta_ev](#dos_edelta_ev) | [dos_sigma](#dos_sigma) | [dos_scale](#dos_scale) | [dos_emin_ev](#dos_emin_ev) | [dos_emax_ev](#dos_emax_ev) | [dos_nche](#dos_nche)
- [Exact exchange](#exact-exchange) (Under tests)

  [exx_hybrid_alpha](#exx_hybrid_alpha) | [exx_hse_omega](#exx_hse_omega) | [exx_separate_loop](#exx_separate_loop) | [exx_hybrid_step](#exx_hybrid_step) | [exx_lambda](#exx_lambda) | [exx_pca_threshold](#exx_pca_threshold) | [exx_c_threshold](#exx_c_threshold) | [exx_v_threshold](#exx_v_threshold) | [exx_c_grad_threshold](#exx_c_grad_threshold) | [exx_v_grad_threshold](#exx_v_grad_threshold) | [exx_dm_threshold](#exx_dm_threshold) | [exx_schwarz_threshold](#exx_schwarz_threshold) | [exx_cauchy_threshold](#exx_cauchy_threshold) | [exx_cauchy_grad_threshold](#exx_cauchy_grad_threshold) | [exx_ccp_threshold](#exx_ccp_threshold) | [exx_ccp_rmesh_times](#exx_ccp_rmesh_times) | [exx_distribute_type](#exx_distribute_type) | [exx_opt_orb_lmax](#exx_opt_orb_lmax) | [exx_opt_orb_ecut](#exx_opt_orb_ecut) | [exx_opt_orb_tolerence](#exx_opt_orb_tolerence) | [exx_real_number](#exx_real_number)
- [Molecular dynamics](#molecular-dynamics)

  [md_type](#md_type) | [md_thermostat](#md_thermostat) | [md_nstep](#md_nstep) | [md_restart](#md_restart) | [md_dt](#md_dt) | [md_tfirst, md_tlast](#md_tfirst-md_tlast) | [md_dumpfreq](#md_dumpfreq) | [md_restartfreq](#md_restartfreq) | [md_seed](#md_seed) | [md_tfreq](#md_tfreq) | [md_tchain](#md_tchain) | [md_pmode](#md_pmode) | [md_pcouple](#md_pcouple) | [md_pfirst, md_plast](#md_pfirst-md_plast) | [md_pfreq](#md_pfreq) | [md_pchain](#md_pchain) | [out_force](#out_force) | [out_vel](#out_vel) | [out_virial](#out_virial) | [lj_rcut](#lj_rcut) | [lj_epsilon](#lj_epsilon) | [lj_sigma](#lj_sigma) | [pot_file](#pot_file) | [msst_direction](#msst_direction) | [msst_vel](#msst_vel) | [msst_vis](#msst_vis) | [msst_tscale](#msst_tscale) | [msst_qmass](#msst_qmass) | [md_damp](#md_damp) | [md_tolerance](#md_tolerance) | [md_nraise](#md_nraise)
- [vdW correction](#vdw-correction)

  [vdw_method](#vdw_method) | [vdw_s6](#vdw_s6) | [vdw_s8](#vdw_s8) | [vdw_a1](#vdw_a1) | [vdw_a2](#vdw_a2) | [vdw_d](#vdw_d) | [vdw_abc](#vdw_abc) | [vdw_C6_file](#vdw_c6_file) | [vdw_C6_unit](#vdw_c6_unit) | [vdw_R0_file](#vdw_r0_file) | [vdw_R0_unit](#vdw_r0_unit) | [vdw_cutoff_type](#vdw_cutoff_type) | [vdw_cutoff_radius](#vdw_cutoff_radius) | [vdw_radius_unit](#vdw_radius_unit) | [vdw_cutoff_period](#vdw_cutoff_period) | [vdw_cn_thr](#vdw_cn_thr) | [vdw_cn_thr_unit](#vdw_cn_thr_unit)
- [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)

  [berry_phase](#berry_phase) | [gdir](#gdir) | [towannier90](#towannier90) | [nnkpfile](#nnkpfile) | [wannier_spin](#wannier_spin)
- [TDDFT: time dependent density functional theory](#tddft-time-dependent-density-functional-theory) (Under tests)

  [td_force_dt](#td_force_dt) | [td_vext](#td_vext) | [td_vext_dire](#td_vext_dire) | [td_stype](#td_stype) | [td_ttype](#td_ttype) | [td_tstart](#td_tstart) | [td_tend](#td_tend) | [td_lcut1](#td_lcut1) | [td_lcut2](#td_lcut2) | [td_gauss_freq](#td_gauss_freq) | [td_guass_phase](#td_gauss_phase) | [td_gauss_sigma](#td_gauss_sigma) | [td_gauss_t0](#td_gauss_t0)| [td_gauss_amp](#td_gauss_amp) | [td_trape_freq](#td_trape_freq) | [td_trape_phase](#td_trape_phase) | [td_trape_t1](#td_trape_t1) | [td_trape_t2](#td_trape_t2) | [td_trape_t3](#td_trape_t3) | [td_trape_amp](#td_trape_amp) | [td_trigo_freq1](#td_trigo_freq1) | [td_trigo_freq2](#td_trigo_freq2) | [td_trigo_phase1](#td_trigo_phase1) | [td_trigo_phase2](#td_trigo_phase2) | [td_trigo_amp](#td_trigo_amp) | [td_heavi_t0](#td_heavi_t0) | [td_heavi_amp](#td_heavi_amp) | [out_dipole](#out_dipole) |[out_efield](#out_efield)| [ocp](#ocp) | [ocp_set](#ocp_set) | [td_val_elec_01](#td_val_elec_01) | [td_val_elec_02](#td_val_elec_02) |[td_val_elec_03](#td_val_elec_03)
- [DFT+*U* correction](#dftu-correction) (Under development)

  [dft_plus_u](#dft_plus_u) | [orbital_corr](#orbital_corr) | [hubbard_u](#hubbard_u) | [yukawa_potential](#yukawa_potential) | [yukawa_lambda](#yukawa_lambda) | [omc](#omc)
- [Variables useful for debugging](#variables-useful-for-debugging)

  [nurse](#nurse) | [t_in_h](#t_in_h) | [vl_in_h](#vl_in_h) | [vnl_in_h](#vnl_in_h) | [vh_in_h](#vh_in_h) | [vion_in_h](#vion_in_h) | [test_force](#test_force) | [test_stress](#test_stress) | [colour](#colour) | [test_skip_ewald](#test_skip_ewald)
- [NAOs](#naos)

  [bessel_nao_ecut](#bessel_nao_ecut) | [bessel_nao_tolerence](#bessel_nao_tolerence) | [bessel_nao_rcut](#bessel_nao_rcut) | [bessel_nao_smooth](#bessel_nao_smooth) | [bessel_nao_sigma](#bessel_nao_sigma)
- [DeePKS](#deepks)

  [deepks_out_labels](#deepks_out_labels) | [deepks_scf](#deepks_scf) | [deepks_model](#deepks_model) | [bessel_descriptor_lmax](#bessel_descriptor_lmax) | [bessel_descriptor_ecut](#bessel_descriptor_ecut) | [bessel_descriptor_tolerence](#bessel_descriptor_tolerence) | [bessel_descriptor_rcut](#bessel_descriptor_rcut) | [bessel_descriptor_smooth](#bessel_descriptor_smooth) | [bessel_descriptor_sigma](#bessel_descriptor_sigma) | [deepks_bandgap](#deepks_bandgap) | [deepks_out_unittest](#deepks_out_unittest)
- [OFDFT: orbital free density functional theory](#ofdft-orbital-free-density-functional-theory)

  [of_kinetic](#of_kinetic) | [of_method](#of_method) | [of_conv](#of_conv) | [of_tole](#of_tole) | [of_tolp](#of_tolp) | [of_tf_weight](#of_tf_weight) | [of_vw_weight](#of_vw_weight) | [of_wt_alpha](#of_wt_alpha) | [of_wt_beta](#of_wt_beta) | [of_wt_rho0](#of_wt_rho0) | [of_hold_rho0](#of_hold_rho0) | [of_read_kernel](#of_read_kernel) | [of_kernel_file](#of_kernel_file) | [of_full_pw](#of_full_pw) | [of_full_pw_dim](#of_full_pw_dim)
- [Electric field and dipole correction](#electric-field-and-dipole-correction)

  [efield_flag](#efield_flag) | [dip_cor_flag](#dip_cor_flag) | [efield_dir](#efield_dir) | [efield_pos_max](#efield_pos_max) | [efield_pos_dec](#efield_pos_dec) | [efield_amp](#efield_amp)
- [Gate field (compensating charge)](#gate-field-compensating-charge)

  [gate_flag](#gate_flag) | [zgate](#zgate) | [block](#block) | [block_down](#block_down) | [block_up](#block_up) | [block_height](#block_height)
- [Electronic conductivities](#electronic-conductivities)

  [cal_cond](#cal_cond) | [cond_nche](#cond_nche) | [cond_dw](#cond_dw) | [cond_wcut](#cond_wcut) | [cond_wenlarge](#cond_wenlarge) | [cond_fwhm](#cond_fwhm) | [cond_nonlocal](#cond_nonlocal)
- [Implicit solvation model](#implicit-solvation-model)

  [imp_sol](#imp_sol) | [eb_k](#eb_k) | [tau](#tau) | [sigma_k](#sigma_k) | [nc_k](#nc_k)

[back to top](#full-list-of-input-keywords)

## System variables

These variables are used to control general system parameters.

### suffix

- **Type**: String
- **Description**: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the `suffix` is the name you can pick up for your convenience.
- **Default**: ABACUS

### ntype

- **Type**: Integer
- **Description**: Number of different atom species in this calculation. If this value is not equal to the atom species in the STRU file, ABACUS will stop and quit. If not set or set to 0, ABACUS will automatically set it to the atom species in the STRU file.
- **Default**: 0

### calculation

- **Type**: String
- **Description**: Specify the type of calculation.

  - *scf*: do self-consistent electronic structure calculation
  - *relax*: do structure relaxation calculation, one can use `relax_nmax` to decide how many ionic relaxations you want.
  - *cell-relax*: do variable-cell relaxation calculation.
  - *nscf*: do the non self-consistent electronic structure calculations. For this option, you need a charge density file. For nscf calculations with planewave basis set, pw_diag_thr should be <= 1e-3.
  - *istate*: For LCAO basis. Please see the explanation for variable `nbands_istate`.
  - *ienvelope*: Envelope function for LCAO basis. Please see the explanation for variable `nbands_istate`.
  - *md*: molecular dynamics
  - *test_memory* : checks memory required for the calculation. The number is not quite reliable, please use it with care
  - *test_neighbour* : only performs neighbouring atom search
  - *gen_bessel* : generates projectors (a series of Bessel functions) for DeePKS; see also keywords bessel_descriptor_lmax, bessel_descriptor_rcut and bessel_descriptor_tolerence. A file named `jle.orb` will be generated which contains the projectors. An example is provided in examples/H2O-deepks-pw.
  - *get_S* : only works for multi-k calculation with LCAO basis. Generates and writes the overlap matrix to a file names `SR.csr` in the working directory. The format of the file will be the same as that generated by [out_mat_hs2](#out_mat_hs2).
- **Default**: scf

### esolver_type

- **Type**: String
- **Description**: choose the energy solver.
  - ksdft: Kohn-Sham density functional theory;
  - ofdft: orbital-free density functional theory;
  - sdft: [stochastic density functional theory](#electronic-structure-sdft);
  - tddft: real-time time-dependent density functional theory (TDDFT);
  - lj: Leonard Jones potential;
  - dp: DeeP potential;
- **Default**: ksdft

### symmetry

- **Type**: Integer
- **Description**: takes value 1, 0 or -1.
  - if set to 1, symmetry analysis will be performed to determine the type of Bravais lattice and associated symmetry operations. (point groups only)
  - if set to 0, only time reversal symmetry would be considered in symmetry operations, which implied k point and -k point would be treated as a single k point with twice the weight.
  - if set to -1, no symmetry will be considered.
- **Default**: 0

### kpar

- **Type**: Integer
- **Description**: divide all processors into kpar groups, and k points will be distributed among each group. The value taken should be less than or equal to the number of k points as well as the number of MPI threads.
- **Default**: 1

### bndpar

- **Type**: Integer
- **Description**: divide all processors into bndpar groups, and bands (only stochastic orbitals now) will be distributed among each group. It should be larger than 0.
- **Default**: 1

### latname

- **Type**: String
- **Description**: Specifies the type of Bravias lattice. When set to `none`, the three lattice vectors are supplied explicitly in STRU file. When set to a certain Bravais lattice type, there is no need to provide lattice vector, but a few lattice parameters might be required. For more information regarding this parameter, consult the [page on STRU file](stru.md).
  Available options are (correspondence with ibrav in QE is given in parenthesis):
  - `none`: free structure.
  - `sc`: simple cubic. (1)
  - `fcc`: face-centered cubic. (2)
  - `bcc`: body-centered cubic. (3)
  - `hexagonal`: hexagonal. (4)
  - `trigonal`: trigonal. (5)
  - `st`: simple tetragonal. (6)
  - `bct`: body-centered tetragonal. (7)
  - `so`: orthorhombic. (8)
  - `baco`: base-centered orthorhombic. (9)
  - `fco`: face-centered orthorhombic. (10)
  - `bco`: body-centered orthorhombic. (11)
  - `sm`: simple monoclinic. (12)
  - `bacm`: base-centered monoclinic. (13)
  - `triclinic`: triclinic. (14)
- **Default**: `none`

### init_wfc

- **Type**: String
- **Description**: Only useful for plane wave basis only now. It is the name of the starting wave functions. In the future. we should also make this variable available for localized orbitals set.
  Available options are:
  - `atomic`: from atomic pseudo wave functions. If they are not enough, other wave functions are initialized with random numbers.
  - `atomic+random`: add small random numbers on atomic pseudo-wavefunctions
  - `file`: from file
  - `random`: random numbers
- **Default**:`atomic`

### init_chg

- **Type**: String
- **Description**: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density. If set to `atomic`, the density is starting from the summation of the atomic density of single atoms. If set this to `file`, the density will be read in from a file. Besides, when you do `nspin=1` calculation, you only need the density file SPIN1_CHGCAR. However, if you do `nspin=2` calculation, you also need the density file SPIN2_CHGCAR. The density file should be output with these names if you set out_chg = 1 in INPUT file.
- **Default**: atomic

### init_vel

- **Type**: Boolean
- **Description**: Read the atom velocity from the atom file (STRU) if set to true.
- **Default**: false

### nelec

- **Type**: Real
- **Description**: If >0.0, this denotes the total number of electrons in the system. Must be less than 2*nbands. If set to 0.0, the total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
- **Default**: 0.0

### nupdown

- **Type**: Real
- **Description**: If >0.0, this denotes the difference number of electrons between spin-up and spin-down in the system. The range of value must in [-nelec ~ nelec]. It is one method of constraint DFT, the fermi energy level will separate to E_Fermi_up and E_Fermi_down. If set to 0.0, no constrain apply to system.
- **Default**: 0.0

### dft_functional

- **Type**: String
- **Description**: In our package, the XC functional can either be set explicitly using the `dft_functional` keyword in `INPUT` file. If `dft_functional` is not specified, ABACUS will use the xc functional indicated in the pseudopotential file.
  On the other hand, if dft_functional is specified, it will overwrite the functional from pseudopotentials and performs calculation with whichever functional the user prefers. We further offer two ways of supplying exchange-correlation functional. The first is using 'short-hand' names such as 'LDA', 'PBE', 'SCAN'. A complete list of 'short-hand' expressions can be found in [the source code](../../../source/module_hamilt_general/module_xc/xc_functional.cpp). The other way is only available when ***compiling with LIBXC***, and it allows for supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by a plus sign, for example, 'dft_functional='LDA_X_1D_EXPONENTIAL+LDA_C_1D_CSC'. The list of LIBXC keywords can be found on its [website](https://www.tddft.org/programs/libxc/functionals/). In this way, **we support all the LDA,GGA and mGGA functionals provided by LIBXC**.

  Furthermore, the old INPUT parameter exx_hybrid_type for hybrid functionals has been absorbed into dft_functional. Options are `hf` (pure Hartree-Fock), `pbe0`(PBE0), `hse` (Note: in order to use HSE functional, LIBXC is required). Note also that HSE has been tested while PBE0 has NOT been fully tested yet, and the maximum CPU cores for running exx in parallel is $N(N+1)/2$, with N being the number of atoms. And forces for hybrid functionals are not supported yet.

  If set to `opt_orb`, the program will not perform hybrid functional calculation. Instead, it is going to generate opt-ABFs as discussed in this [article](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c00481).
- **Default**: same as UPF file.

### xc_temperature

- **Type**: Real
- **Description**: specifies temperature when using temperature-dependent XC functionals (KSDT and so on); unit in Rydberg
- **Default** : 0.0

### pseudo_rcut

- **Type**: Real
- **Description**: Cut-off of radial integration for pseudopotentials, in Bohr.
- **Default**: 15

### pseudo_mesh

- **Type**: Integer
- **Description**: If set to 0, then use our own mesh for radial integration of pseudopotentials; if set to 1, then use the mesh that is consistent with quantum espresso.
- **Default**: 0

### mem_saver

- **Type**: Boolean
- **Description**: Used only for nscf calculations. If set to 1, then a memory saving technique will be used for many k point calculations.
- **Default**: 0

### diago_proc

- **Type**: Integer
- **Description**: If set to a positive number, then it specifies the number of threads used for carrying out diagonalization. Must be less than or equal to total number of MPI threads. Also, when cg diagonalization is used, diago_proc must be the same as the total number of MPI threads. If set to 0, then it will be set to the number of MPI threads. Normally, it is fine just leave it to the default value. Only used for pw base.
- **Default**: 0

### nbspline

- **Type**: Integer
- **Description**: If set to a natural number, a Cardinal B-spline interpolation will be used to calculate Structure Factor. `nbspline` represents the order of B-spline basis and a larger one can get more accurate results but cost more.
  It is turned off by default.
- **Default**: -1

### kspacing

- **Type**: Real
- **Description**: Set the smallest allowed spacing between k points, unit in 1/bohr. It should be larger than 0.0, and suggest smaller than 0.25. When you have set this value > 0.0, then the KPT file is unnecessary, and the number of K points nk_i = max(1, int(|b_i|/KSPACING)+1), where b_i is the reciprocal lattice vector. The default value 0.0 means that ABACUS will read the applied KPT file. Notice: if gamma_only is set to be true, kspacing is invalid.
- **Default**: 0.0

### min_dist_coef

- **Type**: Real
- **Description**: a factor related to the allowed minimum distance between two atoms. At the beginning, ABACUS will check the structure, and if the distance of two atoms is shorter than min_dist_coef*(standard covalent bond length), we think this structure is unreasonable. If you want to calculate some structures in extreme conditions like high pressure, you should set this parameter as a smaller value or even 0.
- **Default**: 0.2

### symmetry_prec

- **Type**: Real
- **Description**: The accuracy for symmetry judgment. The unit is Bohr.
- **Default**: 1.0e-5

### device

- **Type**: String
- **Description**: Specifies the computing device for ABACUS.

  Available options are:

  - `cpu`: for CPUs via Intel, AMD, or Other supported CPU devices
  - `gpu`: for GPUs via CUDA.

  Known limitations:

  - `pw basis`: required by the `gpu` acceleration options
  - `cg ks_solver`: required by the `gpu` acceleration options
- **Default**: `cpu`

[back to top](#full-list-of-input-keywords)

## Variables related to input files

These variables are used to control parameters related to input files.

### stru_file

- **Type**: String
- **Description**: This parameter specifies the name of structure file which contains various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
- **Default**: STRU

### kpoint_file

- **Type**: String
- **Description**: This parameter specifies the name of k-points file. Note that if you use atomic orbitals as basis, and you only use gamma point, you don't need to have k-point file in your directory, ABACUS will automatically generate `KPT` file. Otherwise, if you use more than one k-point, please do remember the algorithm in ABACUS is different for gamma only and various k-point dependent simulations. So first you should turn off the k-point algorithm by set `gamma_only = 0` in `INPUT` and then you should setup your own k-points file.
- **Default**: KPT

### pseudo_dir

- **Type**: String
- **Description**: This parameter specifies pseudopotential directory.
- **Default**: ./

### orbital_dir

- **Type**: String
- **Description**: This parameter specifies orbital file directory.
- **Default**: ./

### read_file_dir

- **Type**: String
- **Description**: when the program needs to read files such as electron density(`SPIN1_CHG`) as a starting point, this variable tells the location of the files. For example, './' means the file is located in the working directory.
- **Default**: OUT.$suffix

### wannier_card

- **Type**: String
- **Description**: Relevant when using ABACUS with wannier90. Tells the name of the input file related to wannier90.
- **Default**: "none"

[back to top](#full-list-of-input-keywords)

## Plane wave related variables

These variables are used to control the plane wave related parameters.

### ecutwfc

- **Type**: Real
- **Description**: Energy cutoff for plane wave functions, the unit is **Rydberg**. Note that even for localized orbitals basis, you still need to setup an energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set, etc. Also, because our orbitals are generated by matching localized orbitals to a chosen set of wave functions from a certain energy cutoff, this set of localize orbitals is most accurate under this same plane wave energy cutoff.
- **Default**: 50

### nx, ny, nz

- **Type**: Integer
- **Description**: If set to a positive number, then the three variables specify the numbers of FFT grid points in x, y, z directions, respectively. If set to 0, the number will be calculated from ecutwfc.
- **Default**: 0

### pw_seed

- **Type**: Integer
- **Description**: Only useful for plane wave basis only now. It is the random seed to initialize wave functions. Only positive integers are avilable.
- **Default**:0

### pw_diag_thr

- **Type**: Real
- **Description**: Only used when you use `diago_type = cg` or `diago_type = david`. It indicates the threshold for the first electronic iteration, from the second iteration the pw_diag_thr will be updated automatically. **For nscf calculations with planewave basis set, pw_diag_thr should be <= 1e-3.**
- **Default**: 0.01

### pw_diag_nmax

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = cg` or `ks_solver = dav`. It indicates the maximal iteration number for cg/david method.
- **Default**: 40

### pw_diag_ndim

- **Type**: Integer
- **Description**: Only useful when you use `ks_solver = dav`. It indicates the maximal dimension for the Davidson method.
- **Default**: 4

[back to top](#full-list-of-input-keywords)

## Numerical atomic orbitals related variables

These variables are used to control the numerical atomic orbitals related parameters.

### nb2d

- **Type**: Integer
- **Description**: In LCAO calculations, we arrange the total number of processors in an 2D array, so that we can partition the wavefunction matrix (number of bands*total size of atomic orbital basis) and distribute them in this 2D array. When the system is large, we group processors into sizes of nb2d, so that multiple processors take care of one row block (a group of atomic orbitals) in the wavefunction matrix. If set to 0, nb2d will be automatically set in the program according to the size of atomic orbital basis:
  - if size <= 500 : nb2d = 1
  - if 500 < size <= 1000 : nb2d = 32
  - if size > 1000 : nb2d = 64;
- **Default**: 0

### lmaxmax

- **Type**: Integer
- **Description**: If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.
- **Default**: 2.

### lcao_ecut

- **Type**: Real
- **Description**: Energy cutoff when calculating LCAO two-center integrals. In Ry.
- **Default**: 50

### lcao_dk

- **Type**: Real
- **Description**: Delta k for 1D integration in LCAO
- **Default**: 0.01

### lcao_dr

- **Type**: Real
- **Description**: Delta r for 1D integration in LCAO
- **Default**: 0.01

### lcao_rmax

- **Type**: Real
- **Description**: Max R for 1D two-center integration table
- **Default**: 30

### search_radius

- **Type**: Real
- **Description**: Set the search radius for finding neighbouring atoms. If set to -1, then the radius will be set to maximum of projector and orbital cut-off.
- **Default**: -1

### search_pbc

- **Type**: Boolean
- **Description**: In searching for neighbouring atoms, if set to 1, then periodic images will also be searched. If set to 0, then periodic images will not be searched.
- **Default**: 1

### bx, by, bz

- **Type**: Integer
- **Description**: In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed.
- **Default**: 2

[back to top](#full-list-of-input-keywords)

## Electronic structure

These variables are used to control the electronic structure and geometry relaxation
calculations.

### basis_type

- **Type**: String
- **Description**: This is an important parameter to choose basis set in ABACUS.
  - *pw*: Using plane-wave basis set only.
  - *lcao_in_pw*: Expand the localized atomic set in plane-wave basis.
  - lcao: Using localized atomic orbital sets.
- **Default**: pw

### ks_solver

- **Type**: String
- **Description**: It's about the choice of diagonalization methods for the hamiltonian matrix expanded in a certain basis set.

  For plane-wave basis,

  - cg: cg method.
  - dav: the Davidson algorithm. (Currently not working with Intel MKL library).

  For atomic orbitals basis,

  - genelpa: This method should be used if you choose localized orbitals.
  - scalapack-gvx: scalapack can also be used for localized orbitals.
  - cusolver: this method needs building with the cusolver component for lcao and at least one gpu is available.

  If you set ks_solver=`genelpa` for basis_type=`pw`, the program will be stopped with an error message:

  ```
  genelpa can not be used with plane wave basis.
  ```

  Then the user has to correct the input file and restart the calculation.
- **Default**: `cg` (pw) or `genelpa` (lcao)

### nbands

- **Type**: Integer
- **Description**: Number of Kohn-Sham orbitals to calculate. It is recommended you setup this value, especially when you use smearing techniques, more bands should be included.
- **Default**:
  - nspin=1: 1.2\*occupied_bands, occupied_bands + 10)
  - nspin=2: max(1.2\*nelec_spin, nelec_spin + 10) , nelec_spin = max(nelec_spin_up, nelec_spin_down)
  - nspin=4: 1.2\*nelec, nelec + 20)

### nbands_istate

- **Type**: Integer
- **Description**: Only used when `calculation = ienvelope` or `calculation = istate`, this variable indicates how many bands around the Fermi level you would like to calculate. `ienvelope` means to calculate the envelope functions of wave functions $\Psi_{i}=\Sigma_{\mu}C_{i\mu}\Phi_{\mu}$, where $\Psi_{i}$ is the ith wave function with the band index $i$ and $\Phi_{\mu}$ is the localized atomic orbital set. `istate` means to calculate the density of each wave function $|\Psi_{i}|^{2}$. Specifically, suppose we have highest occupied bands at 100th wave functions. And if you set this variable to 5, it will print five wave functions from 96th to 105th. But before all this can be carried out, the wave functions coefficients  should be first calculated and written into a file by setting the flag `out_wfc_lcao = 1`.
- **Default**: 5

### nspin

- **Type**: Integer
- **Description**: Number of spin components of wave functions. There are only two choices now: 1 or 2, meaning non spin or collinear spin. For the case of [noncollinear polarized](../scf/spin.md#noncollinear-spin-polarized-calculations), nspin will be automatically set to 4 without being specified in user input.
- **Default**: 1

### smearing_method

- **Type**: String
- **Description**: It indicates which occupation and smearing method is used in the calculation.
  - fixed: use fixed occupations.
  - gauss or gaussian: use Gaussian smearing method.
  - mp: use methfessel-paxton smearing method; recommended for metals.
  - fd: Fermi-Dirac smearing method: $f=1/\{1+\exp[(E-\mu)/kT]\}$ and smearing_sigma below is the temperature $T$ (in Ry).
- **Default**: fixed

### smearing_sigma

- **Type**: Real
- **Description**: energy range for smearing, the unit is Rydberg.
- **Default**: 0.001

### smearing_sigma_temp

- **Type**: Real
- **Description**: energy range for smearing, and is the same as smearing_sigma, but the unit is K. smearing_sigma = 1/2 * kB * smearing_sigma_temp.

### mixing_type

- **Type**: String
- **Description**: Charge mixing methods. We offer the following 3 options:
  - plain: Just simple mixing.
  - pulay: Standard Pulay method.
  - broyden: Broyden method.
- **Default**: pulay

### mixing_beta

- **Type**: Real
- **Description**: mixing parameter: 0 means no new charge
- **Default**: 0.7

### mixing_ndim

- **Type**: Integer
- **Description**: It indicates the mixing dimensions in Pulay, Pulay method uses the density from previous mixing_ndim steps and do a charge mixing based on this density.
- **Default**: 8

### mixing_gg0

- **Type**: Real
- **Description**: When set to a positive number, the high frequency wave vectors will be suppressed by multiplying a scaling factor $\frac{k^2}{k^2+gg0^2}$; if set to 0, then no Kerker scaling is performed. Setting mixing_gg0 = 1.5 is normally a good starting point.
- **Default**: 0.0

### mixing_tau
- **Type**: Boolean
- **Description**: Only relevant for meta-GGA calculations. If set to true, then the kinetic energy density will also be mixed. It seems for general cases, SCF converges fine even without this mixing. However, if there is difficulty in converging SCF for meta-GGA, it might be helpful to turn this on.
- **Default**: False

### mixing_dftu
- **Type**: Boolean
- **Description**: Only relevant for DFT+U calculations. If set to true, then the occupation matrices will also be mixed by plain mixing. From experience this is not very helpful if the +U calculation does not converge.
- **Default**: False

### gamma_only

- **Type**: Integer
- **Description**: It is an important parameter **only to be used in localized orbitals set**.
  If you set gamma_only = 1, ABACUS uses gamma only, the algorithm is faster and you don't need to specify the k-points file. If you set gamma_only = 0, more than one k-point is used and the ABACUS is slower compared to the gamma only algorithm.

  > Note: If gamma_only is set to 1, the KPT file will be overwritten. So make sure to turn off gamma_only for multi-k calculations.
  >
- **Default**: 0

### printe

- **Type**: Integer
- **Description**: Print out energy for each band for every printe step
- **Default**: 100

### scf_nmax

- **Type**: Integer
- **Description**: This variable indicates the maximal iteration number for electronic iterations.
- **Default**: 100

### scf_thr

- **Type**: Real
- **Description**: An important parameter in ABACUS. It's the threshold for electronic iteration. It represents the charge density error between two sequential densities from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
- **Default**: 1.0e-9

### chg_extrap

- **Type**: String
- **Description**: Methods to do extrapolation of density when ABACUS is doing geometry relaxations.
  - atomic: atomic extrapolation
  - first-order: first-order extrapolation
  - second-order: second-order extrapolation
- **Default**: atomic

### lspinorb

- **Type**: Boolean
- **Description**: whether to consider spin-orbital coupling effect in the calculation. When set to 1, `nspin` is also automatically set to 4.
- **Default**: 0

### noncolin

- **Type**: Boolean
- **Description**: whether to allow non-collinear polarization, in which case the coupling between spin up and spin down will be taken into account. If set to 1, `nspin` is also automatically set to 4.
- **Default**: 0

### soc_lambda

- **Type**: Real
- **Description**: Relevant for soc calculations. Sometimes, for some real materials, both scalar-relativistic and full-relativistic can not describe the exact spin-orbit coupling. Artificial modulation may help in such cases.

  `soc_lambda`, which has value range [0.0, 1.0] , is used for modulate SOC effect.

  In particular, `soc_lambda 0.0` refers to scalar-relativistic case and `soc_lambda 1.0` refers to full-relativistic case.
- **Default**: 1.0

[back to top](#full-list-of-input-keywords)

## Electronic structure (SDFT)

These variables are used to control the parameters of stochastic DFT (SDFT),  mix stochastic-deterministic DFT (MDFT), or complete-basis Chebyshev method (CT). We suggest using SDFT to calculate high-temperature systems and we only support [smearing_method](#smearing_method) "fd".

### method_sto

- **Type**: Integer
- **Description**:
  - Different method to do SDFT.
  - 1: SDFT calculates $T_n(\hat{h})\ket{\chi}$ twice, where $T_n(x)$ is the n-th order Chebyshev polynomial and $\hat{h}=\frac{\hat{H}-\bar{E}}{\Delta E}$ owning eigenvalue $\in(-1,1)$. This method cost less memory but is slower.
  - 2: SDFT calculates $T_n(\hat{h})\ket{\chi}$ once but need much more memory. This method is much faster. Besides, it calculates $N_e$ with $\bra{\chi}\sqrt{\hat f}\sqrt{\hat f}\ket{\chi}$, which needs a smaller [nche_sto](#nche_sto). However, when memory is not enough, only method 1 can be used.
  - other: use 2
- **Default**: 2

### nbands_sto

- **Type**: Integer
- **Description**:
  - nbands_sto>0: Number of stochastic orbitals to calculate in SDFT and MDFT.  More bands obtain more precise results or smaller stochastic errors ($ \propto 1/\sqrt{N_{\chi}}$);
  - nbands_sto=0: Complete basis will be used to replace stochastic orbitals with the Chebyshev method (CT) and it will get the results the same as KSDFT without stochastic errors.
  - If you want to do MDFT. [nbands](#nbands) which represents the number of KS orbitals should be set.
- **Default**: 256

### nche_sto

- **Type**: Integer
- **Description**: Chebyshev expansion orders for SDFT, MDFT, CT methods.
- **Default**:100

### emin_sto

- **Type**: Real
- **Description**: Trial energy to guess the lower bound of eigen energies of the Hamitonian Operator $\hat{H}$. The unit is Ry.
- **Default**:0.0

### emax_sto

- **Type**: Real
- **Description**: Trial energy to guess the upper bound of eigen energies of the Hamitonian Operator $\hat{H}$. The unit is Ry.
- **Default**:0.0

### seed_sto

- **Type**: Integer
- **Description**: The random seed to generate stochastic orbitals.
  - seed_sto>=0: Stochastic orbitals have the form of $\exp(i2\pi\theta(G))$, where $\theta$ is a uniform distribution in $(0,1)$. If seed_sto = 0, the seed is decided by time(NULL).
  - seed_sto<=-1: Stochastic orbitals have the form of $\pm1$ with equal probability. If seed_sto = -1, the seed is decided by time(NULL).
- **Default**:0

### initsto_freq

- **Type**: Integer
- **Description**: Frequency (once each initsto_freq steps) to generate new stochastic orbitals when running md.
  - positive integer: Update stochastic orbitals
  - 0:                Never change stochastic orbitals.
- **Default**:0

### npart_sto

- **Type**: Integer
- **Description**: Make memory cost to 1/npart_sto times of the previous one when running post process of SDFT like DOS with method_sto = 2.
- **Default**:1

[back to top](#full-list-of-input-keywords)

## Geometry relaxation

These variables are used to control the geometry relaxation.

### relax_nmax

- **Type**: Integer
- **Description**: The maximal number of ionic iteration steps, the minimum value is 1.
- **Default**: 1

### cal_force

- **Type**: Boolean
- **Description**: If set to 1, calculate the force at the end of the electronic iteration. 0 means the force calculation is turned off. It is automatically set to 1 if `calculation` is `cell-relax`, `relax`, or `md`.
- **Default**: 0

### force_thr

- **Type**: Real
- **Description**: The threshold of the force convergence, it indicates the largest force among all the atoms, the unit is Ry=Bohr
- **Default**: 0.001 Ry/Bohr = 0.0257112 eV/Angstrom

### force_thr_ev

- **Type**: Real
- **Description**: The threshold of the force convergence, has the same function as force_thr, just the unit is different, it is eV/Angstrom, you can choose either one as you like. The recommended value for using atomic orbitals is 0.04 eV/Angstrom.
- **Default**: 0.0257112 eV/Angstrom

### force_thr_ev2

- **Type**: Real
- **Description**: The calculated force will be set to 0 when it is smaller than force_thr_ev2.
- **Default**: 0.0 eV/Angstrom

### relax_bfgs_w1

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into the paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.01

### relax_bfgs_w2

- **Type**: Real
- **Description**: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into the paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
- **Default**: 0.5

### relax_bfgs_rmax

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the maximal movement of all the atoms. The sum of the movements from all atoms can be increased during the optimization steps. However, it will not be larger than relax_bfgs_rmax Bohr.
- **Default**: 0.8

### relax_bfgs_rmin

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the minimal movement of all the atoms. When the movement of all the atoms is smaller than relax_bfgs_rmin Bohr, and the force convergence is still not achieved, the calculation will break down.
- **Default**: 1e-5

### relax_bfgs_init

- **Type**: Real
- **Description**: This variable is for geometry optimization. It indicates the initial movement of all the atoms. The sum of the movements from all atoms is relax_bfgs_init Bohr.
- **Default**: 0.5

### cal_stress

- **Type**: Integer
- **Description**: If set to 1, calculate the stress at the end of the electronic iteration. 0 means the stress calculation is turned off. It is automatically set to 1 if `calculation` is `cell-relax`.
- **Default**: 0

### stress_thr

- **Type**: Real
- **Description**: The threshold of the stress convergence, it indicates the largest stress among all the directions, the unit is KBar,
- **Default**: 0.01

### press1, press2, press3

- **Type**: Real
- **Description**: the external pressures along three axes, the compressive stress is taken to be positive, and the unit is KBar.
- **Default**: 0

### fixed_axes

- **Type**: String
- **Description**: which axes are fixed when do cell relaxation. Possible choices are:
  - None : default; all can relax
  - volume : relaxation with fixed volume
  - shape : fix shape but change volume (i.e. only lattice constant changes)
  - a : fix a axis during relaxation
  - b : fix b axis during relaxation
  - c : fix c axis during relaxation
  - ab : fix both a and b axes during relaxation
  - ac : fix both a and c axes during relaxation
  - bc : fix both b and c axes during relaxation

> Note : fixed_axes = "shape" and "volume" are only available for [relax_new](#relax_new) = 1

- **Default**: None

### fixed_ibrav

- **Type**: Boolean
- **Description**: when set to true, the lattice type will be preserved during relaxation. Must be used along with [relax_new](#relax_new) set to true, and a specific [latname](#latname) must be provided

> Note: it is possible to use fixed_ibrav with fixed_axes, but please make sure you know what you are doing. For example, if we are doing relaxation of a simple cubic lattice (latname = "sc"), and we use fixed_ibrav along with fixed_axes = "volume", then the cell is never allowed to move and as a result, the relaxation never converges.

- **Default**: False

### fixed_atoms

- **Type**: Boolean
- **Description**: when set to true, the direct coordinates of atoms will be preserved during variable-cell relaxation. If set to false, users can still fix certain components of certain atoms by using the `m` keyword in `STRU` file. For the latter option, check the end of this [instruction](stru.md).
- **Default**: False

### relax_method

- **Type**: String
- **Description**: The method to do geometry optimizations, note that there are two implementations of the CG method, see [relax_new](#relax_new):
  - bfgs: using BFGS algorithm.
  - sd: using steepest-descent algorithm.
  - cg: using cg algorithm.
- **Default**: cg

### relax_cg_thr

- **Type**: Real
- **Description**: When move-method is set to 'cg-bfgs', a mixed cg-bfgs algorithm is used. The ions first move according to cg method, then switched to bfgs when the maximum of force on atoms is reduced below cg-threshold. The unit is eV/Angstrom.
- **Default**: 0.5

### relax_new

- **Type**: Boolean
- **Description**: At around the end of 2022 we made a new implementation of the CG method for relax and cell-relax calculations. But the old implementation was also kept. To use the new method, set relax_new to true. To use the old one, set it to false.
- **Default**: True

### relax_scale_force

- **Type**: Real
- **Description**: This parameter is only relavant when `relax_new` is set to True. It controls the size of the first CG step. A smaller value means the first step along a new CG direction is smaller. This might be helpful for large systems, where it is safer to take a smaller initial step to prevent the collapse of the whole configuration.
- **Default**: 0.5

### cell_factor

- **Type**: Real
- **Description**: Used in the construction of the pseudopotential tables. It should exceed the maximum linear contraction of the cell during a simulation.
- **Default**: 1.2

[back to top](#full-list-of-input-keywords)

## Variables related to output information

These variables are used to control the output of properties.

### out_mul

- **Type**: Boolean
- **Description**: If set to 1, ABACUS will output the Mulliken population analysis result. The name of the output file is mulliken.txt
- **Default**: 0

### out_freq_elec

- **Type**: Integer
- **Description**: If set to >1, it represents the frequency of electronic iters to output charge density (if [out_chg](#out_chg) is turned on) and wavefunction (if [out_wfc_pw](#out_wfc_pw) or [out_wfc_r](#out_wfc_r) is turned on). If set to 0, ABACUS will output them only when converged in SCF. Used for the restart of SCF.
- **Default**: 0

### out_freq_ion

- **Type**: Integer
- **Description**: If set to >1, it represents the frequency of ionic steps to output charge density (if [out_chg](#out_chg) is turned on) and wavefunction (if [out_wfc_pw](#out_wfc_pw) or [out_wfc_r](#out_wfc_r) is turned on). If set to 0, ABACUS will output them only when ionic steps reach its maximum step. Used for the restart of MD or Relax.
- **Default**: 0

### out_chg

- **Type**: Boolean
- **Description**: If set to 1, ABACUS will output the charge density on real space grid. The name of the density file is SPIN1_CHGCAR and SPIN2_CHGCAR (if nspin = 2). Suppose each density on grid has coordinate (x; y; z). The circle order of the density on real space grid is: z is the outer loop, then y and finally x (x is moving fastest).
- **Default**: 0

### out_pot

- **Type**: Integer
- **Description**: If set to 1, ABACUS will output the local potential on real space grid. The name of the file is SPIN1_POT and SPIN2_POT (if nspin = 2). If set to 2, ABACUS will output the electrostatic potential on real space grid. The name of the file is ElecStaticPot and ElecStaticP ot_AV E (along the z-axis).

> Note : output = 1 is currently broken as of v3.0.2

- **Default**: 0

### out_dm

- **Type**: Boolean
- **Description**: If set to 1, ABACUS will output the density matrix of localized orbitals, only useful for localized orbitals set. The name of the output file is SPIN1_DM and SPIN2_DM in the output directory.
- **Default**: 0

### out_wfc_pw

- **Type**: Integer
- **Description**: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs the coefficients of wave functions into text files. The file names are WAVEFUNC$K.txt, where $K is the index of k point. When set this variable to 2, results are stored in binary files. The file names are WAVEFUNC$K.dat.
- **Default**: 0

### out_wfc_r

- **Type**: Boolean
- **Description**: Only used in **planewave basis** and **ienvelope calculation in localized orbitals** set. When set this variable to 1, it outputs real-space wave functions into  `OUT.suffix/wfc_realspace/`. The file names are wfc_realspace$K$B, where $K is the index of k point, $B is the index of band.
- **Default**: 0

### out_wfc_lcao

- **Type**: Boolean
- **Description**: **Only used in localized orbitals set**. If set to 1, ABACUS will output the wave functions coefficients.
- **Default**: 0

### out_dos

- **Type**: Integer
- **Description**: Controls whether to output the density of state (DOS). The unit of DOS is `(number of states)/(eV * unitcell)`. For more information, refer to the [worked example](../elec_properties/dos.md).
- **Default**: 0

### out_band

- **Type**: Boolean
- **Description**: Controls whether to output the band structure. For more information, refer to the [worked example](../elec_properties/band.md)
- **Default**: 0

### out_proj_band

- **Type**: Boolean
- **Description**: Controls whether to output the projected band structure. For more information, refer to the [worked example](../elec_properties/band.md)
- **Default**: 0

### out_stru

- **Type**: Boolean
- **Description**: If set to 1, then the structure files will be written after each ion step
- **Default**: 0

### out_bandgap

- **Type**: Boolean
- **Description**: If set to 1, the bandgap will be printed out at each SCF step. 
- **Default**: 0

### out_level

- **Type**: String
- **Description**: Controls the level of output. `ie` means write output at electron level; `i` means write additional output at ions level.
- **Default**: ie

### out_alllog

- **Type**: Boolean
- **Description**: determines whether to write log from all ranks in an MPI run. If set to be 1, then each rank will write detained running information to a file named running_${calculation}\_(${rank}+1).log. If set to 0, log will only be written from rank 0 into a file named running_${calculation}.log.
- **Default**: 0

### out_mat_hs

- **Type**: Boolean
- **Description**: For LCAO calculations, if out_mat_hs is set to 1, ABACUS will print the upper triangular part of the Hamiltonian matrices and overlap matrices for each k point into a series of files in the directory `OUT.${suffix}`. The files are named `data-$k-H` and `data-$k-S`, where `$k` is a composite index consisting of the k point index as well as the spin index.

  For nspin = 1 and nspin = 4 calculations, there will be only one spin component, so `$k` runs from 0 up to $nkpoints - 1$. For nspin = 2, `$k` runs from $2*nkpoints - 1$. In the latter case, the files are arranged into blocks of up and down spins. For example, if there are 3 k points, then we have the following correspondence:

  data-0-H : 1st k point, spin up
  data-1-H : 2nd k point, spin up
  data-2-H : 3rd k point, spin up
  data-3-H : 1st k point, spin down
  data-4-H : 2nd k point, spin down
  data-5-H : 3rd k point, spin down

  As for information on the k points, one may look for the `SETUP K-POINTS` section in the running log.

  The first number of the first line in each file gives the size of the matrix, namely, the number of atomic basis functions in the system.

  The rest of the file contains the upper triangular part of the specified matrices. For multi-k calculations, the matrices are Hermitian and the matrix elements are complex; for gamma-only calculations, the matrices are symmetric and the matrix elements are real.
- **Default**: 0

### out_mat_r

- **Type**: Boolean
- **Description**: For LCAO calculations, if out_mat_pos_r is set to 1, ABACUS will calculate and print the matrix representation of the position matrix, namely $\langle \chi_\mu|\hat{r}|\chi_\nu\rangle$ in a file named `data-rR-tr` in the directory `OUT.${suffix}`.

  The file starts with "Matrix Dimension of r(R): " followed by the dimension of the matrix. The rest of the format is arranged into blocks, such as:

  ```
  -5 -5 -5    //R (lattice vector)
  ...
  -5 -5 -4    //R (lattice vector)
  ...
  -5 -5 -3    //R (lattice vector)
  ```

  Each block here contains the matrix for the corresponding cell. There are three columns in each block, giving the matrix elements in x, y, z directions, respectively. There are altogether nbasis * nbasis lines in each block, which emulates the matrix elements.

  > Note: This functionality is not available for gamma_only calculations. If you want to use it in gamma_only calculations, you should turn off gamma_only, and explicitly specifies that gamma point is the only k point in the KPT file.
  >
- **Default**: 0

### out_mat_hs2

- **Type**: Boolean
- **Description**: For LCAO calculations, if out_mat_hs2 is set to 1, ABACUS will generate files containing the Hamiltonian matrix H(R) and overlap matrix S(R).

  For single-point SCF calculations, if nspin = 1 or nspin = 4, two files `data-HR-sparse_SPIN0.csr` and `data-SR-sparse_SNPIN0.csr` are generated, which contain the Hamiltonian matrix H(R) and overlap matrix S(R) respectively. For nspin = 2, three files `data-HR-sparse_SPIN0.csr` and `data-HR-sparse_SPIN1.csr` and `data-SR-sparse_SPIN0.csr` are created, where the first two contain H(R) for spin up and spin down, respectively.

  As for molecular dynamics calculations, the matrices will be printed for every `n` step, where `n` is set by another input parameter `out_hs2_interval`. Output files will be put in a separate directory, matrix_HS, and will be named `$x`_data-HR-sparse_SPIN0.csr, etc. Here `$x` is the number of MD step. For example, if we are running a 10-step MD with out_hs2_interval = 3, then `$x` will be 0,3,6 and 9.

  Each file starts with two lines, the first gives the dimension of the matrix, while the latter indicates how many different `R` are in the file.

  The rest of the files are arranged in blocks. Each block starts with a line giving the lattice vector `R` and the number of nonzero matrix elements, such as:

  ```
  -3 1 1 1020
  ```

  which means there are 1020 nonzero elements in the (-3,1,1) cell.

  If there is no nonzero matrix element, then the next block starts immediately on the next line. Otherwise, there will be 3 extra lines in the block, which gives the matrix in CSR format. According to Wikipedia:

  The CSR format stores a sparse m × n matrix M in row form using three (one-dimensional) arrays (V, COL_INDEX, ROW_INDEX). Let NNZ denote the number of nonzero entries in M. (Note that zero-based indices shall be used here.)

  - The arrays V and COL_INDEX are of length NNZ, and contain the non-zero values and the column indices of those values respectively.
  - The array ROW_INDEX is of length m + 1 and encodes the index in V and COL_INDEX where the given row starts. This is equivalent to ROW_INDEX[j] encoding the total number of nonzeros above row j. The last element is NNZ , i.e., the fictitious index in V immediately after the last valid index NNZ - 1.

  > Note: This functionality is not available for gamma_only calculations. If you want to use it in gamma_only calculations, you should turn off gamma_only, and explicitly specifies that gamma point is the only k point in the KPT file.

- **Default**: 0

### out_mat_t
- **Type**: Boolean
- **Description**: For LCAO calculations, if out_mat_t is set to 1, ABACUS will generate files containing the kinetic energy matrix T(R). The format will be the same as the Hamiltonian matrix H(R) and overlap matrix S(R) as mentioned in [out_mat_hs2](#out_mat_hs2). The name of the files will be `data-TR-sparse_SPIN0.csr` and so on. Also controled by [out_hs2_interval](#out_hs2_interval).
- **Default**: 0

### out_mat_dh
- **Type**: Boolean
- **Description**: For LCAO calculations, if out_mat_dh is set to 1, ABACUS will generate files containing the derivatives of the Hamiltonian matrix. The format will be the same as the Hamiltonian matrix H(R) and overlap matrix S(R) as mentioned in [out_mat_hs2](#out_mat_hs2). The name of the files will be `data-dHRx-sparse_SPIN0.csr` and so on. Also controled by [out_hs2_interval](#out_hs2_interval).
- **Default**: 0

### out_hs2_interval

- **Type**: Integer
- **Description**: Only relevant for printing H(R), S(R), T(R), dH(R) matrices during MD. It controls the interval at which to print. Check input parameter [out_mat_hs2](#out_mat_hs2) for more information.
- **Default**: 1

### out_element_info

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

### restart_save

- **Type**: Boolean
- **Description**: Only for LCAO, store charge density file and H matrix file every scf step for restart.
- **Default**: 0

### restart_load

- **Type**: Boolean
- **Description**: Only for LCAO, used for restart, only if that:
  - set restart_save as true and do scf calculation before.
  - please ensure the suffix is the same as calculation before and density file and H matrix file exist.
    Restart from stored density file and H matrix file.
- **Default**: 0

### dft_plus_dmft

- **Type**: Boolean
- **Description**: Whether to generate output to be used in dmft. It seems this functionality is not working anymore.
- **Default**: 0

### rpa

- **Type**: Boolean
- **Description**: Generate output files used in rpa calculation.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Density of states

These variables are used to control the calculation of DOS.

### dos_edelta_ev

- **Type**: Real
- **Description**: controls the step size in writing DOS (in eV).
- **Default**: 0.01

### dos_sigma

- **Type**: Real
- **Description**: controls the width of Gaussian factor when obtaining smeared DOS (in eV).
- **Default**: 0.07

### dos_scale

- **Type**: Real
- **Description**: the energy range of dos output is given by (emax-emin)*(1+dos_scale), centered at (emax+emin)/2. This parameter will be used when dos_emin and dos_emax are not set.
- **Default**: 0.01

### dos_emin_ev

- **Type**: Real
- **Description**: minimal range for dos (in eV). If we set it, "dos_scale" will be ignored.
- **Default**: minimal eigenenergy of $\hat{H}$

### dos_emax_ev

- **Type**: Real
- **Description**: maximal range for dos (in eV). If we set it, "dos_scale" will be ignored.
- **Default**: maximal eigenenergy of $\hat{H}$

### dos_nche

- **Type**: Integer
- **Description**: orders of Chebyshev expansions when using SDFT to calculate DOS
- **Default**: 100

[back to top](#full-list-of-input-keywords)

## NAOs

These variables are used to control the generation of numerical atomic orbitals (NAOs). NAOs is the linear combination of bessel functions.

### bessel_nao_ecut

- **Type**: Real
- **Description**: energy cutoff of bessel functions.
- **Default**: same as ecutwfc

### bessel_nao_tolerence

- **Type**: Real
- **Description**: tolerance when searching for the zeros of bessel functions.
- **Default**: 1.0e-12

### bessel_nao_rcut

- **Type**: Real
- **Description**: cutoff radius of bessel functions.
- **Default**: 6.0

### bessel_nao_smooth

- **Type**: Boolean
- **Description**: whether the bessel functions smooth at radius cutoff. 
- **Default**: 1

### bessel_nao_sigma

- **Type**: Real
- **Description**: energy range for smooth. See also `bessel_nao_smooth`.
- **Default**: 0.1

[back to top](#full-list-of-input-keywords)

## DeePKS

These variables are used to control the usage of DeePKS method (a comprehensive data-driven approach to improve the accuracy of DFT).
Warning: this function is not robust enough for the current version. Please try the following variables at your own risk:

### deepks_out_labels

- **Type**: Boolean
- **Description**: when set to 1, ABACUS will calculate and output descriptor for DeePKS training. In `LCAO` calculation, a path of *.orb file is needed to be specified under `NUMERICAL_DESCRIPTOR`in `STRU`file. For example:

  ```
  NUMERICAL_ORBITAL
  H_gga_8au_60Ry_2s1p.orb
  O_gga_7au_60Ry_2s2p1d.orb

  NUMERICAL_DESCRIPTOR
  jle.orb
  ```
- **Default**: 0

### deepks_scf

- **Type**: Boolean
- **Description**: only when deepks is enabled in `LCAO` calculation can this variable set to 1. Then, a trained, traced model file is needed for self-consistent field iteration in DeePKS method.
- **Default**: 0

### deepks_model

- **Type**: String
- **Description**: the path of the trained, traced NN model file (generated by deepks-kit). used when deepks_scf is set to 1.
- **Default**: None

### bessel_descriptor_lmax

- **Type**: Integer
- **Description**: the projectors used in DeePKS are bessel functions. To generate such projectors, set calculation type to `gen_bessel` and run ABACUS. The lmax of Bessel functions is specified using bessel_descriptor_lmax. See also [calculation](#calculation).
- **Default**: 2

### bessel_descriptor_ecut

- **Type**: Real
- **Description**: energy cutoff of bessel functions. See also `bessel_descriptor_lmax`.
- **Default**: same as ecutwfc

### bessel_descriptor_tolerence

- **Type**: Real
- **Description**: tolerance when searching for the zeros of bessel functions. See also `bessel_descriptor_lmax`.
- **Default**: 1.0e-12

### bessel_descriptor_rcut

- **Type**: Real
- **Description**: cutoff radius of bessel functions. See also `bessel_descriptor_lmax`.
- **Default**: 6.0

### bessel_descriptor_smooth

- **Type**: Boolean
- **Description**: whether the bessel functions smooth at radius cutoff. See also `bessel_descriptor_lmax`.
- **Default**: 1

### bessel_descriptor_sigma

- **Type**: Real
- **Description**: energy range for smooth. See also `bessel_descriptor_smooth`.
- **Default**: 0.1

### deepks_bandgap

- **Type**: Boolean
- **Description**: whether to include deepks bandgap correction.
- **Default**: False

### deepks_out_unittest

- **Type**: Boolean
- **Description**: this is used to generate some files for constructing DeePKS unit test. Not relevant when running actual calculations. When set to 1, ABACUS needs to be run with only 1 process.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## OFDFT: orbital free density functional theory

### of_kinetic

* **Type**: string
* **Description**: the type of kinetic energy density functional, including tf, vw, wt, and tf+.
* **Default**: wt

### of_method

- **Type**: string
- **Description**: the optimization method used in OFDFT.
  - cg1: Polak-Ribiere. Standard CG algorithm.
  - cg2: Hager-Zhang (generally faster than cg1).
  - tn: Truncated Newton algorithm.
- **Default**:tn

### of_conv

- **Type**: string
- **Description**: criterion used to check the convergence of OFDFT.
  - energy: total energy changes less than 'of_tole'.
  - potential: the norm of potential is less than 'of_tolp'.
  - both: both energy and potential must satisfy the convergence criterion.
- **Default**: energy

### of_tole

- **Type**: Double
- **Description**: tolerance of the energy change (in Ry) for determining the convergence.
- **Default**: 2e-6

### of_tolp

- **Type**: Double
- **Description**: tolerance of potential (in a.u.) for determining the convergence.
- **Default**: 1e-5

### of_tf_weight

- **Type**: Double
- **Description**: weight of TF KEDF.
- **Default**: 1

### of_vw_weight

- **Type**: Double
- **Description**: weight of vW KEDF.
- **Default**: 1

### of_wt_alpha

- **Type**: Double
- **Description**: parameter alpha of WT KEDF.
- **Default**: $5/6$

### of_wt_beta

- **Type**: Double
- **Description**: parameter beta of WT KEDF.
- **Default**: $5/6$

### of_wt_rho0

- **Type**: Double
- **Description**: the average density of system, in Bohr^-3.
- **Default**: 0

### of_hold_rho0

- **Type**: Boolean
- **Description**: If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be set to 1 automatically if of_wt_rho0 is not zero.
- **Default**: 0

### of_read_kernel

- **Type**: Boolean
- **Description**: If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from formula. Only usable for WT KEDF.
- **Default**: 0

### of_kernel_file

- **Type**: String
- **Description**: The name of WT kernel file.
- **Default**: WTkernel.txt

### of_full_pw

- **Type**: Boolean
- **Description**: If set to 1, ecut will be ignored while collecting planewaves, so that all planewaves will be used in FFT.
- **Default**: 1

### of_full_pw_dim

- **Type**: Integer
- **Description**: If of_full_pw = 1, the dimension of FFT will be restricted to be (0) either odd or even; (1) odd only; (2) even only.
  Note that even dimensions may cause slight errors in FFT. It should be ignorable in ofdft calculation, but it may make Cardinal B-**spline** interpolation unstable, so set `of_full_pw_dim = 1` if `nbspline != -1`.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electric field and dipole correction

These variables are relevant to electric field and dipole correction

### efield_flag

- **Type**: Boolean
- **Description**: If set to true, a saw-like potential simulating an electric field
  is added to the bare ionic potential.
- **Default**: false

### dip_cor_flag

- **Type**: Boolean
- **Description**: If dip_cor_flag == true and efield_flag == true,  a dipole correction is also
  added to the bare ionic potential. If you want no electric field, parameter efield_amp  should be zero. Must be used ONLY in a slab geometry for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
- **Default**: false

### efield_dir

- **Type**: Integer
- **Description**: The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir = 0, 1 or 2. Used only if efield_flag == true.
- **Default**: 2

### efield_pos_max

- **Type**: Real
- **Description**: Position of the maximum of the saw-like potential along crystal axis efield_dir, within the  unit cell, 0 < efield_pos_max < 1. Used only if efield_flag == true.
- **Default**: 0.5

### efield_pos_dec

- **Type**: Real
- **Description**: Zone in the unit cell where the saw-like potential decreases, 0 < efield_pos_dec < 1. Used only if efield_flag == true.
- **Default**: 0.1

### efield_amp

- **Type**: Real
- **Description**: Amplitude of the electric field, in ***Hartree*** a.u.; 1 a.u. = 51.4220632*10^10 V/m. Used only if efield_flag == true. The saw-like potential increases with slope efield_amp  in the region from (efield_pos_max+efield_pos_dec-1) to (efield_pos_max), then decreases until (efield_pos_max+efield_pos_dec), in units of the crystal vector efield_dir. Important: the change of slope of this potential must be located in the empty region, or else unphysical forces will result.
- **Default**: 0.0

[back to top](#full-list-of-input-keywords)

## Gate field (compensating charge)

These variables are relevant to gate field (compensating charge)

### gate_flag

- **Type**: Boolean
- **Description**: In the case of charged cells, setting gate_flag == true represents the addition of compensating charge by a charged plate, which is placed at **zgate**. Note that the direction is specified by **efield_dir**.
- **Default**: false

### zgate

- **Type**: Real
- **Description**: Specify the position of the charged plate in units of the unit cell (0 <= **zgate** < 1).
- **Default**: 0.5

### block

- **Type**: Boolean
- **Description**: Add a potential barrier to the total potential to avoid electrons spilling into the vacuum region for electron doping. Potential barrier is from **block_down** to **block_up** and has a height of **block_height**. If **dip_cor_flag** == true, **efield_pos_dec** is used for a smooth increase and decrease of the potential barrier.
- **Default**: false

### block_down

- **Type**: Real
- **Description**: Lower beginning of the potential barrier in units of the unit cell size (0 <= **block_down** < **block_up** < 1).
- **Default**: 0.45

### block_up

- **Type**: Real
- **Description**: Upper beginning of the potential barrier in units of the unit cell size (0 <= **block_down** < **block_up** < 1).
- **Default**: 0.55

### block_height

- **Type**: Real
- **Description**: Height of the potential barrier in Rydberg.
- **Default**: 0.1

[back to top](#full-list-of-input-keywords)

## Exact Exchange

These variables are relevant when using hybrid functionals

### exx_hybrid_alpha

- **Type**: Real
- **Description**: fraction of Fock exchange in hybrid functionals, so that $E_{X}=\alpha E_{X}+(1-\alpha)E_{X,\text{LDA/GGA}}$
- **Default**: 1 if dft_functional==hf else 0.25

### exx_hse_omega

- **Type**: Real
- **Description**: range-separation parameter in HSE functional, such that $1/r=\text{erfc}(\omega r)/r+\text{erf}(\omega r)/r$.
- **Default**: 0.11

### exx_separate_loop

- **Type**: Boolean
- **Description**: There are two types of iterative approaches provided by ABACUS to evaluate Fock exchange. If this parameter is set to 0, it will start with a GGA-Loop, and then Hybrid-Loop, in which EXX Hamiltonian $H_{exx}$ is updated with electronic iterations. If this parameter is set to 1, a two-step method is employed, i.e. in the inner iterations, density matrix is updated, while in the outer iterations, $H_{exx}$ is calculated based on density matrix that converges in the inner iteration. (Currently not used)
- **Default**: 1

### exx_hybrid_step

- **Type**: Integer
- **Description**: This variable indicates the maximal electronic iteration number in the evaluation of Fock exchange.
- **Default**: 100

### exx_lambda

- **Type**: Real
- **Description**: It is used to compensate for divergence points at G=0 in the evaluation of Fock exchange using *lcao_in_pw* method.
- **Default**: 0.3

### exx_pca_threshold

- **Type**: Real
- **Description**: To accelerate the evaluation of four-center integrals ($ik|jl$), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$. The size of the ABF (i.e. number of $P_{k}$) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculation becomes faster. However, this comes at the cost of computational accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). Smaller components (less than exx_c_threshold) of the $C^{k}_{ij}$ matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). With the approximation $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$, the four-center integral in Fock exchange is expressed as $(ik|jl)=\Sigma_{a,b}C^{a}_{ij}V_{ab}C^{b}_{kl}$, where $V_{ab}=(P_{a}|P_{b})$ is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_dm_threshold

- **Type**: Real
- **Description**: The Fock exchange can be expressed as $\Sigma_{k,l}(ik|jl)D_{kl}$ where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_grad_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). $\nabla C^{k}_{ij}$ is used in force and stress. Smaller components (less than exx_c_grad_threshold) of the $\nabla C^{k}_{ij}$ matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_grad_threshold

- **Type**: Real
- **Description**: See also the entry [exx_pca_threshold](#exx_pca_threshold). With the approximation $\Phi_{i}\Phi_{j}\sim C^{k}_{ij}P_{k}$, the four-center integral in Fock exchange is expressed as $(ik|jl)=\Sigma_{a,b}C^{a}_{ij}V_{ab}C^{b}_{kl}$, where $V_{ab}=(P_{a}|P_{b})$ is a double-center integral. $\nabla V_{ab}$ is used in force and stress. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_schwarz_threshold

- **Type**: Real
- **Description**: In practice the four-center integrals are sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each integral before carrying out explicit evaluations. Those that are smaller than exx_schwarz_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-5.  (Currently not used)
- **Default**: 0

### exx_cauchy_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-7.
- **Default**: 1E-7

### exx_cauchy_grad_threshold

- **Type**: Real
- **Description**: In practice the Fock exchange matrix in force and stress is sparse, and using Cauchy-Schwartz inequality, we can find an upper bound of each matrix element before carrying out explicit evaluations. Those that are smaller than exx_cauchy_grad_threshold will be truncated. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-7.
- **Default**: 1E-7

### exx_ccp_threshold

- **Type**: Real
- **Description**: It is related to the cutoff of on-site Coulomb potentials. (Currently not used)
- **Default**: 1e-8

### exx_ccp_rmesh_times

- **Type**: Real
- **Description**: This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. For HSE, setting it to 1 is enough. But for PBE0, a much larger number must be used.
- **Default**: 1.5 if dft_functional==hse else 10

### exx_distribute_type

- **Type**: String
- **Description**: When running in parallel, the evaluation of Fock exchange is done by distributing atom pairs on different threads, then gather the results. exx_distribute_type governs the mechanism of distribution. Available options are `htime`, `order`, `kmean1` and `kmeans2`. `order` is where atom pairs are simply distributed by their orders. `hmeans` is a distribution where the balance in time is achieved on each processor, hence if the memory is sufficient, this is the recommended method. `kmeans1` and `kmeans2` are two methods where the k-means clustering method is used to reduce memory requirement. They might be necessary for very large systems. (Currently not used)
- **Default**: `htime`

### exx_opt_orb_lmax

- **Type**: Integer
- **Description**: See also the entry [dft_functional](#dft_functional). This parameter is only relevant when dft_functional=`opt_orb`. The radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. exx_opt_orb_lmax gives the maximum l of the spherical Bessel functions. A reasonable choice is 2.
- **Default**: 0

### exx_opt_orb_ecut

- **Type**: Real
- **Description**: See also the entry [dft_functional](#dft_functional). This parameter is only relevant when dft_functional=`opt_orb`. A plane wave basis is used to optimize the radial ABFs. This parameter thus gives the cut-off of plane wave expansion, in Ry. A reasonable choice is 60.
- **Default**: 0

### exx_opt_orb_tolerence

- **Type**: Real
- **Description**: See also the entry [dft_functional](#dft_functional). This parameter is only relevant when dft_functional=`opt_orb`. exx_opt_orb_tolerence determines the threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.
- **Default**: 0

### exx_real_number

- **Type**: Boolen
- **Description**: If set to 1, it will enforce LIBRI to use `double` data type, otherwise, it will enforce LIBRI to use `complex` data type. The default value depends on the [gamma_only](#gamma_only) option.
- **Default**: 1 if gamma_only else 0

[back to top](#full-list-of-input-keywords)

## Molecular dynamics

These variables are used to control the molecular dynamics calculations.

### md_type

- **Type**: Integer
- **Description**: control the algorithm to integrate the equation of motion for md. When `md_type` is set to 0, `md_thermostat` is used to specify the thermostat based on the velocity Verlet algorithm.

  - -1: FIRE method to relax;
  - 0: velocity Verlet algorithm (default: NVE ensemble);
  - 1: Nose-Hoover style non-Hamiltonian equations of motion;
  - 2: NVT ensemble with Langevin method;
  - 4: MSST method;

  ***Note: when md_type is set to 1, md_tfreq is required to stablize temperature. It is an empirical parameter whose value is system-dependent, ranging from 1/(40\*md_dt) to 1/(100\*md_dt). An improper choice of its value might lead to failure of job.***
- **Default**: 1

### md_thermostat

- **Type**: String
- **Description**: specify the thermostat based on the velocity Verlet algorithm (useful when `md_type` is set to 0).

  - nve: NVE ensemble.
  - anderson: NVT ensemble with Anderson thermostat, see the parameter `md_nraise`.
  - berendsen: NVT ensemble with Berendsen thermostat, see the parameter `md_nraise`.
  - rescaling: NVT ensemble with velocity Rescaling method 1, see the parameter `md_tolerance`.
  - rescale_v: NVT ensemble with velocity Rescaling method 2, see the parameter `md_nraise`.
- **Default**: NVE

### md_nstep

- **Type**: Integer
- **Description**: the total number of md steps.
- **Default**: 10

### md_restart

- **Type**: Boolean
- **Description**: to control whether restart md.
  - 0: When set to 0, ABACUS will calculate md normally.
  - 1: When set to 1, ABACUS will calculate md from the last step in your test before.
- **Default**: 0

### md_dt

- **Type**: Real
- **Description**: This is the time step(fs) used in md simulation.
- **Default**: 1.0

### md_tfirst, md_tlast

- **Type**: Real
- **Description**: This is the temperature (K) used in md simulation. The default value of md_tlast is md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.
- **Default**: No default

### md_dumpfreq

- **Type**: Integer
- **Description**: This is the frequency to dump md information.
- **Default**: 1

### md_restartfreq

- **Type**: Integer
- **Description**: This is the frequency to output restart information.
- **Default**: 5

### md_seed

- **Type**: Integer
- **Description**:
  - md_seed < 0: No srand() in MD initialization.
  - md_seed >= 0: srand(md_seed) in MD initialization.
- **Default**: -1

### md_tfreq

- **Type**: Real
- **Description**: control the frequency of the temperature oscillations during the simulation. If it is too large, the temperature will fluctuate violently; if it is too small, the temperature will take a very long time to equilibrate with the atomic system.
- **Default**: 1/40/md_dt

### md_tchain

- **Type**: Integer
- **Description**: number of thermostats coupled with the particles in the Nose Hoover Chain method.
- **Default**: 1

### md_pmode

- **Type**: String
- **Description**: specify the NVT or NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
  - none: NVT ensemble.
  - iso: NPT ensemble with isotropic cetl fluctuations.
  - aniso: NPT ensemble with anisotropic cetl fluctuations.
  - tri: NPT ensemble with non-orthogonal (triclinic) simulation box.
- **Default**: none

### md_pcouple

- **Type**: String
- **Description**: the coupled lattice vectors will scale proportionally.
  - none: three lattice vectors scale independently.
  - xyz: lattice vectors x, y, and z scale proportionally.
  - xy: lattice vectors x and y scale proportionally.
  - xz: lattice vectors x and z scale proportionally.
  - yz: lattice vectors y and z scale proportionally.
- **Default**: none

### md_pfirst, md_plast

- **Type**: Real
- **Description**: This is the target pressure (KBar) used in npt ensemble simulation, the default value of `md_plast` is `md_pfirst`. If `md_plast` is set to be different from `md_pfirst`, ABACUS will automatically change the target pressure from `md_pfirst` to `md_plast`.
- **Default**: No default

### md_pfreq

- **Type**: Real
- **Description**: control the frequency of the pressure oscillations during the NPT ensemble simulation. If it is too large, the pressure will fluctuate violently; if it is too small, the pressure will take a very long time to equilibrate with the atomic system.
- **Default**: 1/400/md_dt

### md_pchain

- **Type**: Integer
- **Description**: number of thermostats coupled with the barostat in the Nose Hoover Chain method.
- **Default**: 1

### out_force

- **Type**: Boolean
- **Description**: Output atomic forces into the file `MD_dump` or not. If `true`, forces will be written, otherwise forces will not be written.
- **Default**: false

### out_vel

- **Type**: Boolean
- **Description**: Output atomic velocities into the file `MD_dump` or not. If `true`, velocities will be written, otherwise velocities will not be written.
- **Default**: false

### out_virial

- **Type**: Boolean
- **Description**: Output lattice virial into the file `MD_dump` or not. If `true`, lattice virial will be written, otherwise lattice virial will not be written.
- **Default**: false

### lj_rcut

- **Type**: Real
- **Description**: Cut-off radius for Leonard Jones potential (angstrom).
- **Default**: 8.5 (for He)

### lj_epsilon

- **Type**: Real
- **Description**: The value of epsilon for Leonard Jones potential (eV).
- **Default**: 0.01032 (for He)

### lj_sigma

- **Type**: Real
- **Description**: The value of sigma for Leonard Jones potential (angstrom).
- **Default**: 3.405 (for He)

### pot_file

- **Type**: String
- **Description**: The filename of potential files for CMD such as DP.
- **Default**: graph.pb

### msst_direction

- **Type**: Integer
- **Description**: the direction of shock wave for MSST.
- **Default**: 2 (z direction)

### msst_vel

- **Type**: Real
- **Description**: the velocity of shock wave (Angstrom/fs) for MSST.
- **Default**: 0.0

### msst_vis

- **Type**: Real
- **Description**: artificial viscosity (mass/length/time) for MSST.
- **Default**: 0.0

### msst_tscale

- **Type**: Real
- **Description**: reduction in initial temperature (0~1) used to compress volume in MSST.
- **Default**: 0.01

### msst_qmass

- **Type**: Real
- **Description**: Inertia of extended system variable. Used only when md_type is 4, you should set a number that is larger than 0. Note that Qmass of NHC is set by md_tfreq.
- **Default**: No default

### md_damp

- **Type**: Real
- **Description**: damping parameter (fs) used to add force in Langevin method.
- **Default**: 1.0

### md_tolerance

- **Type**: Real
- **Description**: Tolerance for velocity rescaling. Velocities are rescaled if the current and target temperature differ more than `md_tolerance` (Kelvin).
- **Default**: 100.0

### md_nraise

- **Type**: Integer
- **Description**:
  - Anderson: the "collision frequency" parameter is given as 1/`md_nraise`;
  - Berendsen: the "rise time" parameter is given in units of the time step: tau = `md_nraise`*`md_dt`, so `md_dt`/tau = 1/`md_nraise`;
  - Rescale_v: every `md_nraise` steps the current temperature is rescaled to the target temperature;
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## DFT+*U* correction

These variables are used to control DFT+U correlated parameters

### dft_plus_u

- **Type**: Boolean
- **Description**: If set to 1, ABCUS will calculate plus U correction, which is especially important for correlated electron.
- **Default**: 0

### orbital_corr

- **Type**: Integer
- **Description**: $l_1,l_2,l_3,\ldots$ for atom type 1,2,3 respectively.(usually 2 for d electrons and 3 for f electrons) .Specify which orbits need plus U correction for each atom. If set to -1, the correction would not be calculated for this atom.
- **Default**: None

### hubbard_u

- **Type**: Real
- **Description**: Hubbard Coulomb interaction parameter U(ev) in plus U correction, which should be specified for each atom unless Yukawa potential is used.
> Note : since we only implemented the simplified scheme by Duradev, the 'U' here is actually Ueff which is given by hubbard U minus hund J.
- **Default**: 0.0

### yukawa_potential

- **Type**: Boolean
- **Description**: whether to use the local screen Coulomb potential method to calculate the values of U and J. If this is set to 1, hubbard_u does not need to be specified.
- **Default**: 0

### yukawa_lambda

- **Type**: Real
- **Description**: The screen length of Yukawa potential. Relevant if `yukawa_potential` is set to 1. If left to default, we will calculate the screen length as an average of the entire system. It's better to stick to the default setting unless there is a very good reason.
- **Default**: calculated on the fly.

### omc

- **Type**: Integer
- **Description**: The parameter controls what form of occupation matrix control we are using. If set to 0, then no occupation matrix control is performed, and the onsite density matrix will be calculated from wavefunctions in each SCF step. If set to 1, then the first SCF step will use an initial density matrix read from a file named `initial_onsite.dm`, but for later steps, the onsite density matrix will be updated. If set to 2, the same onsite density matrix from `initial_onsite.dm` will be used throughout the entire calculation.
> Note : The easiest way to create `initial_onsite.dm` is to run a DFT+U calculation, look for a file named `onsite.dm` in the OUT.prefix directory, and make replacements there. The format of the file is rather straight-forward.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## vdW correction

These variables are used to control vdW-corrected related parameters.

### vdw_method

- **Type**: String
- **Description**: If set to d2 ,d3_0 or d3_bj, ABACUS will calculate corresponding vdW correction, which is DFT-D2, DFT-D3(0) or DFTD3(BJ) method. And this correction includes energy and forces. `none` means that no vdW-corrected method has been used.
- **Default**: none

### vdw_s6

- **Type**: Real
- **Description**: This scale factor is to optimize the interaction energy deviations. For DFT-D2, it is found to be 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). For DFT-D3, recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default value of this parameter in ABACUS is set to be the recommended value for PBE.
- **Default**: 0.75 if vdw_method is chosen to be d2; 1.0 if vdw_method is chosen to be d3_0 or d3_bj

### vdw_s8

- **Type**: Real
- **Description**: This scale factor is only relevant for D3(0) and D3(BJ) methods. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default value of this parameter in ABACUS is set to be the recommended value for PBE.
- **Default**: 0.722 if vdw_method is chosen to be d3_0; 0.7875 if vdw_method is chosen to be d3_bj

### vdw_a1

- **Type**: Real
- **Description**: This damping function parameter is relevant for D3(0) and D3(BJ) methods. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default value of this parameter in ABACUS is set to be the recommended value for PBE.
- **Default**: 1.217 if vdw_method is chosen to be d3_0; 0.4289 if vdw_method is chosen to be d3_bj

### vdw_a2

- **Type**: Real
- **Description**: This damping function parameter is only relevant for the DFT-D3(BJ) approach. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default value of this parameter in ABACUS is set to be the recommended value for PBE.
- **Default**: 1.0 if vdw_method is chosen to be d3_0; 4.4407 if vdw_method is chosen to be d3_bj

### vdw_d

- **Type**: Real
- **Description**: The variable is to control the damping rate of damping function of DFT-D2.
- **Default**: 20

### vdw_abc

- **Type**: Integer
- **Description**: The variable is to control whether three-body terms are calculated for DFT-D3 approaches, including D3(0) and D3(BJ). If set to 1, ABACUS will calculate three-body term, otherwise, the three-body term is not included.
- **Default**: 0

### vdw_C6_file

- **Type**: String
- **Description**: This variable is relevant if the user wants to manually set the $C_6$ parameters in D2 method. It gives the name of the file which contains the list of $C_6$ parameters for each element.

  If not set, ABACUS will use the default $C_6$ Parameters stored in the program. The default values of $C_6$ for elements 1_H up to 86_Rn can be found by searching for `C6_default` in the [source code](https://github.com/deepmodeling/abacus-develop/blob/develop/source/module_hamilt_general/module_vdw/vdwd2_parameters.cpp). The unit is Jnm6/mol.

  Otherwise, if user wants to manually set the $C_6$ Parameters, they should provide a file containing the $C_6$ parameters to be used. An example is given by:

  ```
  H  0.1
  Si 9.0
  ```

  Namely, each line contains the element name and the corresponding $C_6$ parameter.
- **Default**: default

### vdw_C6_unit

- **Type**: String
- **Description**: This variable is relevant if the user wants to manually set the $C_6$ parameters in D2 method. It specified the unit of the supplied $C_6$ parameters. Allowed values are: `Jnm6/mol` (meaning J·nm^{6}/mol) and `eVA`(meaning eV·Angstrom)
- **Default**: Jnm6/mol

### vdw_R0_file

- **Type**: String
- **Description**: This variable is relevant if the user wants to manually set the $R_0$ parameters in D2 method.
  If not set, ABACUS will use the default $C_6$ Parameters stored in the program. The default values of $C_6$ for elements 1_H up to 86_Rn can be found by searching for `R0_default` in the [source code](https://github.com/deepmodeling/abacus-develop/blob/develop/source/module_hamilt_general/module_vdw/vdwd2_parameters.cpp). The unit is Angstrom.

  Otherwise, if the user wants to manually set the $C_6$ Parameters, they should provide a file containing the $C_6$ parameters to be used. An example is given by:

  ```
  Li 1.0
  Cl 2.0
  ```

  Namely, each line contains the element name and the corresponding $R_0$ parameter.
- **Default**: default

### vdw_R0_unit

- **Type**: String
- **Description**: This variable is relevant if the user wants to manually set the $R_0$ parameters in D2 method. It specified the unit of the supplied $C_6$ parameters. Allowed values are: `A`(meaning Angstrom) and `Bohr`.
- **Default**: A

### vdw_cutoff_type

- **Type**: String
- **Description**: When applying Van-der-Waals correction in periodic systems, a cutoff radius needs to be supplied to avoid infinite  summation. In ABACUS, we restrict the range of correction to a supercell centered around the unit cell at origin.

  In ABACUS, we provide two ways to determine the extent of the supercell.

  When `vdw_cutoff_type` is set to `radius`, the supercell is chosen such that it is contained in a sphere centered at the origin. The radius of the sphere is specified by `vdw_cutoff_radius`.

  When `vdw_cutoff_type` is set to `period`, the extent of the supercell is specified explicitly using keyword `vdw_cutoff_period`.
- **Default**: radius

### vdw_cutoff_radius

- **Type**: Real
- **Description**: If `vdw_cutoff_type` is set to `radius`, this variable specifies the radius of the cutoff sphere. For DFT-D2, the default value is 56.6918, while for DFT-D3, the default value is 95.
- **Default**: 56.6918 if vdw_method is chosen to be d2; 95 if vdw_method is chosen to be d3_0 or d3_bj

### vdw_radius_unit

- **Type**: String
- **Description**: If `vdw_cutoff_type` is set to `radius`, this variable specifies the unit of `vdw_cutoff_radius`. Two values are allowed: `A`(meaning Angstrom) and `Bohr`.
- **Default**: Bohr

### vdw_cutoff_period

- **Type**: Integer Integer Integer
- **Description**: If vdw_cutoff_type is set to `period`, the three integers supplied here will explicitly specify the extent of the supercell in the directions of the three basis lattice vectors.
- **Default**: 3 3 3

### vdw_cn_thr

- **Type**: Real
- **Description**: Only relevant for D3 correction. The cutoff radius when calculating coordination numbers.
- **Default**: 40

### vdw_cn_thr_unit

- **Type**: String
- **Description**: Unit of the coordination number cutoff. Two values are allowed: `A`(meaning Angstrom) and `Bohr`.
- **Default**: Bohr

[back to top](#full-list-of-input-keywords)

## Berry phase and wannier90 interface

These variables are used to control berry phase and wannier90 interface parameters.

### berry_phase

- **Type**: Boolean
- **Description**: 1, calculate berry phase; 0, not calculate berry phase.
- **Default**: 0

### gdir

- **Type**: Integer
- **Description**:
  - 1: calculate the polarization in the direction of the lattice vector a_1 that is defined in STRU file.
  - 2: calculate the polarization in the direction of the lattice vector a_2 that is defined in STRU file.
  - 3: calculate the polarization in the direction of the lattice vector a_3 that is defined in STRU file.
- **Default**: 3

### towannier90

- **Type**: Integer
- **Description**: 1, generate files for wannier90 code; 0, no generate.
- **Default**: 0

### nnkpfile

- **Type**: String
- **Description**: the file name when you run “wannier90 -pp ...”.
- **Default**: seedname.nnkp

### wannier_spin

- **Type**: String
- **Description**: If nspin is set to 2,
  - up: calculate spin up for wannier function.
  - down: calculate spin down for wannier function.
- **Default**: up

[back to top](#full-list-of-input-keywords)

## TDDFT: time dependent density functional theory

### td_edm

- **Type**: int
- **Description**: the method to calculate the energy density matrix.
  - 0: new method (use the original formula).
  - 1: old method (use the formula for ground state).
- **Default**: 0

### td_print_eij

- **Type**: double
- **Description**: print the Eij(<\psi_i|H|\psi_j>) elements which are larger than td_print_eij. if td_print_eij <0, don't print Eij 
- **Default**: -1

### td_force_dt

- **Type**: Real
- **Description**: Time-dependent evolution force changes time step. (fs)
- **Default**: 0.02

### td_vext

- **Type**: Boolean
- **Description**:
  - 1: add a laser material interaction (extern laser field).
  - 0: no extern laser field.
- **Default**: 0

### td_vext_dire

- **Type**: String
- **Description**:
  If `td_vext` is true, the td_vext_dire is a string to set the number of electric field, like "1 2" representing external electric field is added to the x and y axis at the same time. Parameters of electric field can also be written as a string, like `td_gauss_phase 0 1.5707963267948966` representing the Gauss field in the x and y directions has a phase delay of Pi/2. See below for more parameters of electric field.
  - 1: the direction of external light field is along x axis.
  - 2: the direction of external light field is along y axis.
  - 3: the direction of external light field is along z axis.
- **Default**: 1

### td_stype

- **Type**: Integer
- **Description**:
  type of electric field in space domain
  - 0: length gauge.
  - 1: velocity gauge.
- **Default**: 0

### td_ttype

- **Type**: String
- **Description**:
  type of electric field in time domain
  - 0: Gaussian type function.
  - 1: Trapezoid function.
  - 2: Trigonometric function.
  - 3: Heaviside function.
  - 4: HHG function.
- **Default**: 0

### td_tstart

- **Type**: Integer
- **Description**:
  nubmer of step where electric field start
- **Default**: 1

### td_tend

- **Type**: Integer
- **Description**:
  nubmer of step where electric field end
- **Default**: 100

### td_lcut1

- **Type**: Double
- **Description**:
  cut1 of interval in length gauge
  E = E0 , cut1<x<cut2
  E = -E0/(cut1+1-cut2) , x<cut1 or cut2<x<1
- **Default**: 0.05

### td_lcut2

- **Type**: Double
- **Description**:
  cut2 of interval in length gauge
- **Default**: 0.05

### td_gauss_freq

- **Type**: String
- **Description**:
  frequency of Gauss type elctric field  (fs^-1)
  amp*cos(2pi*f(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 22.13

### td_gauss_phase

- **Type**: String
- **Description**:
  phase of Gauss type elctric field  
  amp*cos(2pi*f(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 0.0

### td_gauss_sigma

- **Type**: String
- **Description**:
  sigma of Gauss type elctric field  (fs)
  amp*cos(2pi*f(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 30.0

### td_gauss_t0

- **Type**: String
- **Description**:
  step number of time center of Gauss type elctric field  
  amp*cos(2pi*f(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 100

### td_gauss_amp

- **Type**: String
- **Description**:
  amplitude of Gauss type elctric field  (V/A)
  amp*cos(2pi*f(t-t0)+phase)exp(-(t-t0)^2/2sigma^2)
- **Default**: 0.25

### td_trape_freq

- **Type**: String
- **Description**:
  frequency of Trapezoid type elctric field  (fs^-1)
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 1.60

### td_trape_phase

- **Type**: String
- **Description**:
  phase of Trapezoid type elctric field  
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 0.0

### td_trape_t1

- **Type**: String
- **Description**:
  step number of time interval 1 of Trapezoid type elctric field  
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 1875

### td_trape_t2

- **Type**: String
- **Description**:
  step number of time interval 2 of Trapezoid type elctric field  
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 5625

### td_trape_t3

- **Type**: String
- **Description**:
  step number of time interval 3 of Trapezoid type elctric field  
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 7500

### td_trape_amp

- **Type**: String
- **Description**:
  amplitude of Trapezoid type elctric field  (V/A)
  E = amp*cos(2pi*f*t+phase) t/t1 , t<t1
  E = amp*cos(2pi*f*t+phase) , t1<t<t2
  E = amp*cos(2pi*f*t+phase) (1-(t-t2)/(t3-t2)) , t2<t<t3
  E = 0 , t>t3
- **Default**: 2.74

### td_trigo_freq1

- **Type**: String
- **Description**:
  frequence 1 of Trigonometric type elctric field  (fs^-1)
  amp*cos(2*pi*f1*t+phase1)*sin(2*pi*f2*t+phase2)^2
- **Default**: 1.164656

### td_trigo_freq2

- **Type**: String
- **Description**:
  frequence 2 of Trigonometric type elctric field  (fs^-1)
  amp*cos(2*pi*f1*t+phase1)*sin(2*pi*f2*t+phase2)^2
- **Default**: 0.029116

### td_trigo_phase1

- **Type**: String
- **Description**:
  phase 1 of Trigonometric type elctric field  
  amp*cos(2*pi*f1*t+phase1)*sin(2*pi*f2*t+phase2)^2
- **Default**: 0.0

### td_trigo_phase2

- **Type**: String
- **Description**:
  phase 2 of Trigonometric type elctric field  
  amp*cos(2*pi*f1*t+phase1)*sin(2*pi*f2*t+phase2)^2
- **Default**: 0.0

### td_trigo_amp

- **Type**: String
- **Description**:
  amplitude of Trigonometric type elctric field  (V/A)
  amp*cos(2*pi*f1*t+phase1)*sin(2*pi*f2*t+phase2)^2
- **Default**: 2.74

### td_heavi_t0

- **Type**: String
- **Description**:
  step number of switch time of Heaviside type elctric field 
  E = amp , t<t0
  E = 0.0 , t>t0
- **Default**: 100
### td_heavi_amp

- **Type**: String
- **Description**:
  amplitude of Heaviside type elctric field  (V/A)
  E = amp , t<t0
  E = 0.0 , t>t0
- **Default**: 2.74

### out_dipole

- **Type**: Integer
- **Description**:
  - 1: Output dipole.
  - 0: do not Output dipole.
- **Default**: 0

### out_efield

- **Type**: Integer
- **Description**:
  - 1: Output efield.
  - 0: do not Output efield.
- **Default**: 0


### ocp

- **Type**: Boolean
- **Description**: choose whether calculating constrained DFT or not.
  - For PW and LCAO codes. if set to 1, occupations of bands will be setting of "ocp_set".
  - For TDDFT in LCAO codes. if set to 1, occupations will be constrained since second ionic step.
  - For OFDFT, this feature can't be used.
- **Default**:0

### ocp_set

- **Type**: string
- **Description**: If ocp is true, the ocp_set is a string to set the number of occupancy, like 1 10 * 1 0 1 representing the 13 band occupancy, 12th band occupancy 0 and the rest 1, the code is parsing this string into an array through a regular expression.
- **Default**: none

### td_val_elec_01

- **Type**: Integer
- **Description**: Only useful when calculating the dipole. Specifies the number of valence electron associated with the first element.
- **Default**: 1.0

### td_val_elec_02

- **Type**: Integer
- **Description**: Only useful when calculating the dipole. Specifies the number of valence electron associated with the second element.
- **Default**: 1.0

### td_val_elec_03

- **Type**: Integer
- **Description**: Only useful when calculating the dipole. Specifies the number of valence electron associated with the third element.
- **Default**: 1.0

[back to top](#full-list-of-input-keywords)

## Variables useful for debugging

### nurse

- **Type**: Boolean
- **Description**: If set to 1, the Hamiltonian matrix and S matrix in each iteration will be written in output.
- **Default**: 0

### t_in_h

- **Type**: Boolean
- **Description**: If set to 0, then kinetic term will not be included in obtaining the Hamiltonian.
- **Default**: 1

### vl_in_h

- **Type**: Boolean
- **Description**: If set to 0, then local pseudopotential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

### vnl_in_h

- **Type**: Boolean
- **Description**:  If set to 0, then non-local pseudopotential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

### vh_in_h

- **Type**: Boolean
- **Description**:  If set to 0, then Hartree potential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

### vion_in_h

- **Type**: Boolean
- **Description**:  If set to 0, then local ionic potential term will not be included in obtaining the Hamiltonian.
- **Default**: 1

### test_force

- **Type**: Boolean
- **Description**: If set to 1, then detailed components in forces will be written to output.
- **Default**: 0

### test_stress

- **Type**: Boolean
- **Description**: If set to 1, then detailed components in stress will be written to output.
- **Default**: 0

### colour

- **Type**: Boolean
- **Description**: If set to 1, output to terminal will have some color.
- **Default**: 0

### test_skip_ewald

- **Type**: Boolean
- **Description**: If set to 1, then ewald energy will not be calculated.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electronic conductivities

Frequency-dependent electronic conductivities can be calculated with Kubo-Greenwood formula[Phys. Rev. B 83, 235120 (2011)].

Onsager coefficiencies:

$L_{mn}(\omega)=(-1)^{m+n}\frac{2\pi e^2\hbar^2}{3m_e^2\omega\Omega}$

$\times\sum_{ij\alpha\mathbf{k}}W(\mathbf{k})\left(\frac{\epsilon_{i\mathbf{k}}+\epsilon_{j\mathbf{k}}}{2}-\mu\right)^{m+n-2} \times |\langle\Psi_{i\mathbf{k}}|\nabla_\alpha|\Psi_{j\mathbf{k}}\rangle|^2$

$\times[f(\epsilon_{i\mathbf{k}})-f(\epsilon_{j\mathbf{k}})]\delta(\epsilon_{j\mathbf{k}}-\epsilon_{i\mathbf{k}}-\hbar\omega).$

They can also be computed by $j$-$j$ correlation function.

$L_{mn}=\frac{2e^{m+n-2}}{3\Omega\hbar\omega}\Im[\tilde{C}_{mn}(\omega)]$

$\tilde{C}_{mn}=\int_0^\infty C_{mn}(t)e^{-i\omega t}e^{-\frac{1}{2}(\Delta E)^2t^2}dt$

$C_{mn}(t)=-2\theta(t)\Im\left\{Tr\left[\sqrt{\hat f}\hat{j}_m(1-\hat{f})e^{i\frac{\hat{H}}{\hbar}t}\hat{j}_ne^{-i\frac{\hat{H}}{\hbar}t}\sqrt{\hat f}\right]\right\}$,

where $j_1$ is electric flux and $j_2$ is thermal flux.

Frequency-dependent electric conductivities: $\sigma(\omega)=L_{11}(\omega)$.

Frequency-dependent thermal conductivities: $\kappa(\omega)=\frac{1}{e^2T}\left(L_{22}-\frac{L_{12}^2}{L_{11}}\right)$.

DC electric conductivities: $\sigma = \lim_{\omega\to 0}\sigma(\omega)$.

Thermal conductivities: $\kappa = \lim_{\omega\to 0}\kappa(\omega)$.

### cal_cond

- **Type**: Boolean
- **Description**: If set to 1, electronic conductivities will be calculated. Only supported in calculations of SDFT and KSDFT_PW.
- **Default**: 0

### cond_nche

- **Type**: Integer
- **Description**: Chebyshev expansion orders for stochastic Kubo Greenwood. Only used when the calculation is SDFT.
- **Default**: 20

### cond_dw

- **Type**: Real
- **Description**: Frequency interval ($d\omega$) for frequency-dependent conductivities. The unit is eV.
- **Default**: 0.1

### cond_wcut

- **Type**: Real
- **Description**: Cutoff frequency for frequency-dependent conductivities. The unit is eV.
- **Default**: 10.0

### cond_wenlarge

- **Type**: Integer
- **Description**: Control the t interval: dt = $\frac{\pi}{\omega_{cut}\times\omega enlarge}$
- **Default**: 10

### cond_fwhm

- **Type**: Integer
- **Description**: We use gaussian functions to approximate $\delta(E)\approx \frac{1}{\sqrt{2\pi}\Delta E}e^{-\frac{E^2}{2{\Delta E}^2}}$. FWHM for conductivities, $FWHM=2*\sqrt{2\ln2}\cdot \Delta E$. The unit is eV.
- **Default**: 0.3

### cond_nonlocal

- **Type**: Boolean
- **Description**: Conductivities need to calculate velocity matrix $\bra{\psi_i}\hat{v}\ket{\psi_j}$ and $m\hat{v}=\hat{p}+\frac{im}{\hbar}[\hat{V}_{NL},\hat{r}]$. If `cond_nonlocal` is false, $m\hat{v}\approx\hat{p}$.
- **Default**: True

[back to top](#full-list-of-input-keywords)

## Implicit solvation model

These variables are used to control the usage of implicit solvation model. This approach treats the solvent as a continuous medium instead of individual “explicit” solvent molecules, which means that the solute embedded in an implicit solvent and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath.

### imp_sol

- **Type**: Boolean
- **Description**: If set to 1, an implicit solvation correction is considered.
- **Default**: 0

### eb_k

- **Type**: Real
- **Description**: The relative permittivity of the bulk solvent, 80 for water. Used only if `imp_sol` == true.
- **Default**: 80

### tau

- **Type**: Real
- **Description**: The effective surface tension parameter, which describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent that are not captured by the electrostatic terms. The unit is $Ry/Bohr^{2}$.
- **Default**: 1.0798e-05

### sigma_k

- **Type**: Real
- **Description**: We assume a diffuse cavity that is implicitly determined by the electronic structure of the solute.
  `sigma_k` is the parameter that describes the width of the diffuse cavity.
- **Default**: 0.6

### nc_k

- **Type**: Real
- **Description**: It determines at what value of the electron density the dielectric cavity forms.
  The unit is $Bohr^{-3}$.
- **Default**: 0.00037

[back to top](#full-list-of-input-keywords)

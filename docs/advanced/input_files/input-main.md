# Full List of INPUT Keywords

<!-- This file is auto-generated from parameters.yaml -->
<!-- Do not edit manually - changes will be overwritten -->

<!-- Table of Contents -->
- [Full List of INPUT Keywords](#full-list-of-input-keywords)
  - [System variables](#system-variables)
    - [suffix](#suffix)
    - [ntype](#ntype)
    - [calculation](#calculation)
    - [esolver\_type](#esolver_type)
    - [symmetry](#symmetry)
    - [symmetry\_prec](#symmetry_prec)
    - [symmetry\_autoclose](#symmetry_autoclose)
    - [cal\_force](#cal_force)
    - [kpar](#kpar)
    - [bndpar](#bndpar)
    - [latname](#latname)
    - [init\_wfc](#init_wfc)
    - [init\_chg](#init_chg)
    - [init\_vel](#init_vel)
    - [mem\_saver](#mem_saver)
    - [cal\_stress](#cal_stress)
    - [diago\_proc](#diago_proc)
    - [nbspline](#nbspline)
    - [kspacing](#kspacing)
    - [min\_dist\_coef](#min_dist_coef)
    - [device](#device)
    - [precision](#precision)
    - [timer\_enable\_nvtx](#timer_enable_nvtx)
    - [cell\_factor](#cell_factor)
    - [dm\_to\_rho](#dm_to_rho)
    - [chg\_extrap](#chg_extrap)
    - [nb2d](#nb2d)
    - [cal\_symm\_repr](#cal_symm_repr)
  - [Input files](#input-files)
    - [stru\_file](#stru_file)
    - [kpoint\_file](#kpoint_file)
    - [pseudo\_dir](#pseudo_dir)
    - [orbital\_dir](#orbital_dir)
    - [read\_file\_dir](#read_file_dir)
    - [restart\_load](#restart_load)
    - [spillage\_outdir](#spillage_outdir)
  - [Plane wave related variables](#plane-wave-related-variables)
    - [ecutwfc](#ecutwfc)
    - [ecutrho](#ecutrho)
    - [nx](#nx)
    - [ny](#ny)
    - [nz](#nz)
    - [ndx](#ndx)
    - [ndy](#ndy)
    - [ndz](#ndz)
    - [pw\_seed](#pw_seed)
    - [diag\_subspace](#diag_subspace)
    - [erf\_ecut](#erf_ecut)
    - [fft\_mode](#fft_mode)
    - [erf\_height](#erf_height)
    - [erf\_sigma](#erf_sigma)
    - [pw\_diag\_thr](#pw_diag_thr)
    - [diago\_smooth\_ethr](#diago_smooth_ethr)
    - [use\_k\_continuity](#use_k_continuity)
    - [pw\_diag\_nmax](#pw_diag_nmax)
    - [pw\_diag\_ndim](#pw_diag_ndim)
    - [diago\_cg\_prec](#diago_cg_prec)
  - [Numerical atomic orbitals related variables](#numerical-atomic-orbitals-related-variables)
    - [lmaxmax](#lmaxmax)
    - [lcao\_ecut](#lcao_ecut)
    - [lcao\_dk](#lcao_dk)
    - [lcao\_dr](#lcao_dr)
    - [lcao\_rmax](#lcao_rmax)
    - [search\_radius](#search_radius)
    - [bx](#bx)
    - [by](#by)
    - [bz](#bz)
    - [elpa\_num\_thread](#elpa_num_thread)
    - [num\_stream](#num_stream)
  - [Electronic structure](#electronic-structure)
    - [basis\_type](#basis_type)
    - [ks\_solver](#ks_solver)
    - [nbands](#nbands)
    - [nelec](#nelec)
    - [nelec\_delta](#nelec_delta)
    - [nupdown](#nupdown)
    - [dft\_functional](#dft_functional)
    - [xc\_temperature](#xc_temperature)
    - [xc\_exch\_ext](#xc_exch_ext)
    - [xc\_corr\_ext](#xc_corr_ext)
    - [pseudo\_rcut](#pseudo_rcut)
    - [pseudo\_mesh](#pseudo_mesh)
    - [nspin](#nspin)
    - [smearing\_method](#smearing_method)
    - [smearing\_sigma](#smearing_sigma)
    - [smearing\_sigma\_temp](#smearing_sigma_temp)
    - [mixing\_type](#mixing_type)
    - [mixing\_beta](#mixing_beta)
    - [mixing\_beta\_mag](#mixing_beta_mag)
    - [mixing\_ndim](#mixing_ndim)
    - [mixing\_restart](#mixing_restart)
    - [mixing\_dmr](#mixing_dmr)
    - [mixing\_gg0](#mixing_gg0)
    - [mixing\_gg0\_mag](#mixing_gg0_mag)
    - [mixing\_gg0\_min](#mixing_gg0_min)
    - [mixing\_angle](#mixing_angle)
    - [mixing\_tau](#mixing_tau)
    - [mixing\_dftu](#mixing_dftu)
    - [gamma\_only](#gamma_only)
    - [scf\_nmax](#scf_nmax)
    - [scf\_thr](#scf_thr)
    - [scf\_ene\_thr](#scf_ene_thr)
    - [scf\_thr\_type](#scf_thr_type)
    - [scf\_os\_stop](#scf_os_stop)
    - [scf\_os\_thr](#scf_os_thr)
    - [scf\_os\_ndim](#scf_os_ndim)
    - [sc\_os\_ndim](#sc_os_ndim)
    - [lspinorb](#lspinorb)
    - [noncolin](#noncolin)
    - [soc\_lambda](#soc_lambda)
    - [dfthalf\_type](#dfthalf_type)
  - [Electronic structure (SDFT)](#electronic-structure-sdft)
    - [method\_sto](#method_sto)
    - [nbands\_sto](#nbands_sto)
    - [nche\_sto](#nche_sto)
    - [emin\_sto](#emin_sto)
    - [emax\_sto](#emax_sto)
    - [seed\_sto](#seed_sto)
    - [initsto\_ecut](#initsto_ecut)
    - [initsto\_freq](#initsto_freq)
    - [npart\_sto](#npart_sto)
  - [Geometry relaxation](#geometry-relaxation)
    - [relax\_method](#relax_method)
    - [relax\_new](#relax_new)
    - [relax\_scale\_force](#relax_scale_force)
    - [relax\_nmax](#relax_nmax)
    - [relax\_cg\_thr](#relax_cg_thr)
    - [force\_thr](#force_thr)
    - [force\_thr\_ev](#force_thr_ev)
    - [force\_zero\_out](#force_zero_out)
    - [relax\_bfgs\_w1](#relax_bfgs_w1)
    - [relax\_bfgs\_w2](#relax_bfgs_w2)
    - [relax\_bfgs\_rmax](#relax_bfgs_rmax)
    - [relax\_bfgs\_rmin](#relax_bfgs_rmin)
    - [relax\_bfgs\_init](#relax_bfgs_init)
    - [stress\_thr](#stress_thr)
    - [press1](#press1)
    - [press2](#press2)
    - [press3](#press3)
    - [fixed\_axes](#fixed_axes)
    - [fixed\_ibrav](#fixed_ibrav)
    - [fixed\_atoms](#fixed_atoms)
  - [Output information](#output-information)
    - [out\_freq\_ion](#out_freq_ion)
    - [out\_freq\_td](#out_freq_td)
    - [out\_freq\_elec](#out_freq_elec)
    - [out\_chg](#out_chg)
    - [out\_pot](#out_pot)
    - [out\_dmk](#out_dmk)
    - [out\_dmr](#out_dmr)
    - [out\_wfc\_pw](#out_wfc_pw)
    - [out\_wfc\_lcao](#out_wfc_lcao)
    - [out\_dos](#out_dos)
    - [out\_ldos](#out_ldos)
    - [out\_band](#out_band)
    - [out\_proj\_band](#out_proj_band)
    - [out\_stru](#out_stru)
    - [out\_level](#out_level)
    - [out\_mat\_hs](#out_mat_hs)
    - [out\_mat\_hs2](#out_mat_hs2)
    - [out\_mat\_tk](#out_mat_tk)
    - [out\_mat\_r](#out_mat_r)
    - [out\_mat\_t](#out_mat_t)
    - [out\_mat\_dh](#out_mat_dh)
    - [out\_mat\_ds](#out_mat_ds)
    - [out\_mat\_xc](#out_mat_xc)
    - [out\_mat\_xc2](#out_mat_xc2)
    - [out\_mat\_l](#out_mat_l)
    - [out\_xc\_r](#out_xc_r)
    - [out\_eband\_terms](#out_eband_terms)
    - [out\_mul](#out_mul)
    - [out\_app\_flag](#out_app_flag)
    - [out\_ndigits](#out_ndigits)
    - [out\_element\_info](#out_element_info)
    - [restart\_save](#restart_save)
    - [rpa](#rpa)
    - [out\_pchg](#out_pchg)
    - [out\_wfc\_norm](#out_wfc_norm)
    - [out\_wfc\_re\_im](#out_wfc_re_im)
    - [if\_separate\_k](#if_separate_k)
    - [out\_elf](#out_elf)
    - [out\_spillage](#out_spillage)
    - [out\_alllog](#out_alllog)
  - [Density of states](#density-of-states)
    - [dos\_edelta\_ev](#dos_edelta_ev)
    - [dos\_sigma](#dos_sigma)
    - [dos\_scale](#dos_scale)
    - [dos\_emin\_ev](#dos_emin_ev)
    - [dos\_emax\_ev](#dos_emax_ev)
    - [dos\_nche](#dos_nche)
    - [stm\_bias](#stm_bias)
    - [ldos\_line](#ldos_line)
  - [NAOs](#naos)
    - [bessel\_nao\_ecut](#bessel_nao_ecut)
    - [bessel\_nao\_tolerence](#bessel_nao_tolerence)
    - [bessel\_nao\_rcut](#bessel_nao_rcut)
    - [bessel\_nao\_smooth](#bessel_nao_smooth)
    - [bessel\_nao\_sigma](#bessel_nao_sigma)
  - [DeePKS](#deepks)
    - [deepks\_out\_labels](#deepks_out_labels)
    - [deepks\_out\_freq\_elec](#deepks_out_freq_elec)
    - [deepks\_out\_base](#deepks_out_base)
    - [deepks\_scf](#deepks_scf)
    - [deepks\_equiv](#deepks_equiv)
    - [deepks\_model](#deepks_model)
    - [bessel\_descriptor\_lmax](#bessel_descriptor_lmax)
    - [bessel\_descriptor\_ecut](#bessel_descriptor_ecut)
    - [bessel\_descriptor\_tolerence](#bessel_descriptor_tolerence)
    - [bessel\_descriptor\_rcut](#bessel_descriptor_rcut)
    - [bessel\_descriptor\_smooth](#bessel_descriptor_smooth)
    - [bessel\_descriptor\_sigma](#bessel_descriptor_sigma)
    - [deepks\_bandgap](#deepks_bandgap)
    - [deepks\_band\_range](#deepks_band_range)
    - [deepks\_v\_delta](#deepks_v_delta)
    - [deepks\_out\_unittest](#deepks_out_unittest)
  - [OFDFT: orbital free density functional theory](#ofdft-orbital-free-density-functional-theory)
    - [of\_kinetic](#of_kinetic)
    - [of\_method](#of_method)
    - [of\_conv](#of_conv)
    - [of\_tole](#of_tole)
    - [of\_tolp](#of_tolp)
    - [of\_tf\_weight](#of_tf_weight)
    - [of\_vw\_weight](#of_vw_weight)
    - [of\_wt\_alpha](#of_wt_alpha)
    - [of\_wt\_beta](#of_wt_beta)
    - [of\_wt\_rho0](#of_wt_rho0)
    - [of\_hold\_rho0](#of_hold_rho0)
    - [of\_lkt\_a](#of_lkt_a)
    - [of\_xwm\_rho\_ref](#of_xwm_rho_ref)
    - [of\_xwm\_kappa](#of_xwm_kappa)
    - [of\_read\_kernel](#of_read_kernel)
    - [of\_kernel\_file](#of_kernel_file)
    - [of\_full\_pw](#of_full_pw)
    - [of\_full\_pw\_dim](#of_full_pw_dim)
  - [ML-KEDF: machine learning based kinetic energy density functional for OFDFT](#ml-kedf-machine-learning-based-kinetic-energy-density-functional-for-ofdft)
    - [of\_ml\_gene\_data](#of_ml_gene_data)
    - [of\_ml\_device](#of_ml_device)
    - [of\_ml\_feg](#of_ml_feg)
    - [of\_ml\_nkernel](#of_ml_nkernel)
    - [of\_ml\_kernel](#of_ml_kernel)
    - [of\_ml\_kernel\_scaling](#of_ml_kernel_scaling)
    - [of\_ml\_yukawa\_alpha](#of_ml_yukawa_alpha)
    - [of\_ml\_kernel\_file](#of_ml_kernel_file)
    - [of\_ml\_gamma](#of_ml_gamma)
    - [of\_ml\_p](#of_ml_p)
    - [of\_ml\_q](#of_ml_q)
    - [of\_ml\_tanhp](#of_ml_tanhp)
    - [of\_ml\_tanhq](#of_ml_tanhq)
    - [of\_ml\_chi\_p](#of_ml_chi_p)
    - [of\_ml\_chi\_q](#of_ml_chi_q)
    - [of\_ml\_gammanl](#of_ml_gammanl)
    - [of\_ml\_pnl](#of_ml_pnl)
    - [of\_ml\_qnl](#of_ml_qnl)
    - [of\_ml\_xi](#of_ml_xi)
    - [of\_ml\_tanhxi](#of_ml_tanhxi)
    - [of\_ml\_tanhxi\_nl](#of_ml_tanhxi_nl)
    - [of\_ml\_tanh\_pnl](#of_ml_tanh_pnl)
    - [of\_ml\_tanh\_qnl](#of_ml_tanh_qnl)
    - [of\_ml\_tanhp\_nl](#of_ml_tanhp_nl)
    - [of\_ml\_tanhq\_nl](#of_ml_tanhq_nl)
    - [of\_ml\_chi\_xi](#of_ml_chi_xi)
    - [of\_ml\_chi\_pnl](#of_ml_chi_pnl)
    - [of\_ml\_chi\_qnl](#of_ml_chi_qnl)
    - [of\_ml\_local\_test](#of_ml_local_test)
    - [ml\_exx](#ml_exx)
  - [TDOFDFT: time dependent orbital free density functional theory](#tdofdft-time-dependent-orbital-free-density-functional-theory)
    - [of\_cd](#of_cd)
    - [of\_mcd\_alpha](#of_mcd_alpha)
  - [Electric field and dipole correction](#electric-field-and-dipole-correction)
    - [efield\_flag](#efield_flag)
    - [dip\_cor\_flag](#dip_cor_flag)
    - [efield\_dir](#efield_dir)
    - [efield\_pos\_max](#efield_pos_max)
    - [efield\_pos\_dec](#efield_pos_dec)
    - [efield\_amp](#efield_amp)
  - [Gate field (compensating charge)](#gate-field-compensating-charge)
    - [gate\_flag](#gate_flag)
    - [zgate](#zgate)
    - [block](#block)
    - [block\_down](#block_down)
    - [block\_up](#block_up)
    - [block\_height](#block_height)
  - [Exact Exchange (Common)](#exact-exchange-common)
    - [exx\_fock\_alpha](#exx_fock_alpha)
    - [exx\_erfc\_alpha](#exx_erfc_alpha)
    - [exx\_erfc\_omega](#exx_erfc_omega)
    - [exx\_separate\_loop](#exx_separate_loop)
    - [exx\_hybrid\_step](#exx_hybrid_step)
    - [exx\_mixing\_beta](#exx_mixing_beta)
  - [Exact Exchange (LCAO in PW)](#exact-exchange-lcao-in-pw)
    - [exx\_fock\_lambda](#exx_fock_lambda)
  - [Exact Exchange (LCAO)](#exact-exchange-lcao)
    - [exx\_pca\_threshold](#exx_pca_threshold)
    - [exx\_c\_threshold](#exx_c_threshold)
    - [exx\_cs\_inv\_thr](#exx_cs_inv_thr)
    - [exx\_v\_threshold](#exx_v_threshold)
    - [exx\_dm\_threshold](#exx_dm_threshold)
    - [exx\_c\_grad\_threshold](#exx_c_grad_threshold)
    - [exx\_v\_grad\_threshold](#exx_v_grad_threshold)
    - [exx\_c\_grad\_r\_threshold](#exx_c_grad_r_threshold)
    - [exx\_v\_grad\_r\_threshold](#exx_v_grad_r_threshold)
    - [exx\_ccp\_rmesh\_times](#exx_ccp_rmesh_times)
    - [exx\_opt\_orb\_lmax](#exx_opt_orb_lmax)
    - [exx\_opt\_orb\_ecut](#exx_opt_orb_ecut)
    - [exx\_opt\_orb\_tolerence](#exx_opt_orb_tolerence)
    - [exx\_real\_number](#exx_real_number)
    - [exx\_singularity\_correction](#exx_singularity_correction)
    - [rpa\_ccp\_rmesh\_times](#rpa_ccp_rmesh_times)
    - [exx\_symmetry\_realspace](#exx_symmetry_realspace)
    - [out\_ri\_cv](#out_ri_cv)
  - [Exact Exchange (PW)](#exact-exchange-pw)
    - [exxace](#exxace)
    - [exx\_gamma\_extrapolation](#exx_gamma_extrapolation)
    - [ecutexx](#ecutexx)
    - [exx\_thr\_type](#exx_thr_type)
    - [exx\_ene\_thr](#exx_ene_thr)
  - [Molecular dynamics](#molecular-dynamics)
    - [md\_type](#md_type)
    - [md\_nstep](#md_nstep)
    - [md\_dt](#md_dt)
    - [md\_thermostat](#md_thermostat)
    - [md\_tfirst](#md_tfirst)
    - [md\_tlast](#md_tlast)
    - [md\_prec\_level](#md_prec_level)
    - [md\_restart](#md_restart)
    - [md\_restartfreq](#md_restartfreq)
    - [md\_dumpfreq](#md_dumpfreq)
    - [dump\_force](#dump_force)
    - [dump\_vel](#dump_vel)
    - [dump\_virial](#dump_virial)
    - [md\_seed](#md_seed)
    - [md\_tfreq](#md_tfreq)
    - [md\_tchain](#md_tchain)
    - [md\_pmode](#md_pmode)
    - [ref\_cell\_factor](#ref_cell_factor)
    - [md\_pcouple](#md_pcouple)
    - [md\_pfirst](#md_pfirst)
    - [md\_plast](#md_plast)
    - [md\_pfreq](#md_pfreq)
    - [md\_pchain](#md_pchain)
    - [lj\_rule](#lj_rule)
    - [lj\_eshift](#lj_eshift)
    - [lj\_rcut](#lj_rcut)
    - [lj\_epsilon](#lj_epsilon)
    - [lj\_sigma](#lj_sigma)
    - [pot\_file](#pot_file)
    - [dp\_rescaling](#dp_rescaling)
    - [dp\_fparam](#dp_fparam)
    - [dp\_aparam](#dp_aparam)
    - [msst\_direction](#msst_direction)
    - [msst\_vel](#msst_vel)
    - [msst\_vis](#msst_vis)
    - [msst\_tscale](#msst_tscale)
    - [msst\_qmass](#msst_qmass)
    - [md\_damp](#md_damp)
    - [md\_tolerance](#md_tolerance)
    - [md\_nraise](#md_nraise)
    - [cal\_syns](#cal_syns)
    - [dmax](#dmax)
  - [DFT+U correction](#dftu-correction)
    - [dft\_plus\_u](#dft_plus_u)
    - [dft\_plus\_dmft](#dft_plus_dmft)
    - [orbital\_corr](#orbital_corr)
    - [hubbard\_u](#hubbard_u)
    - [yukawa\_potential](#yukawa_potential)
    - [yukawa\_lambda](#yukawa_lambda)
    - [uramping](#uramping)
    - [omc](#omc)
    - [onsite\_radius](#onsite_radius)
  - [Spin-Constrained DFT](#spin-constrained-dft)
    - [sc\_mag\_switch](#sc_mag_switch)
    - [decay\_grad\_switch](#decay_grad_switch)
    - [sc\_thr](#sc_thr)
    - [nsc](#nsc)
    - [nsc\_min](#nsc_min)
    - [sc\_scf\_nmin](#sc_scf_nmin)
    - [alpha\_trial](#alpha_trial)
    - [sccut](#sccut)
    - [sc\_drop\_thr](#sc_drop_thr)
    - [sc\_scf\_thr](#sc_scf_thr)
  - [vdW correction](#vdw-correction)
    - [vdw\_method](#vdw_method)
    - [vdw\_s6](#vdw_s6)
    - [vdw\_s8](#vdw_s8)
    - [vdw\_a1](#vdw_a1)
    - [vdw\_a2](#vdw_a2)
    - [vdw\_d](#vdw_d)
    - [vdw\_abc](#vdw_abc)
    - [vdw\_c6\_file](#vdw_c6_file)
    - [vdw\_c6\_unit](#vdw_c6_unit)
    - [vdw\_r0\_file](#vdw_r0_file)
    - [vdw\_r0\_unit](#vdw_r0_unit)
    - [vdw\_cutoff\_type](#vdw_cutoff_type)
    - [vdw\_cutoff\_radius](#vdw_cutoff_radius)
    - [vdw\_radius\_unit](#vdw_radius_unit)
    - [vdw\_cutoff\_period](#vdw_cutoff_period)
    - [vdw\_cn\_thr](#vdw_cn_thr)
    - [vdw\_cn\_thr\_unit](#vdw_cn_thr_unit)
  - [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)
    - [berry\_phase](#berry_phase)
    - [gdir](#gdir)
    - [towannier90](#towannier90)
    - [nnkpfile](#nnkpfile)
    - [wannier\_method](#wannier_method)
    - [wannier\_spin](#wannier_spin)
    - [out\_wannier\_mmn](#out_wannier_mmn)
    - [out\_wannier\_amn](#out_wannier_amn)
    - [out\_wannier\_eig](#out_wannier_eig)
    - [out\_wannier\_unk](#out_wannier_unk)
    - [out\_wannier\_wvfn\_formatted](#out_wannier_wvfn_formatted)
  - [RT-TDDFT: Real-Time Time-Dependent Density Functional Theory](#rt-tddft-real-time-time-dependent-density-functional-theory)
    - [estep\_per\_md](#estep_per_md)
    - [td\_dt](#td_dt)
    - [td\_edm](#td_edm)
    - [td\_print\_eij](#td_print_eij)
    - [td\_propagator](#td_propagator)
    - [td\_vext](#td_vext)
    - [td\_vext\_dire](#td_vext_dire)
    - [td\_stype](#td_stype)
    - [td\_ttype](#td_ttype)
    - [td\_tstart](#td_tstart)
    - [td\_tend](#td_tend)
    - [td\_lcut1](#td_lcut1)
    - [td\_lcut2](#td_lcut2)
    - [td\_gauss\_freq](#td_gauss_freq)
    - [td\_gauss\_phase](#td_gauss_phase)
    - [td\_gauss\_sigma](#td_gauss_sigma)
    - [td\_gauss\_t0](#td_gauss_t0)
    - [td\_gauss\_amp](#td_gauss_amp)
    - [td\_trape\_freq](#td_trape_freq)
    - [td\_trape\_phase](#td_trape_phase)
    - [td\_trape\_t1](#td_trape_t1)
    - [td\_trape\_t2](#td_trape_t2)
    - [td\_trape\_t3](#td_trape_t3)
    - [td\_trape\_amp](#td_trape_amp)
    - [td\_trigo\_freq1](#td_trigo_freq1)
    - [td\_trigo\_freq2](#td_trigo_freq2)
    - [td\_trigo\_phase1](#td_trigo_phase1)
    - [td\_trigo\_phase2](#td_trigo_phase2)
    - [td\_trigo\_amp](#td_trigo_amp)
    - [td\_heavi\_t0](#td_heavi_t0)
    - [td\_heavi\_amp](#td_heavi_amp)
    - [init\_vecpot\_file](#init_vecpot_file)
    - [ocp](#ocp)
    - [ocp\_set](#ocp_set)
    - [out\_dipole](#out_dipole)
    - [out\_current](#out_current)
    - [out\_current\_k](#out_current_k)
    - [out\_efield](#out_efield)
    - [out\_vecpot](#out_vecpot)
  - [Variables useful for debugging](#variables-useful-for-debugging)
    - [nurse](#nurse)
    - [t\_in\_h](#t_in_h)
    - [vl\_in\_h](#vl_in_h)
    - [vnl\_in\_h](#vnl_in_h)
    - [vh\_in\_h](#vh_in_h)
    - [vion\_in\_h](#vion_in_h)
    - [test\_force](#test_force)
    - [test\_stress](#test_stress)
    - [test\_skip\_ewald](#test_skip_ewald)
  - [Electronic conductivities](#electronic-conductivities)
    - [cal\_cond](#cal_cond)
    - [cond\_che\_thr](#cond_che_thr)
    - [cond\_dw](#cond_dw)
    - [cond\_wcut](#cond_wcut)
    - [cond\_dt](#cond_dt)
    - [cond\_dtbatch](#cond_dtbatch)
    - [cond\_smear](#cond_smear)
    - [cond\_fwhm](#cond_fwhm)
    - [cond\_nonlocal](#cond_nonlocal)
  - [Implicit solvation model](#implicit-solvation-model)
    - [imp\_sol](#imp_sol)
    - [eb\_k](#eb_k)
    - [tau](#tau)
    - [sigma\_k](#sigma_k)
    - [nc\_k](#nc_k)
  - [Quasiatomic Orbital (QO) analysis](#quasiatomic-orbital-qo-analysis)
    - [qo\_switch](#qo_switch)
    - [qo\_basis](#qo_basis)
    - [qo\_strategy](#qo_strategy)
    - [qo\_screening\_coeff](#qo_screening_coeff)
    - [qo\_thr](#qo_thr)
  - [PEXSI](#pexsi)
    - [pexsi\_npole](#pexsi_npole)
    - [pexsi\_inertia](#pexsi_inertia)
    - [pexsi\_nmax](#pexsi_nmax)
    - [pexsi\_comm](#pexsi_comm)
    - [pexsi\_storage](#pexsi_storage)
    - [pexsi\_ordering](#pexsi_ordering)
    - [pexsi\_row\_ordering](#pexsi_row_ordering)
    - [pexsi\_nproc](#pexsi_nproc)
    - [pexsi\_symm](#pexsi_symm)
    - [pexsi\_trans](#pexsi_trans)
    - [pexsi\_method](#pexsi_method)
    - [pexsi\_nproc\_pole](#pexsi_nproc_pole)
    - [pexsi\_temp](#pexsi_temp)
    - [pexsi\_gap](#pexsi_gap)
    - [pexsi\_delta\_e](#pexsi_delta_e)
    - [pexsi\_mu\_lower](#pexsi_mu_lower)
    - [pexsi\_mu\_upper](#pexsi_mu_upper)
    - [pexsi\_mu](#pexsi_mu)
    - [pexsi\_mu\_thr](#pexsi_mu_thr)
    - [pexsi\_mu\_expand](#pexsi_mu_expand)
    - [pexsi\_mu\_guard](#pexsi_mu_guard)
    - [pexsi\_elec\_thr](#pexsi_elec_thr)
    - [pexsi\_zero\_thr](#pexsi_zero_thr)
  - [Linear Response TDDFT](#linear-response-tddft)
    - [ri\_hartree\_benchmark](#ri_hartree_benchmark)
    - [aims\_nbasis](#aims_nbasis)
  - [Linear Response TDDFT (Under Development Feature)](#linear-response-tddft-under-development-feature)
    - [xc\_kernel](#xc_kernel)
    - [lr\_init\_xc\_kernel](#lr_init_xc_kernel)
    - [lr\_solver](#lr_solver)
    - [lr\_thr](#lr_thr)
    - [nocc](#nocc)
    - [nvirt](#nvirt)
    - [lr\_nstates](#lr_nstates)
    - [lr\_unrestricted](#lr_unrestricted)
    - [abs\_wavelen\_range](#abs_wavelen_range)
    - [out\_wfc\_lr](#out_wfc_lr)
    - [abs\_gauge](#abs_gauge)
    - [abs\_broadening](#abs_broadening)
  - [Reduced Density Matrix Functional Theory](#reduced-density-matrix-functional-theory)
    - [rdmft](#rdmft)
    - [rdmft\_power\_alpha](#rdmft_power_alpha)

## System variables

### suffix

- **Type**: String
- **Description**: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the suffix is the name you can pick up for your convenience.
- **Default**: ABACUS

### ntype

- **Type**: Integer
- **Description**: Number of different atom species in the calculation.
- **Default**: 0

### calculation

- **Type**: String
- **Description**: Specify the type of calculation.

  - scf: perform self-consistent electronic structure calculations
  - nscf: perform non-self-consistent electronic structure calculations. A charge density file is required
  - relax: perform structure relaxation calculations, the relax_nmax parameter depicts the maximal number of ionic iterations
  - cell-relax: perform cell relaxation calculations
  - md: perform molecular dynamics simulations
  - get_pchg: obtain partial (band-decomposed) charge densities (for LCAO basis only). See out_pchg for more information
  - get_wf: obtain real space wave functions (for LCAO basis only). See out_wfc_norm and out_wfc_re_im for more information
  - get_s: obtain the overlap matrix formed by localized orbitals (for LCAO basis with multiple k points). the file name is SR.csr with file format being the same as that generated by out_mat_hs2
  - gen_bessel: generates projectors, i.e., a series of Bessel functions, for the DeePKS method (for LCAO basis only)
  - gen_opt_abfs: generate opt-ABFs as discussed in this article
  - test_memory: obtain a rough estimation of memory consumption for the calculation
  - test_neighbour: obtain information of neighboring atoms (for LCAO basis only), please specify a positive search_radius manually
- **Default**: scf

### esolver_type

- **Type**: String
- **Description**: Choose the energy solver.
  - ksdft: Kohn-Sham density functional theory
  - ofdft: orbital-free density functional theory
  - tdofdft: time-dependent orbital-free density functional theory
  - sdft: stochastic density functional theory
  - tddft: real-time time-dependent density functional theory (RT-TDDFT)
  - lj: Leonard Jones potential
  - dp: DeeP potential
  - nep: Neuroevolution Potential
  - ks-lr: Kohn-Sham density functional theory + LR-TDDFT (Under Development Feature)
  - lr: LR-TDDFT with given KS orbitals (Under Development Feature)
- **Default**: ksdft

### symmetry

- **Type**: String
- **Description**: Takes value 1, 0 or -1.
  - -1: No symmetry will be considered. It is recommended to set -1 for non-colinear + soc calculations, where time reversal symmetry is broken sometimes.
  - 0: Only time reversal symmetry would be considered in symmetry operations, which implied k point and -k point would be treated as a single k point with twice the weight.
  - 1: Symmetry analysis will be performed to determine the type of Bravais lattice and associated symmetry operations. (point groups, space groups, primitive cells, and irreducible k-points)

  > Note: When symmetry is enabled (value 1), k-points are reduced to the irreducible Brillouin zone (IBZ). For explicit k-point lists with custom weights (see KPT file), the custom weights are preserved during symmetry reduction. For Monkhorst-Pack grids, uniform weights are used.
- **Default**: default

### symmetry_prec

- **Type**: Real
- **Description**: The accuracy for symmetry analysis. Typically, the default value is good enough, but if the lattice parameters or atom positions in STRU file are not accurate enough, this value should be enlarged.
  > Note: if calculation==cell_relax, this value can be dynamically changed corresponding to the variation of accuracy of the lattice parameters and atom positions during the relaxation.
- **Default**: 1.0e-6
- **Unit**: Bohr

### symmetry_autoclose

- **Type**: Boolean
- **Availability**: *symmetry==1*
- **Description**: Control how to deal with error in symmetry analysis due to inaccurate lattice parameters or atom positions in STRU file, especially useful when calculation==cell-relax
  - False: quit with an error message
  - True: automatically set symmetry to 0 and continue running without symmetry analysis
- **Default**: True

### cal_force

- **Type**: Boolean
- **Description**: If set to True, calculate the force at the end of the electronic iteration.
- **Default**: False

### kpar

- **Type**: Integer
- **Description**: Divide all processors into kpar groups, and k points will be distributed among each group. The value taken should be less than or equal to the number of k points as well as the number of MPI processes.
- **Default**: 1

### bndpar

- **Type**: Integer
- **Description**: Divide all processors into bndpar groups, and bands (only stochastic orbitals now) will be distributed among each group. It should be larger than 0.
- **Default**: 1

### latname

- **Type**: String
- **Description**: Specifies the type of Bravias lattice. When set to none, the three lattice vectors are supplied explicitly in STRU file.

  Available options are:

  - none: free structure
  - sc: simple cubic
  - fcc: face-centered cubic
  - bcc: body-centered cubic
  - hexagonal: hexagonal
  - trigonal: trigonal
  - st: simple tetragonal
  - bct: body-centered tetragonal
  - so: orthorhombic
  - baco: base-centered orthorhombic
  - fco: face-centered orthorhombic
  - bco: body-centered orthorhombic
  - sm: simple monoclinic
  - bacm: base-centered monoclinic
  - triclinic: triclinic
- **Default**: none

### init_wfc

- **Type**: String
- **Description**: The type of the starting wave functions.

  Available options are:

  - atomic: from atomic pseudo wave functions. If they are not enough, other wave functions are initialized with random numbers.
  - atomic+random: add small random numbers on atomic pseudo-wavefunctions
  - file: from binary files wf*.dat, which are output by setting out_wfc_pw to 2.
  - random: random numbers
  - nao: from numerical atomic orbitals. If they are not enough, other wave functions are initialized with random numbers.
  - nao+random: add small random numbers on numerical atomic orbitals

  > Note: Only the file option is useful for the lcao basis set, which is mostly used when calculation is set to get_wf and get_pchg.
- **Default**: atomic

### init_chg

- **Type**: String
- **Description**: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density.

  - atomic: the density is starting from the summation of the atomic density of single atoms.
  - file: the density will be read in from a binary file charge-density.dat first. If it does not exist, the charge density will be read in from cube files.
  - wfc: the density will be calculated by wavefunctions and occupations.
  - dm: the density will be calculated by real space density matrix(DMR) of LCAO base.
  - hr: the real space Hamiltonian matrix(HR) will be read in from file hrs1_nao.csr in directory read_file_dir.
  - auto: Abacus first attempts to read the density from a file; if not found, it defaults to using atomic density.
- **Default**: atomic

### init_vel

- **Type**: Boolean
- **Description**: - True: read the atom velocity (atomic unit : 1 a.u. = 21.877 Angstrom/fs) from the atom file (STRU) and determine the initial temperature md_tfirst. If md_tfirst is unset or less than zero, init_vel is autoset to be true.
  - False: assign value to atom velocity using Gaussian distributed random numbers.
- **Default**: False

### mem_saver

- **Type**: Integer
- **Availability**: *Used only for nscf calculations with plane wave basis set.*
- **Description**: Save memory when performing nscf calculations.
  - 0: no memory saving techniques are used.
  - 1: a memory saving technique will be used for many k point calculations.
- **Default**: 0

### cal_stress

- **Type**: Boolean
- **Description**: If set to True, calculate the stress at the end of the electronic iteration.
- **Default**: False

### diago_proc

- **Type**: Integer
- **Availability**: *Used only for plane wave basis set.*
- **Description**: - 0: it will be set to the number of MPI processes.
  - &gt;0: it specifies the number of processes used for carrying out diagonalization. Must be less than or equal to total number of MPI processes.
- **Default**: 0

### nbspline

- **Type**: Integer
- **Description**: If set to a natural number, a Cardinal B-spline interpolation will be used to calculate Structure Factor. nbspline represents the order of B-spline basis and a larger one can get more accurate results but cost more. It is turned off by default.
- **Default**: -1

### kspacing

- **Type**: Vector of Real (1 or 3 values)
- **Description**: Set the smallest allowed spacing between k points, unit in 1/bohr. It should be larger than 0.0, and suggest smaller than 0.25. When you have set this value &gt; 0.0, then the KPT file is unnecessary. The default value 0.0 means that ABACUS will read the applied KPT file.

  > Note: If gamma_only is set to be true, kspacing is invalid.
- **Default**: 0.0

### min_dist_coef

- **Type**: Real
- **Description**: A factor related to the allowed minimum distance between two atoms. At the beginning, ABACUS will check the structure, and if the distance of two atoms is shorter than min_dist_coef*(standard covalent bond length), we think this structure is unreasonable.
- **Default**: 0.2

### device

- **Type**: String
- **Description**: Specifies the computing device for ABACUS.

  Available options are:

  - cpu: for CPUs via Intel, AMD, or Other supported CPU devices
  - gpu: for GPUs via CUDA or ROCm.

  > Note: ks_solver must also be set to the algorithms supported. lcao_in_pw currently does not support gpu.
- **Default**: cpu

### precision

- **Type**: String
- **Availability**: *Used only for plane wave basis set.*
- **Description**: Specifies the precision when performing scf calculation.
  - single: single precision
  - double: double precision
- **Default**: double

### timer_enable_nvtx

- **Type**: Boolean
- **Description**: Controls whether NVTX profiling labels are emitted by the timer. This feature is only effective on CUDA platforms.

  - True: Enable NVTX profiling labels in the timer.
  - False: Disable NVTX profiling labels in the timer.
- **Default**: False

### cell_factor

- **Type**: Real
- **Description**: Used in the construction of the pseudopotential tables. For cell-relax calculations, this is automatically set to 2.0.
- **Default**: 1.2

### dm_to_rho

- **Type**: Boolean
- **Description**: Reads density matrix in npz format and calculates electron density.
- **Default**: False

### chg_extrap

- **Type**: String
- **Description**: Charge extrapolation method for MD and relaxation calculations.
- **Default**: default

### nb2d

- **Type**: Integer
- **Description**: In LCAO calculations, the Hamiltonian and overlap matrices are distributed across 2D processor grid. This parameter controls the 2D block size for distribution.
- **Default**: 0

### cal_symm_repr

- **Type**: Integer \[Integer\](optional)
- **Description**: Whether to print the matrix representation of symmetry operation to running log file. If the first value is given as 1, then all matrix representations will be printed. The second optional parameter controls the precision (number of digits) to print, default is 3, which is enough for a quick check.
- **Default**: 1 3

[back to top](#full-list-of-input-keywords)

## Input files

### stru_file

- **Type**: String
- **Description**: The name of the structure file containing various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
- **Default**: STRU

### kpoint_file

- **Type**: String
- **Description**: The name of the k-point file that includes the k-point information of Brillouin zone.
- **Default**: KPT

### pseudo_dir

- **Type**: String
- **Description**: The directory of pseudopotential files. This parameter is combined with the pseudopotential filenames in the STRU file to form the complete pseudopotential file paths.
- **Default**: ""

### orbital_dir

- **Type**: String
- **Description**: The directory to save numerical atomic orbitals. This parameter is combined with orbital filenames in the STRU file to form the complete orbital file paths.
- **Default**: ""

### read_file_dir

- **Type**: String
- **Description**: Location of files, such as the electron density (chgs1.cube), required as a starting point.
- **Default**: OUT.$suffix

### restart_load

- **Type**: Boolean
- **Availability**: *Used only when numerical atomic orbitals are employed as basis set.*
- **Description**: If restart_save is set to true and an electronic iteration is finished, calculations can be restarted from the charge density file, which are saved in the former calculation.
- **Default**: False

### spillage_outdir

- **Type**: String
- **Availability**: *Used only for plane wave basis set.*
- **Description**: The directory to save the spillage files.
- **Default**: "./"

[back to top](#full-list-of-input-keywords)

## Plane wave related variables

### ecutwfc

- **Type**: Real
- **Description**: Energy cutoff for plane wave functions. Note that even for localized orbitals basis, you still need to setup an energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set.
  > Note: ecutwfc and ecutrho can be set simultaneously. If only one parameter is set, abacus will automatically set another parameter based on the 4-time relationship.
- **Default**: 50 for PW basis, 100 for LCAO basis
- **Unit**: Ry

### ecutrho

- **Type**: Real
- **Description**: Energy cutoff for charge density and potential. For norm-conserving pseudopotential you should stick to the default value, you can reduce it by a little but it will introduce noise especially on forces and stress.
- **Default**: 4*ecutwfc
- **Unit**: Ry

### nx

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points in x direction. If set to 0, the number will be calculated from ecutrho.

  > Note: You must specify all three dimensions (nx, ny, nz) for this setting to be used.
- **Default**: 0

### ny

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points in y direction. If set to 0, the number will be calculated from ecutrho.

  > Note: You must specify all three dimensions (nx, ny, nz) for this setting to be used.
- **Default**: 0

### nz

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points in z direction. If set to 0, the number will be calculated from ecutrho.

  > Note: You must specify all three dimensions (nx, ny, nz) for this setting to be used.
- **Default**: 0

### ndx

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points for the dense part of charge density in x direction. If set to 0, the number will be calculated from ecutwfc.

  > Note: You must specify all three dimensions (ndx, ndy, ndz) for this setting to be used. These parameters must be used combined with nx, ny, nz. If nx, ny, nz are unset, ndx, ndy, ndz are used as nx, ny, nz.
- **Default**: 0

### ndy

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points for the dense part of charge density in y direction. If set to 0, the number will be calculated from ecutwfc.

  > Note: You must specify all three dimensions (ndx, ndy, ndz) for this setting to be used. These parameters must be used combined with nx, ny, nz. If nx, ny, nz are unset, ndx, ndy, ndz are used as nx, ny, nz.
- **Default**: 0

### ndz

- **Type**: Integer
- **Description**: If set to a positive number, specifies the number of FFT grid points for the dense part of charge density in z direction. If set to 0, the number will be calculated from ecutwfc.

  > Note: You must specify all three dimensions (ndx, ndy, ndz) for this setting to be used. These parameters must be used combined with nx, ny, nz. If nx, ny, nz are unset, ndx, ndy, ndz are used as nx, ny, nz.
- **Default**: 0

### pw_seed

- **Type**: Integer
- **Availability**: *Only used for plane wave basis.*
- **Description**: Specify the random seed to initialize wave functions. Only positive integers are available.
- **Default**: 0

### diag_subspace

- **Type**: Integer
- **Description**: The method to diagonalize subspace in dav_subspace method.
  - 0: by LAPACK
  - 1: by GenELPA
  - 2: by ScaLAPACK
- **Default**: 0

### erf_ecut

- **Type**: Real
- **Description**: Used in variable-cell molecular dynamics (or in stress calculation). See erf_sigma for details.
- **Default**: 0.0
- **Unit**: Ry

### fft_mode

- **Type**: Integer
- **Description**: Set the mode of FFTW.
  - 0: FFTW_ESTIMATE
  - 1: FFTW_MEASURE
  - 2: FFTW_PATIENT
  - 3: FFTW_EXHAUSTIVE
- **Default**: 0

### erf_height

- **Type**: Real
- **Description**: Used in variable-cell molecular dynamics (or in stress calculation). See erf_sigma for details.
- **Default**: 0.0
- **Unit**: Ry

### erf_sigma

- **Type**: Real
- **Description**: In order to recover the accuracy of a constant energy cutoff calculation, the kinetic functional is modified, which is used in variable-cell molecular dynamics (or in stress calculation).
- **Default**: 0.1
- **Unit**: Ry

### pw_diag_thr

- **Type**: Real
- **Description**: Only used when you use ks_solver = cg/dav/dav_subspace/bpcg. It indicates the threshold for the first electronic iteration, from the second iteration the pw_diag_thr will be updated automatically. For nscf calculations with planewave basis set, pw_diag_thr should be &lt;= 1e-3.
- **Default**: 0.01

### diago_smooth_ethr

- **Type**: Boolean
- **Description**: If TRUE, the smooth threshold strategy, which applies a larger threshold (10e-5) for the empty states, will be implemented in the diagonalization methods. (This strategy should not affect total energy, forces, and other ground-state properties, but computational efficiency will be improved.) If FALSE, the smooth threshold strategy will not be applied.
- **Default**: false

### use_k_continuity

- **Type**: Boolean
- **Availability**: *Used only for plane wave basis set.*
- **Description**: If TRUE, the wavefunctions at k-point will be initialized from the converged wavefunctions at the nearest k-point, which can speed up the SCF convergence. Only works for PW basis.
- **Default**: false

### pw_diag_nmax

- **Type**: Integer
- **Description**: Only useful when you use ks_solver = cg/dav/dav_subspace/bpcg. It indicates the maximal iteration number for cg/david/dav_subspace/bpcg method.
- **Default**: 40

### pw_diag_ndim

- **Type**: Integer
- **Description**: Only useful when you use ks_solver = dav or ks_solver = dav_subspace. It indicates dimension of workspace(number of wavefunction packets, at least 2 needed) for the Davidson method. A larger value may yield a smaller number of iterations in the algorithm but uses more memory and more CPU time in subspace diagonalization.
- **Default**: 4

### diago_cg_prec

- **Type**: Integer
- **Description**: Preconditioner type for conjugate gradient diagonalization method.
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## Numerical atomic orbitals related variables

### lmaxmax

- **Type**: Integer
- **Description**: If not equals to 2, then the maximum l channels on LCAO is set to lmaxmax. If 2, then the number of l channels will be read from the LCAO data sets. Normally no input should be supplied for this variable so that it is kept as its default.
- **Default**: 2.

### lcao_ecut

- **Type**: Real
- **Description**: Energy cutoff (in Ry) for two-center integrals in LCAO. The two-center integration table are obtained via a k space integral whose upper limit is about sqrt(lcao_ecut).
- **Default**: ecutwfc

### lcao_dk

- **Type**: Real
- **Description**: the interval of k points for two-center integrals. The two-center integration table are obtained via a k space integral on a uniform grid with spacing lcao_dk.
- **Default**: 0.01
- **Unit**: Bohr

### lcao_dr

- **Type**: Real
- **Description**: r spacing of the integration table of two-center integrals.
- **Default**: 0.01
- **Unit**: Bohr

### lcao_rmax

- **Type**: Real
- **Description**: Maximum distance for the two-center integration table.
- **Default**: 30
- **Unit**: Bohr

### search_radius

- **Type**: Real
- **Description**: Searching radius in finding the neighbouring atoms. By default the radius will be automatically determined by the cutoffs of orbitals and nonlocal beta projectors.
- **Default**: -1
- **Unit**: Bohr

### bx

- **Type**: Integer
- **Description**: In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.
- **Default**: 0

### by

- **Type**: Integer
- **Description**: In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.
- **Default**: 0

### bz

- **Type**: Integer
- **Description**: In the matrix operation of grid integral, bx/by/bz grids (in x, y, z directions) are treated as a whole as a matrix element. A different value will affect the calculation speed. The default is 0, which means abacus will automatically calculate these values.
- **Default**: 0

### elpa_num_thread

- **Type**: Integer
- **Description**: Number of threads used in one elpa calculation.

  If the number is below 0 or 0 or beyond the max number of threads, all elpa calculation will be using all mpi threads
- **Default**: -1

### num_stream

- **Type**: Integer
- **Description**: The number of CUDA streams used in LCAO calculations with GPU acceleration.
- **Default**: 4

[back to top](#full-list-of-input-keywords)

## Electronic structure

### basis_type

- **Type**: String
- **Description**: Choose the basis set.
  - pw: Using plane-wave basis set only.
  - lcao: Using localized atomic orbital sets.
  - lcao_in_pw: Expand the localized atomic set in plane-wave basis, non-self-consistent field calculation not tested.
- **Default**: pw

### ks_solver

- **Type**: String
- **Description**: Choose the diagonalization methods for the Hamiltonian matrix expanded in a certain basis set.

  For plane-wave basis,

  - cg: The conjugate-gradient (CG) method.
  - bpcg: The BPCG method, which is a block-parallel Conjugate Gradient (CG) method, typically exhibits higher acceleration in a GPU environment.
  - dav: The Davidson algorithm.
  - dav_subspace: The Davidson algorithm without orthogonalization operation, this method is the most recommended for efficiency. pw_diag_ndim can be set to 2 for this method.

  For numerical atomic orbitals basis,

  - lapack: Use LAPACK to diagonalize the Hamiltonian, only used for serial version
  - genelpa: Use GEN-ELPA to diagonalize the Hamiltonian.
  - scalapack_gvx: Use Scalapack to diagonalize the Hamiltonian.
  - cusolver: Use CUSOLVER to diagonalize the Hamiltonian, at least one GPU is needed.
  - cusolvermp: Use CUSOLVER to diagonalize the Hamiltonian, supporting multi-GPU devices. Note that you should set the number of MPI processes equal to the number of GPUs.
  - elpa: The ELPA solver supports both CPU and GPU. By setting the device to GPU, you can launch the ELPA solver with GPU acceleration (provided that you have installed a GPU-supported version of ELPA, which requires you to manually compile and install ELPA, and the ABACUS should be compiled with -DUSE_ELPA=ON and -DUSE_CUDA=ON). The ELPA solver also supports multi-GPU acceleration.

  If you set ks_solver=genelpa for basis_type=pw, the program will stop with an error message:

  ``text genelpa can not be used with plane wave basis. ``

  Then the user has to correct the input file and restart the calculation.

### nbands

- **Type**: Integer
- **Description**: The number of Kohn-Sham orbitals to calculate. It is recommended to setup this value, especially when smearing techniques are utilized, more bands should be included.

### nelec

- **Type**: Real
- **Description**: - 0.0: The total number of electrons will be calculated by the sum of valence electrons (i.e. assuming neutral system).
  - &gt;0.0: this denotes the total number of electrons in the system. Must be less than 2*nbands.
- **Default**: 0.0

### nelec_delta

- **Type**: Real
- **Description**: The total number of electrons will be calculated by nelec+nelec_delta.
- **Default**: 0.0

### nupdown

- **Type**: Real
- **Description**: - 0.0: no constrain apply to system.
  - &gt;0.0: The different number of electrons between spin-up and spin-down channels. The range of value must be in [-nelec ~ nelec]. It is one type of constrainted DFT method, two Fermi energies will be calculated.
- **Default**: 0.0

### dft_functional

- **Type**: String
- **Description**: In our package, the XC functional can either be set explicitly using the dft_functional keyword in INPUT file. If dft_functional is not specified, ABACUS will use the xc functional indicated in the pseudopotential file. On the other hand, if dft_functional is specified, it will overwrite the functional from pseudopotentials and performs calculation with whichever functional the user prefers. We further offer two ways of supplying exchange-correlation functional. The first is using 'short-hand' names. A complete list of 'short-hand' expressions can be found in the source code. Supported density functionals are:
  - LDA functionals
  - LDA (equivalent with PZ and SLAPZNOGXNOGC), PWLDA
  - GGA functionals
  - PBE (equivalent with SLAPWPBXPBC), PBESOL, REVPBE, WC, BLYP, BP(referred to BP86), PW91, HCTH, OLYP, BLYP_LR
  - meta-GGA functionals
  - SCAN (require LIBXC)
  - Hybrid functionals
  - PBE0, HF
  - If LIBXC is available, additional short-hand names of hybrid functionals are supported: HSE(referred to HSE06), B3LYP, LC_PBE, LC_WPBE, LRC_WPBE, LRC_WPBEH, CAM_PBEH, WP22, CWP22, MULLER (equivalent with POWER)
  - Hybrid meta-GGA functionals
  - SCAN0 (require LIBXC)

  The other way is only available when compiling with LIBXC, and it allows for supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by a plus sign, for example, dft_functional='LDA_X_1D_EXPONENTIAL+LDA_C_1D_CSC'.
- **Default**: Used the same as DFT functional as specified in the pseudopotential files.

### xc_temperature

- **Type**: Real
- **Description**: Specifies temperature when using temperature-dependent XC functionals (KSDT and so on).
- **Default**: 0.0
- **Unit**: Ry

### xc_exch_ext

- **Type**: Integer followed by Real values
- **Description**: Customized parameterization on the exchange part of XC functional. The first value should be the LibXC ID of the original functional, and latter values are external parameters. Default values are those of Perdew-Burke-Ernzerhof (PBE) functional. For more information on LibXC ID of functionals, please refer to LibXC. For parameters of functionals of interest, please refer to the source code of LibXC, such as PBE functional interface in LibXC: gga_x_pbe.c.

  > Note: Solely setting this keyword will take no effect on XC functionals. One should also set dft_functional to the corresponding functional to apply the customized parameterization. Presently this feature can only support parameterization on one exchange functional.
- **Default**: 101 0.8040 0.2195149727645171

### xc_corr_ext

- **Type**: Integer followed by Real values
- **Description**: Customized parameterization on the correlation part of XC functional. The first value should be the LibXC ID of the original functional, and latter values are external parameters. Default values are those of Perdew-Burke-Ernzerhof (PBE) functional. For more information on LibXC ID of functionals, please refer to LibXC. For parameters of functionals of interest, please refer to the source code of LibXC, such as PBE functional interface in LibXC: gga_c_pbe.c.

  > Note: Solely setting this keyword will take no effect on XC functionals. One should also set dft_functional to the corresponding functional to apply the customized parameterization. Presently this feature can only support parameterization on one correlation functional.
- **Default**: 130 0.06672455060314922 0.031090690869654895034 1.0

### pseudo_rcut

- **Type**: Real
- **Description**: Cut-off of radial integration for pseudopotentials.
- **Default**: 15
- **Unit**: Bohr

### pseudo_mesh

- **Type**: Boolean
- **Description**: - 0: Use a mesh for radial integration of pseudopotentials.
  - 1: Use the mesh that is consistent with quantum espresso
- **Default**: 0

### nspin

- **Type**: Integer
- **Description**: The number of spin components of wave functions.
  - 1: Spin degeneracy
  - 2: Collinear spin polarized.
  - 4: For the case of noncollinear polarized, nspin will be automatically set to 4 without being specified by the user.
- **Default**: 1

### smearing_method

- **Type**: String
- **Description**: It indicates which occupation and smearing method is used in the calculation.
  - fixed: fixed occupations (available for non-coductors only)
  - gauss or gaussian: Gaussian smearing method.
  - mp: methfessel-paxton smearing method; recommended for metals.
  - mp2: 2-nd methfessel-paxton smearing method; recommended for metals.
  - mv or cold: marzari-vanderbilt smearing method.
  - fd: Fermi-Dirac smearing method: and smearing_sigma below is the temperature (in Ry).
- **Default**: gauss

### smearing_sigma

- **Type**: Real
- **Description**: Energy range for smearing.
- **Default**: 0.015
- **Unit**: Ry

### smearing_sigma_temp

- **Type**: Real
- **Description**: Energy range for smearing, smearing_sigma = 1/2 kB smearing_sigma_temp.
- **Default**: 2 * smearing_sigma / kB.
- **Unit**: K

### mixing_type

- **Type**: String
- **Description**: Charge mixing methods.
  - plain: Just simple mixing.
  - pulay: Standard Pulay method. P. Pulay Chemical Physics Letters, (1980)
  - broyden: Simplified modified Broyden method. D.D. Johnson Physical Review B (1988)

  In general, the convergence of the Broyden method is slightly faster than that of the Pulay method.
- **Default**: broyden

### mixing_beta

- **Type**: Real
- **Description**: In general, the formula of charge mixing can be written as rho_new = rho_old + mixing_beta * drho, where rho_new represents the new charge density after charge mixing, rho_old represents the charge density in previous step, drho is obtained through various mixing methods, and mixing_beta is set by this parameter. A lower value of 'mixing_beta' results in less influence of drho on rho_new, making the self-consistent field (SCF) calculation more stable. However, it may require more steps to achieve convergence. We recommend the following options:
  - 0.8: nspin=1
  - 0.4: nspin=2 and nspin=4
  - 0: keep charge density unchanged, usually used for restarting with init_chg=file or testing.
  - 0.1 or less: if convergence of SCF calculation is difficult to reach, please try 0 &lt; mixing_beta &lt; 0.1.

  Note: For low-dimensional large systems, the setup of mixing_beta=0.1, mixing_ndim=20, and mixing_gg0=1.0 usually works well.
- **Default**: 0.8 for nspin=1, 0.4 for nspin=2 and nspin=4.

### mixing_beta_mag

- **Type**: Real
- **Description**: Mixing parameter of magnetic density.
- **Default**: 4*mixing_beta, but the maximum value is 1.6.

### mixing_ndim

- **Type**: Integer
- **Description**: It indicates the mixing dimensions in Pulay or Broyden. Pulay and Broyden method use the density from previous mixing_ndim steps and do a charge mixing based on this density.

  For systems that are difficult to converge, one could try increasing the value of 'mixing_ndim' to enhance the stability of the self-consistent field (SCF) calculation.
- **Default**: 8

### mixing_restart

- **Type**: Real
- **Description**: If the density difference between input and output drho is smaller than mixing_restart, SCF will restart at next step which means SCF will restart by using output charge density from perivos iteration as input charge density directly, and start a new mixing. Notice that mixing_restart will only take effect once in one SCF.
- **Default**: 0

### mixing_dmr

- **Type**: Boolean
- **Availability**: *Only for mixing_restart &gt;= 0.0*
- **Description**: At n-th iteration which is calculated by drho&lt;mixing_restart, SCF will start a mixing for real-space density matrix by using the same coefficiences as the mixing of charge density.
- **Default**: false

### mixing_gg0

- **Type**: Real
- **Description**: Whether to perfom Kerker scaling for charge density.
  - &gt;0: The high frequency wave vectors will be suppressed by multiplying a scaling factor. Setting mixing_gg0 = 1.0 is normally a good starting point. Kerker preconditioner will be automatically turned off if mixing_beta &lt;= 0.1.
  - 0: No Kerker scaling is performed.

  For systems that are difficult to converge, particularly metallic systems, enabling Kerker scaling may aid in achieving convergence.
- **Default**: 1.0

### mixing_gg0_mag

- **Type**: Real
- **Description**: Whether to perfom Kerker preconditioner of magnetic density. Note: we do not recommand to open Kerker preconditioner of magnetic density unless the system is too hard to converge.
- **Default**: 0.0

### mixing_gg0_min

- **Type**: Real
- **Description**: The minimum kerker coefficient.
- **Default**: 0.1

### mixing_angle

- **Type**: Real
- **Availability**: *Only relevant for non-colinear calculations nspin=4.*
- **Description**: Normal broyden mixing can give the converged result for a given magnetic configuration. If one is not interested in the energies of a given magnetic configuration but wants to determine the ground state by relaxing the magnetic moments' directions, one cannot rely on the standard Broyden mixing algorithm. To enhance the ability to find correct magnetic configuration for non-colinear calculations, ABACUS implements a promising mixing method proposed by J. Phys. Soc. Jpn. 82 (2013) 114706. Here, mixing_angle is the angle mixing parameter. In fact, only mixing_angle=1.0 is implemented currently.
  - &lt;=0: Normal broyden mixing
  - &gt;0: Angle mixing for the modulus with mixing_angle=1.0
- **Default**: -10.0

### mixing_tau

- **Type**: Boolean
- **Availability**: *Only relevant for meta-GGA calculations.*
- **Description**: Whether to mix the kinetic energy density.
  - True: The kinetic energy density will also be mixed. It seems for general cases, SCF converges fine even without this mixing. However, if there is difficulty in converging SCF for meta-GGA, it might be helpful to turn this on.
  - False: The kinetic energy density will not be mixed.
- **Default**: False

### mixing_dftu

- **Type**: Boolean
- **Availability**: *Only relevant for DFT+U calculations.*
- **Description**: Whether to mix the occupation matrices.
  - True: The occupation matrices will also be mixed by plain mixing. From experience this is not very helpful if the +U calculation does not converge.
  - False: The occupation matrices will not be mixed.
- **Default**: False

### gamma_only

- **Type**: Boolean
- **Availability**: *Only used in localized orbitals set*
- **Description**: Whether to use gamma_only algorithm.
  - 0: more than one k-point is used and the ABACUS is slower compared to the gamma only algorithm.
  - 1: ABACUS uses gamma only, the algorithm is faster and you don't need to specify the k-points file.

  Note: If gamma_only is set to 1, the KPT file will be overwritten. So make sure to turn off gamma_only for multi-k calculations.
- **Default**: 0

### scf_nmax

- **Type**: Integer
- **Description**: This variable indicates the maximal iteration number for electronic iterations.
- **Default**: 100

### scf_thr

- **Type**: Real
- **Description**: It's the density threshold for electronic iteration. It represents the charge density error between two sequential densities from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
- **Default**: 1.0e-9 (plane-wave basis), or 1.0e-7 (localized atomic orbital basis).
- **Unit**: Ry if scf_thr_type=1, dimensionless if scf_thr_type=2

### scf_ene_thr

- **Type**: Real
- **Description**: It's the energy threshold for electronic iteration. It represents the total energy error between two sequential densities from electronic iterations.
- **Default**: -1.0. If the user does not set this parameter, it will not take effect.
- **Unit**: eV

### scf_thr_type

- **Type**: Integer
- **Description**: Choose the calculation method of convergence criterion.
  - 1: the criterion is defined in reciprocal space, which is used in SCF of PW basis with unit Ry.
  - 2: the criterion is defined in real space, where is the number of electron, which is used in SCF of LCAO with unit dimensionless.
- **Default**: 1 (plane-wave basis), or 2 (localized atomic orbital basis).

### scf_os_stop

- **Type**: Boolean
- **Description**: For systems that are difficult to converge, the SCF process may exhibit oscillations in charge density, preventing further progress toward the specified convergence criteria and resulting in continuous oscillation until the maximum number of steps is reached; this greatly wastes computational resources. To address this issue, this function allows ABACUS to terminate the SCF process early upon detecting oscillations, thus reducing subsequent meaningless calculations. The detection of oscillations is based on the slope of the logarithm of historical drho values. To this end, Least Squares Method is used to calculate the slope of the logarithmically taken drho for the previous scf_os_ndim iterations. If the calculated slope is larger than scf_os_thr, stop the SCF.

  - 0: The SCF will continue to run regardless of whether there is oscillation or not.
  - 1: If the calculated slope is larger than scf_os_thr, stop the SCF.
- **Default**: false

### scf_os_thr

- **Type**: Real
- **Description**: The slope threshold to determine if the SCF is stuck in a charge density oscillation. If the calculated slope is larger than scf_os_thr, stop the SCF.
- **Default**: -0.01

### scf_os_ndim

- **Type**: Integer
- **Description**: To determine the number of old iterations' drho used in slope calculations.
- **Default**: mixing_ndim

### sc_os_ndim

- **Type**: Integer
- **Description**: To determine the number of old iterations to judge oscillation, it occured, more accurate lambda with DeltaSpin method would be calculated, only for PW base.
- **Default**: 5

### lspinorb

- **Type**: Boolean
- **Description**: Whether to consider spin-orbit coupling (SOC) effect in the calculation.
  - True: Consider spin-orbit coupling effect. When enabled:
  - nspin is automatically set to 4 (noncollinear spin representation)
  - Symmetry is automatically disabled (SOC breaks inversion symmetry)
  - Requires full-relativistic pseudopotentials with has_so=true in the UPF header
  - False: Do not consider spin-orbit coupling effect.
  - Common Error: "no soc upf used for lspinorb calculation" - ensure you are using full-relativistic pseudopotentials
- **Default**: False

### noncolin

- **Type**: Boolean
- **Description**: Whether to allow non-collinear magnetic moments, where magnetization can point in arbitrary directions (x, y, z components) rather than being constrained to the z-axis.
  - True: Allow non-collinear polarization. When enabled:
  - nspin is automatically set to 4
  - Wave function dimension is doubled (npol=2), and the number of occupied states is doubled
  - Charge density has 4 components (Pauli spin matrices)
  - Cannot be used with gamma_only=true
  - Can be combined with lspinorb=true for SOC effects with non-collinear magnetism
  - False: Do not allow non-collinear polarization (magnetization constrained to z-axis).
  - Relationship with lspinorb:
  - noncolin=0, lspinorb=1: SOC with z-axis magnetism only (for non-magnetic materials with SOC)
  - noncolin=1, lspinorb=0: Non-collinear magnetism without SOC
  - noncolin=1, lspinorb=1: Both non-collinear magnetism and SOC
- **Default**: False

### soc_lambda

- **Type**: Real
- **Availability**: *Only works when lspinorb=true*
- **Description**: Modulates the strength of spin-orbit coupling effect. Sometimes, for some real materials, both scalar-relativistic and full-relativistic pseudopotentials cannot describe the exact spin-orbit coupling. Artificial modulation may help in such cases.

  soc_lambda, which has value range [0.0, 1.0], is used to modulate SOC effect:

  - soc_lambda 0.0: Scalar-relativistic case (no SOC)
  - soc_lambda 1.0: Full-relativistic case (full SOC)
  - Intermediate values: Partial-relativistic SOC (interpolation between scalar and full)

  Use case: When experimental or high-level theoretical results suggest that the SOC effect is weaker or stronger than what full-relativistic pseudopotentials predict, you can adjust this parameter to match the target behavior.
- **Default**: 1.0

### dfthalf_type

- **Type**: Integer
- **Description**: DFT-1/2 type:
  - 0: DFT-1/2 is off.
  - 1: Shell DFT-1/2 method is used.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electronic structure (SDFT)

### method_sto

- **Type**: Integer
- **Availability**: *esolver_type = sdft*
- **Description**: Different methods to do stochastic DFT
  - 1: Calculate twice, this method cost less memory but is slower.
  - 2: Calculate once but needs much more memory. This method is much faster. Besides, it calculates with a smaller nche_sto. However, when the memory is not enough, only method 1 can be used.
  - other: use 2
- **Default**: 2

### nbands_sto

- **Type**: Integer or string
- **Availability**: *esolver_type = sdft*
- **Description**: The number of stochastic orbitals
  - &gt; 0: Perform stochastic DFT. Increasing the number of bands improves accuracy and reduces stochastic errors; To perform mixed stochastic-deterministic DFT, you should set nbands, which represents the number of KS orbitals.
  - 0: Perform Kohn-Sham DFT.
  - all: All complete basis sets are used to replace stochastic orbitals with the Chebyshev method (CT), resulting in the same results as KSDFT without stochastic errors.
- **Default**: 256

### nche_sto

- **Type**: Integer
- **Availability**: *esolver_type = sdft*
- **Description**: Chebyshev expansion orders for stochastic DFT.
- **Default**: 100

### emin_sto

- **Type**: Real
- **Availability**: *esolver_type = sdft*
- **Description**: Trial energy to guess the lower bound of eigen energies of the Hamiltonian Operator.
- **Default**: 0.0
- **Unit**: Ry

### emax_sto

- **Type**: Real
- **Availability**: *esolver_type = sdft*
- **Description**: Trial energy to guess the upper bound of eigen energies of the Hamiltonian Operator.
- **Default**: 0.0
- **Unit**: Ry

### seed_sto

- **Type**: Integer
- **Availability**: *esolver_type = sdft*
- **Description**: The random seed to generate stochastic orbitals.
  - &gt;= 0: Stochastic orbitals have the form of exp(i*theta), where theta is a uniform distribution in [0, 2*pi).
  - 0: the seed is decided by time(NULL).
  - &lt;= -1: Stochastic orbitals have the form of +1 or -1 with equal probability.
  - -1: the seed is decided by time(NULL).
- **Default**: 0

### initsto_ecut

- **Type**: Real
- **Availability**: *esolver_type = sdft*
- **Description**: Stochastic wave functions are initialized in a large box generated by "4*initsto_ecut". initsto_ecut should be larger than ecutwfc. In this method, SDFT results are the same when using different cores. Besides, coefficients of the same G are the same when ecutwfc is rising to initsto_ecut. If it is smaller than ecutwfc, it will be turned off.
- **Default**: 0.0
- **Unit**: Ry

### initsto_freq

- **Type**: Integer
- **Availability**: *esolver_type = sdft*
- **Description**: Frequency (once each initsto_freq steps) to generate new stochastic orbitals when running md.
  - positive integer: Update stochastic orbitals
  - 0: Never change stochastic orbitals.
- **Default**: 0

### npart_sto

- **Type**: Integer
- **Availability**: *method_sto = 2 and out_dos = 1 or cal_cond = True*
- **Description**: Make memory cost to 1/npart_sto times of the previous one when running the post process of SDFT like DOS or conductivities.
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## Geometry relaxation

### relax_method

- **Type**: Vector of string
- **Description**: The methods to do geometry optimization. The available algorithms depend on the relax_new setting.

  First element (algorithm selection):

  - cg: Conjugate gradient (CG) algorithm. Available for both relax_new = True (default, simultaneous optimization) and relax_new = False (nested optimization). See relax_new for implementation details.
  - bfgs: Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton algorithm. Only available when relax_new = False.
  - lbfgs: Limited-memory BFGS algorithm, suitable for large systems. Only available when relax_new = False.
  - cg_bfgs: Mixed method starting with CG and switching to BFGS when force convergence reaches relax_cg_thr. Only available when relax_new = False.
  - sd: Steepest descent algorithm. Only available when relax_new = False. Not recommended for production use.
  - fire: Fast Inertial Relaxation Engine method, a molecular-dynamics-based relaxation algorithm. Use by setting calculation to md and md_type to fire. Ionic velocities must be set in STRU file. See fire for details.

  Second element (BFGS variant, only when first element is bfgs):

  - 1: Traditional BFGS that updates the Hessian matrix B and then inverts it.
  - 2 or omitted: Default BFGS that directly updates the inverse Hessian (recommended).

  > Note: In the 3.10-LTS version, the type of this parameter is std::string. It can be set to "cg", "bfgs", "cg_bfgs", "bfgs_trad", "lbfgs", "sd", "fire".
- **Default**: cg 1

### relax_new

- **Type**: Boolean
- **Description**: Controls which implementation of geometry relaxation to use. At the end of 2022, a new implementation of the Conjugate Gradient (CG) method was introduced for relax and cell-relax calculations, while the old implementation was kept for backward compatibility.


  - True (default): Use the new CG implementation with the following features:
   - Simultaneous optimization of ionic positions and cell parameters (for cell-relax)
   - Line search algorithm for step size determination
   - Only CG algorithm is available (relax_method must be cg)
   - Supports advanced cell constraints: fixed_axes = "shape", "volume", "a", "b", "c", etc.
   - Supports fixed_ibrav to maintain lattice type
   - More efficient for variable-cell relaxation
   - Step size controlled by relax_scale_force

  - False: Use the old implementation with the following features:
   - Nested optimization procedure: ionic positions optimized first, then cell parameters (for cell-relax)
   - Multiple algorithms available: cg, bfgs, lbfgs, sd, cg_bfgs
   - Limited cell constraints: only fixed_axes = "volume" is supported
   - Traditional approach with separate ionic and cell optimization steps
- **Default**: True

### relax_scale_force

- **Type**: Real
- **Availability**: *Only used when relax_new set to True*
- **Description**: The paramether controls the size of the first conjugate gradient step. A smaller value means the first step along a new CG direction is smaller. This might be helpful for large systems, where it is safer to take a smaller initial step to prevent the collapse of the whole configuration.
- **Default**: 0.5

### relax_nmax

- **Type**: Integer
- **Description**: The maximal number of ionic iteration steps. If set to 0, the code performs a quick "dry run", stopping just after initialization. This is useful to check for input correctness and to have the summary printed.
- **Default**: 1 for SCF, 50 for relax and cell-relax calcualtions

### relax_cg_thr

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method = cg_bfgs*
- **Description**: When relax_method is set to cg_bfgs, a mixed algorithm of conjugate gradient (CG) and Broyden–Fletcher–Goldfarb–Shanno (BFGS) is used. The ions first move according to the CG method, then switch to the BFGS method when the maximum force on atoms is reduced below this threshold.
- **Default**: 0.5
- **Unit**: eV/Angstrom

### force_thr

- **Type**: Real
- **Description**: Threshold of the force convergence. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to force_thr_ev except for the unit, you can choose either you like.
- **Default**: 0.001
- **Unit**: Ry/Bohr (25.7112 eV/Angstrom)

### force_thr_ev

- **Type**: Real
- **Description**: Threshold of the force convergence. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to force_thr except for the unit. You may choose either you like.
- **Default**: 0.0257112
- **Unit**: eV/Angstrom (0.03889 Ry/Bohr)

### force_zero_out

- **Type**: Real
- **Description**: The atomic forces that are smaller than force_zero_out will be treated as zero.
- **Default**: 0.0
- **Unit**: eV/Angstrom

### relax_bfgs_w1

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method is bfgs or cg_bfgs*
- **Description**: Controls the Wolfe condition for the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. This parameter sets the sufficient decrease condition (c1 in Wolfe conditions). For more information, see Phys. Chem. Chem. Phys., 2000, 2, 2177.
- **Default**: 0.01

### relax_bfgs_w2

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method is bfgs or cg_bfgs*
- **Description**: Controls the Wolfe condition for the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. This parameter sets the curvature condition (c2 in Wolfe conditions). For more information, see Phys. Chem. Chem. Phys., 2000, 2, 2177.
- **Default**: 0.5

### relax_bfgs_rmax

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method is bfgs or cg_bfgs*
- **Description**: Maximum allowed total displacement of all atoms during geometry optimization. The sum of atomic displacements can increase during optimization steps but cannot exceed this value.
- **Default**: 0.8
- **Unit**: Bohr

### relax_bfgs_rmin

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method = bfgs 1 (traditional BFGS)*
- **Description**: Minimum allowed total displacement of all atoms. When the total atomic displacement falls below this value and force convergence is not achieved, the calculation will terminate. Note: This parameter is not used in the default BFGS algorithm (relax_method = bfgs 2 or bfgs).
- **Default**: 1e-5
- **Unit**: Bohr

### relax_bfgs_init

- **Type**: Real
- **Availability**: *Only used when relax_new = False and relax_method is bfgs or cg_bfgs*
- **Description**: Initial total displacement of all atoms in the first BFGS step. This sets the scale for the initial movement.
- **Default**: 0.5
- **Unit**: Bohr

### stress_thr

- **Type**: Real
- **Description**: The threshold of the stress convergence. The threshold is compared with the largest component of the stress tensor.
- **Default**: 0.5
- **Unit**: kbar

### press1

- **Type**: Real
- **Description**: The external pressures along three axes. Positive input value is taken as compressive stress.
- **Default**: 0
- **Unit**: kbar

### press2

- **Type**: Real
- **Description**: The external pressures along three axes. Positive input value is taken as compressive stress.
- **Default**: 0
- **Unit**: kbar

### press3

- **Type**: Real
- **Description**: The external pressures along three axes. Positive input value is taken as compressive stress.
- **Default**: 0
- **Unit**: kbar

### fixed_axes

- **Type**: String
- **Availability**: *Only used when calculation is set to cell-relax*
- **Description**: Specifies which cell degrees of freedom are fixed during variable-cell relaxation. The available options depend on the relax_new setting:

  When relax_new = True (default), all options are available:

  - None: Default; all cell parameters can relax freely
  - volume: Relaxation with fixed volume (allows shape changes)
  - shape: Fix shape but allow volume changes (hydrostatic pressure only)
  - a: Fix the a-axis lattice vector during relaxation
  - b: Fix the b-axis lattice vector during relaxation
  - c: Fix the c-axis lattice vector during relaxation
  - ab: Fix both a and b axes during relaxation
  - ac: Fix both a and c axes during relaxation
  - bc: Fix both b and c axes during relaxation

  When relax_new = False, all options are now available:

  - None: Default; all cell parameters can relax freely
  - volume: Relaxation with fixed volume (allows shape changes). Volume is preserved by rescaling the lattice after each update.
  - shape: Fix shape but allow volume changes (hydrostatic pressure only). Stress tensor is replaced with isotropic pressure.
  - a, b, c, ab, ac, bc: Fix specific lattice vectors. Gradients for fixed vectors are set to zero.

  > Note: For VASP users, see the ISIF correspondence table in the geometry optimization documentation. Both implementations now support all constraint types.
- **Default**: None

### fixed_ibrav

- **Type**: Boolean
- **Availability**: *Can be used with both relax_new = True and relax_new = False. A specific latname must be provided.*
- **Description**: - True: the lattice type will be preserved during relaxation. The lattice vectors are reconstructed to match the specified Bravais lattice type after each update.
  - False: No restrictions are exerted during relaxation in terms of lattice type

  > Note: it is possible to use fixed_ibrav with fixed_axes, but please make sure you know what you are doing. For example, if we are doing relaxation of a simple cubic lattice (latname = "sc"), and we use fixed_ibrav along with fixed_axes = "volume", then the cell is never allowed to move and as a result, the relaxation never converges. When both are used, fixed_ibrav is applied first, then fixed_axes = "volume" rescaling is applied.
- **Default**: False

### fixed_atoms

- **Type**: Boolean
- **Description**: - True: The direct coordinates of atoms will be preserved during variable-cell relaxation.
  - False: No restrictions are exerted on positions of all atoms. However, users can still fix certain components of certain atoms by using the m keyword in STRU file. For the latter option, check the end of this instruction.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## Output information

### out_freq_ion

- **Type**: Integer
- **Description**: Controls the output interval in ionic steps. When set to a positive integer, information such as charge density, local potential, electrostatic potential, Hamiltonian matrix, overlap matrix, density matrix, and Mulliken population analysis is printed every n ionic steps.

  > Note: In RT-TDDFT calculations, this parameter is inactive; output frequency is instead controlled by out_freq_td.
- **Default**: 0

### out_freq_td

- **Type**: Integer
- **Description**: Controls the output interval in completed electronic evolution steps during RT-TDDFT calculations. When set to a positive integer n, detailed information (see out_freq_ion) is printed every n electron time-evolution steps (i.e., every STEP OF ELECTRON EVOLVE). For example, if you wish to output information once per ionic step, you should set out_freq_td equal to estep_per_md, since one ionic step corresponds to estep_per_md electronic evolution steps.

  > Note: This parameter is only active in RT-TDDFT mode (esolver_type = tddft). It has no effect in ground-state calculations.
- **Default**: 0

### out_freq_elec

- **Type**: Integer
- **Description**: Output the charge density (only binary format, controlled by out_chg), wavefunction (controlled by out_wfc_pw) per out_freq_elec electronic iterations. Note that they are always output when converged or reach the maximum iterations scf_nmax.
- **Default**: scf_nmax

### out_chg

- **Type**: Integer \[Integer\](optional)
- **Description**: The first integer controls whether to output the charge density on real space grids:
  - 1: Output the charge density (in Bohr^-3) on real space grids into the density files in the folder OUT.{suffix} too, which can be read in NSCF calculation.

  In molecular dynamics simulations, the output frequency is controlled by out_freq_ion.

  > Note: In the 3.10-LTS version, the file names are SPIN1_CHG.cube and SPIN1_CHG_INI.cube, etc.
- **Default**: 0 3

### out_pot

- **Type**: Integer \[Integer\](optional)
- **Description**: - 1: Output the total local potential (i.e., local pseudopotential + Hartree potential + XC potential + external electric field (if exists) + dipole correction potential (if exists) + ...) on real space grids (in Ry) into files in the folder OUT.{suffix}. The files are named as:
   - nspin = 1: pots1.cube;
   - nspin = 2: pots1.cube and pots2.cube;
   - nspin = 4: pots1.cube, pots2.cube, pots3.cube, and pots4.cube
  - 2: Output the electrostatic potential on real space grids into OUT.{suffix}/pot_es.cube. The Python script named tools/average_pot/aveElecStatPot.py can be used to calculate the average electrostatic potential along the z-axis and outputs it into ElecStaticPot_AVE. Please note that the total local potential refers to the local component of the self-consistent potential, excluding the non-local pseudopotential. The distinction between the local potential and the electrostatic potential is as follows: local potential = electrostatic potential + XC potential.
  - 3: Apart from 1, also output the total local potential of the initial charge density. The files are named as:
   - nspin = 1: pots1_ini.cube;
   - nspin = 2: pots1_ini.cube and pots2_ini.cube;
   - nspin = 4: pots1_ini.cube, pots2_ini.cube, pots3_ini.cube, and pots4_ini.cube

  The optional second integer controls the output precision. If not provided, the default precision is 8.

  In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.

  > Note: In the 3.10-LTS version, the file names are SPIN1_POT.cube and SPIN1_POT_INI.cube, etc.
- **Default**: 0

### out_dmk

- **Type**: Boolean \[Integer\](optional)
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to output the density matrix for each k-point into files in the folder OUT.${suffix}. The files are named as:
  - For gamma only case:
   - nspin = 1 and 4: dm_nao.csr;
   - nspin = 2: dms1_nao.csr and dms2_nao.csr for the two spin channels.
  - For multi-k points case:
   - nspin = 1 and 4: dmk1_nao.csr, dmk2_nao.csr, ...;
   - nspin = 2: dmk1s1_nao.csr... and dmk1s2_nao.csr... for the two spin channels.

  > Note: In the 3.10-LTS version, the parameter is named out_dm and the file names are SPIN1_DM and SPIN2_DM, etc.
- **Default**: False

### out_dmr

- **Type**: Boolean \[Integer\](optional)
- **Availability**: *Numerical atomic orbital basis (multi-k points)*
- **Description**: Whether to output the density matrix with Bravias lattice vector R index into files in the folder OUT.${suffix}. The files are named as dmr{s}{spin index}{g}{geometry index}{_nao} + {".csr"}. Here, 's' refers to spin, where s1 means spin up channel while s2 means spin down channel, and the sparse matrix format 'csr' is mentioned in out_mat_hs2. Finally, if out_app_flag is set to false, the file name contains the optional 'g' index for each ionic step that may have different geometries, and if out_app_flag is set to true, the density matrix with respect to Bravias lattice vector R accumulates during ionic steps:
  - nspin = 1: dmrs1_nao.csr;
  - nspin = 2: dmrs1_nao.csr and dmrs2_nao.csr for the two spin channels.

  > Note: In the 3.10-LTS version, the parameter is named out_dm1, and the file names are data-DMR-sparse_SPIN0.csr and data-DMR-sparse_SPIN1.csr, etc.
- **Default**: False

### out_wfc_pw

- **Type**: Integer
- **Availability**: *Output electronic wave functions in plane wave basis, or transform the real-space electronic wave function into plane wave basis (see get_wf option in calculation with NAO basis)*
- **Description**: Whether to output the electronic wavefunction coefficients into files and store them in the folder OUT.${suffix}. The files are named as wf{k}{k-point index}{s}{spin index}{g}{geometry index}{e}{electronic iteration index}{_pw} + {".txt"/".dat"}. Here, the s index refers to spin but the label will not show up for non-spin-polarized calculations, where s1 means spin up channel while s2 means spin down channel, and s4 refers to spinor wave functions that contains both spin channels with spin-orbital coupling or noncollinear calculations enabled. For scf or nscf calculations, g index will not appear, but the g index appears for geometry relaxation and molecular dynamics, where one can use the out_freq_ion command to control. To print out the electroinc wave functions every few SCF iterations, use the out_freq_elec command and the e index will appear in the file name.
  - 0: no output
  - 1: (txt format)
   - non-gamma-only with nspin=1: wfk1_pw.txt, wfk2_pw.txt, ...;
   - non-gamma-only with nspin=2: wfk1s1_pw.txt, wfk1s2_pw.txt, wfk2s1_pw.txt, wfk2s2_pw.txt, ...;
   - non-gamma-only with nspin=4: wfk1s4_pw.txt, wfk2s4_pw.txt, ...;
  - 2: (binary format)
   - non-gamma-only with nspin=1: wfk1_pw.dat, wfk2_pw.dat, ...;
   - non-gamma-only with nspin=2: wfk1s1_pw.dat, wfk1s2_pw.dat, wfk2s1_pw.dat, wfk2s2_pw.dat, ...;
   - non-gamma-only with nspin=4: wfk1s4_pw.dat, wfk2s4_pw.dat, ...;

  > Note: In the 3.10-LTS version, the file names are WAVEFUNC1.dat, WAVEFUNC2.dat, etc.
- **Default**: 0

### out_wfc_lcao

- **Type**: Integer
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to output the electronic wavefunction coefficients into files and store them in the folder OUT.${suffix}. The files are named as wf{s}{spin index}{k(optional)}{k-point index}{g(optional)}{geometry index1}{_nao} + {".txt"/".dat"}. Here, 's' refers to spin, where s1 means spin up channel while s2 means spin down channel, and 's12' refer to spinor wave functions that contains both spin channels with spin-orbital coupling or noncollinear calculations enabled. In addition, if 'gamma_only' is set to 0, then the optinoal k-point sampling index appears with the k-point index attached to the electronic wave function file names. Finally, if out_app_flag is set to false, the file name contains the optional 'g' index for each ionic step that may have different geometries, and if out_app_flag is set to true, the wave functions accumulate during ionic steps. If the out_app_flag is set to false, a new folder named WFC will be created, and the wave function files will be saved into it.
  - 0: no output
  - 1: (txt format)
   - gamma-only: wfs1_nao.txt or wfs2_nao.txt, ...;
   - non-gamma-only: wfs1k1_nao.txt or wfs1k2_nao.txt, ...;
  - 2: (binary format)
   - gamma-only: wfs1_nao.dat or wfs2_nao.dat, ...;
   - non-gamma-only: wfs1k1_nao.dat or wfs1k2_nao.dat, ....

  The corresponding sequence of the orbitals can be seen in Basis Set.

  Also controled by out_freq_ion and out_app_flag.

  > Note: In the 3.10-LTS version, the file names are WFC_NAO_GAMMA1_ION1.txt and WFC_NAO_K1_ION1.txt, etc.
- **Default**: 0

### out_dos

- **Type**: Integer
- **Description**: Whether to output the density of states (DOS). For more information, refer to the dos.md.
  - 0: no output
  - 1: output the density of states (DOS)
   - nspin=1 or 4: doss1g{geom}_{basis}.txt, where geom is the geometry index when cell changes or ions move while basis is either pw or nao.
   - nspin=2: doss1g{geom}_{basis}.txt and doss2g{geom}_{basis}.txt for two spin channles.
  - 2: (LCAO) output the density of states (DOS) and the projected density of states (PDOS)
  - 3: output the Fermi surface file (fermi.bxsf) in BXSF format that can be visualized by XCrySDen
- **Default**: 0

### out_ldos

- **Type**: Integer \[Integer\](optional)
- **Description**: Whether to output the local density of states (LDOS), optionally output precision can be set by a second parameter, default is 3.
  - 0: no output
  - 1: output the partial charge density for given bias (controlled by stm_bias) in cube file format, which can be used to plot scanning tunneling spectroscopys to mimick STM images using the Python script plot.py.
  - 2: output LDOS along a line in real space (controlled by ldos_line). Parameters used to control DOS output are also valid for LDOS.
  - 3: output both two LDOS modes above.
- **Default**: 0

### out_band

- **Type**: Boolean \[Integer\](optional)
- **Description**: Whether to output the eigenvalues of the Hamiltonian matrix (in eV) into the running log during electronic iterations and into a file at the end of calculations. The former can be used with the 'out_freq_elec' parameter while the latter option allows the output precision to be set via a second parameter, with a default value of 8. The output file names are:
   - nspin = 1 or 4: eig.txt;
   - nspin = 2: eigs1.txt and eigs2.txt;
   - For more information, refer to the band.md
- **Default**: False

### out_proj_band

- **Type**: Boolean
- **Description**: Whether to output the projected band structure. For more information, refer to the band.md
- **Default**: False

### out_stru

- **Type**: Boolean
- **Description**: Whether to output structure files per ionic step in geometry relaxation calculations into OUT.{istep}_D, where ${istep} is the ionic step.
- **Default**: False

### out_level

- **Type**: String
- **Description**: Control the output level of information in OUT.{calculation}.log.
  - ie: electronic iteration level, which prints useful information for electronic iterations;
  - i: geometry relaxation level, which prints some information for geometry relaxations additionally;
  - m: molecular dynamics level, which does not print some information for simplicity.
- **Default**: ie

### out_mat_hs

- **Type**: Boolean \[Integer\](optional)
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to print the upper triangular part of the Hamiltonian matrices and overlap matrices for each k-point into files in the directory OUT.${suffix}. The second number controls precision. For more information, please refer to hs_matrix.md. Also controled by out_freq_ion and out_app_flag.
  - For gamma only case:
   - nspin = 1: hks1_nao.txt for the Hamiltonian matrix and sks1_nao.txt for the overlap matrix;
   - nspin = 2: hks1_nao.txt and hks2_nao.txt for the Hamiltonian matrix and sks1_nao.txt for the overlap matrix. Note that the code will not output sks2_nao.txt because it is the same as sks1_nao.txt;
   - nspin = 4: hks12_nao.txt for the Hamiltonian matrix and sks12_nao.txt for the overlap matrix.
  - For multi-k points case:
   - nspin = 1: hks1k1_nao.txt for the Hamiltonian matrix at the 1st k-point, and sks1k1_nao.txt for the overlap matrix for the 1st k-point, ...;
   - nspin = 2: hks1k1_nao.txt and hks2k1_nao.txt for the two spin channels of the Hamiltonian matrix at the 1st k-point, and sks1k1_nao.txt for the overlap matrix for the 1st k-point. Note that the code will not output sks2k1_nao.txt because it is the same as sks1k1_nao.txt, ...;
   - nspin = 4: hks12k1_nao.txt for the Hamiltonian matrix at the 1st k-point, and sks12k1_nao.txt for the overlap matrix for the 1st k-point, ...;

  > Note: In the 3.10-LTS version, the file names are data-0-H and data-0-S, etc.
- **Default**: False 8
- **Unit**: Ry

### out_mat_hs2

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Whether to print files containing the Hamiltonian matrix and overlap matrix into files in the directory OUT.${suffix}. For more information, please refer to hs_matrix.md.

  > Note: In the 3.10-LTS version, the file names are data-HR-sparse_SPIN0.csr and data-SR-sparse_SPIN0.csr, etc.
- **Default**: False
- **Unit**: Ry

### out_mat_tk

- **Type**: Boolean \[Integer\](optional)
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to print the upper triangular part of the kinetic matrices for each k-point into OUT.${suffix}/tks1ki_nao.txt, where i is the index of k points. One may optionally provide a second parameter to specify the precision.

  > Note: In the 3.10-LTS version, the file names are data-TR-sparse_SPIN0.csr, etc.
- **Default**: False [8]
- **Unit**: Ry

### out_mat_r

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Whether to print the matrix representation of the position matrix into a file named rr.csr in the directory OUT.${suffix}. If calculation is set to get_s, the position matrix can be obtained without scf iterations. For more information, please refer to position_matrix.md.

  > Note: In the 3.10-LTS version, the file name is data-rR-sparse.csr.
- **Default**: False
- **Unit**: Bohr

### out_mat_t

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Generate files containing the kinetic energy matrix. The format will be the same as the Hamiltonian matrix and overlap matrix as mentioned in out_mat_hs2. The name of the files will be trs1_nao.csr and so on. Also controled by out_freq_ion and out_app_flag.

  > Note: In the 3.10-LTS version, the file name is data-TR-sparse_SPIN0.csr.
- **Default**: False
- **Unit**: Ry

### out_mat_dh

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Whether to print files containing the derivatives of the Hamiltonian matrix. The format will be the same as the Hamiltonian matrix and overlap matrix as mentioned in out_mat_hs2. The name of the files will be dhrxs1_nao.csr, dhrys1_nao.csr, dhrzs1_nao.csr and so on. Also controled by out_freq_ion and out_app_flag.

  > Note: In the 3.10-LTS version, the file name is data-dHRx-sparse_SPIN0.csr and so on.
- **Default**: False
- **Unit**: Ry/Bohr

### out_mat_ds

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Whether to print files containing the derivatives of the overlap matrix. The format will be the same as the overlap matrix as mentioned in out_mat_dh. The name of the files will be dsrxs1.csr and so on. Also controled by out_freq_ion and out_app_flag. This feature can be used with calculation get_s.

  > Note: In the 3.10-LTS version, the file name is data-dSRx-sparse_SPIN0.csr and so on.
- **Default**: False
- **Unit**: Ry/Bohr

### out_mat_xc

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital (NAO) and NAO-in-PW basis*
- **Description**: Whether to print the upper triangular part of the exchange-correlation matrices in Kohn-Sham orbital representation: for each k point into files in the directory OUT.i_nao.txt, where {suffix}/vxc_out.dat. If EXX is calculated, the local and EXX part of band energy will also be printed in OUT.{suffix}/vxc_exx_out.dat, respectively. All the vxc_out.dat files contains 3 integers (nk, nspin, nband) followed by nk*nspin*nband lines of energy Hartree and eV.

  > Note: In the 3.10-LTS version, the file name is k-$k-Vxc and so on.
- **Default**: False
- **Unit**: Ry

### out_mat_xc2

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital (NAO) basis*
- **Description**: Whether to print the exchange-correlation matrices in numerical orbital representation: in CSR format in the directory OUT.s.

  > Note: In the 3.10-LTS version, the file name is Vxc_R_spin$s and so on.
- **Default**: False
- **Unit**: Ry

### out_mat_l

- **Type**: Boolean \[Integer\](optional)
- **Availability**: *Numerical atomic orbital (NAO) basis*
- **Description**: Whether to print the expectation value of the angular momentum operator , , and in the basis of the localized atomic orbitals. The files are named OUT.{suffix}_Lx.dat, OUT.{suffix}_Ly.dat, and OUT.{suffix}_Lz.dat. The second integer controls the precision of the output.
- **Default**: False 8

### out_xc_r

- **Type**: Integer \[Integer\](optional)
- **Description**: The first integer controls whether to output the exchange-correlation (in Bohr^-3) on real space grids using Libxc to folder OUT.${suffix}:
  - 0: rho, amag, sigma, exc
  - 1: vrho, vsigma
  - 2: v2rho2, v2rhosigma, v2sigma2
  - 3: v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3
  - 4: v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 The meaning of the files is presented in Libxc

  The second integer controls the precision of the charge density output, if not given, will use 3 as default.

  The circle order of the charge density on real space grids is: x is the outer loop, then y and finally z (z is moving fastest).
- **Default**: -1 3

### out_eband_terms

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to print the band energy terms separately in the file OUT.{term}_out.dat. The terms include the kinetic, pseudopotential (local + nonlocal), Hartree and exchange-correlation (including exact exchange if calculated).
- **Default**: False

### out_mul

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to print the Mulliken population analysis result into OUT.${suffix}/mulliken.txt. In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.
- **Default**: False

### out_app_flag

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis (not gamma-only algorithm)*
- **Description**: Whether to output r(R), H(R), S(R), T(R), dH(R), dS(R), and wfc matrices in an append manner during molecular dynamics calculations. Check input parameters out_mat_r, out_mat_hs2, out_mat_t, out_mat_dh, out_mat_hs and out_wfc_lcao for more information.
- **Default**: true

### out_ndigits

- **Type**: Integer
- **Availability**: *out_mat_hs 1 case presently.*
- **Description**: Controls the length of decimal part of output data, such as charge density, Hamiltonian matrix, Overlap matrix and so on.
- **Default**: 8

### out_element_info

- **Type**: Boolean
- **Description**: Whether to print element information into files in the directory OUT.{element_label}, including pseudopotential and orbital information of the element (in atomic Ryberg units).
- **Default**: False

### restart_save

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Whether to save charge density files per ionic step, which are used to restart calculations. According to the value of read_file_dir:
  - auto: These files are saved in folder OUT.{read_file_dir}/restart/.

  If EXX(exact exchange) is calculated (i.e. dft_fuctional==hse/hf/pbe0/scan0 or rpa==True), the Hexx(R) files for each processor will also be saved in the above folder, which can be read in EXX calculation with restart_load==True.
- **Default**: False

### rpa

- **Type**: Boolean
- **Description**: Generate output files used in rpa calculations.

  > Note: If symmetry is set to 1, additional files containing the necessary information for exploiting symmetry in the subsequent rpa calculation will be output: irreducible_sector.txt, symrot_k.txt and symrot_R.txt.
- **Default**: False

### out_pchg

- **Type**: String
- **Availability**: *For both PW and LCAO. When basis_type = lcao, used when calculation = get_pchg.*
- **Description**: Specifies the electronic states to calculate the charge densities with state index for, using a space-separated string of 0s and 1s. Each digit in the string corresponds to a state, starting from the first state. A 1 indicates that the charge density should be calculated for that state, while a 0 means the state will be ignored. The parameter allows a compact and flexible notation (similar to ocp_set), for example the syntax 1 4*0 5*1 0 is used to denote the selection of states: 1 means calculate for the first state, 4*0 skips the next four states, 5*1 means calculate for the following five states, and the final 0 skips the next state. It's essential that the total count of states does not exceed the total number of states (nbands); otherwise, it results in an error, and the process exits. The input string must contain only numbers and the asterisk (*) for repetition, ensuring correct format and intention of state selection. The outputs comprise multiple .cube files following the naming convention pchgi[state]s[spin]k[kpoint].cube.
- **Default**: none

### out_wfc_norm

- **Type**: String
- **Availability**: *For both PW and LCAO. When basis_type = lcao, used when calculation = get_wf.*
- **Description**: Specifies the electronic states to calculate the real-space wave function modulus (norm, or known as the envelope function) with state index. The syntax and state selection rules are identical to out_pchg, but the output is the norm of the wave function. The outputs comprise multiple .cube files following the naming convention wfi[state]s[spin]k[kpoint].cube.
- **Default**: none

### out_wfc_re_im

- **Type**: String
- **Availability**: *For both PW and LCAO. When basis_type = lcao, used when calculation = get_wf.*
- **Description**: Specifies the electronic states to calculate the real and imaginary parts of the wave function with state index. The syntax and state selection rules are identical to out_pchg, but the output contains both the real and imaginary components of the wave function. The outputs comprise multiple .cube files following the naming convention wfi[state]s[spin]k[kpoint][re/im].cube.
- **Default**: none

### if_separate_k

- **Type**: Boolean
- **Availability**: *For both PW and LCAO. When basis_type = pw, used if out_pchg is set. When basis_type = lcao, used only when calculation = get_pchg and gamma_only = 0.*
- **Description**: Specifies whether to write the partial charge densities for all k-points to individual files or merge them. Warning: Enabling symmetry may produce unwanted results due to reduced k-point weights and symmetry operations in real space. Therefore when calculating partial charge densities, if you are not sure what you want exactly, it is strongly recommended to set symmetry = -1. It is noteworthy that your symmetry setting should remain the same as that in the SCF procedure.
- **Default**: false

### out_elf

- **Type**: Integer \[Integer\](optional)
- **Availability**: *Only for Kohn-Sham DFT and Orbital Free DFT.*
- **Description**: Whether to output the electron localization function (ELF) in the folder `OUT.${suffix}`. The files are named as
  - nspin = 1:
    - elf.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i}|^2} - \frac{|\nabla\rho|^2}{8\rho}}{\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}}$;
  - nspin = 2:
    - elf1.cube, elf2.cube: ${\rm{ELF}}_\sigma = \frac{1}{1+\chi_\sigma^2}$, $\chi_\sigma = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i,\sigma}|^2} - \frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}$;
    - elf.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i,\sigma}{f_i |\nabla\psi_{i,\sigma}|^2} - \sum_{\sigma}{\frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}}{\sum_{\sigma}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}}$;
  - nspin = 4 (noncollinear):
    - elf.cube: ELF for total charge density, ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i}|^2} - \frac{|\nabla\rho|^2}{8\rho}}{\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}}$

  The second integer controls the precision of the kinetic energy density output, if not given, will use 3 as default. For purpose restarting from this file and other high-precision involved calculation, recommend to use 10.

  In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.
- **Default**: 0 3

### out_spillage

- **Type**: Integer
- **Availability**: *Only for Kohn-Sham DFT with plane-wave basis.*
- **Description**: This output is only intentively needed by the ABACUS numerical atomic orbital generation workflow. This parameter is used to control whether to output the overlap integrals between truncated spherical Bessel functions (TSBFs) and plane-wave basis expanded wavefunctions (named as OVERLAP_Q), and between TSBFs (named as OVERLAP_Sq), also their first order derivatives. The output files are named starting with orb_matrix. A value of 2 would enable the output.
- **Default**: 0

### out_alllog

- **Type**: Boolean
- **Description**: Whether to print information into individual logs from all ranks in an MPI run.
  - True: Information from each rank will be written into individual files named OUT.{calculation}_{suffix}/running_${calculation}.log.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## Density of states

### dos_edelta_ev

- **Type**: Real
- **Description**: The step size in writing Density of States (DOS)
- **Default**: 0.01
- **Unit**: eV

### dos_sigma

- **Type**: Real
- **Description**: The width of the Gaussian factor when obtaining smeared Density of States (DOS)
- **Default**: 0.07
- **Unit**: eV

### dos_scale

- **Type**: Real
- **Description**: Defines the energy range of DOS output as (emax-emin)*(1+dos_scale), centered at (emax+emin)/2. This parameter will be used when dos_emin and dos_emax are not set.
- **Default**: 0.01
- **Unit**: eV

### dos_emin_ev

- **Type**: Real
- **Description**: The minimal range for Density of States (DOS)
  - If set, "dos_scale" will be ignored.
- **Default**: Minimal eigenenergy of
- **Unit**: eV

### dos_emax_ev

- **Type**: Real
- **Description**: The maximal range for Density of States (DOS)
  - If set, "dos_scale" will be ignored.
- **Default**: Maximal eigenenergy of
- **Unit**: eV

### dos_nche

- **Type**: Integer
- **Description**: The order of Chebyshev expansions when using Stochastic Density Functional Theory (SDFT) to calculate DOS.
- **Default**: 100

### stm_bias

- **Type**: Real Real(optional) Integer(optional)
- **Description**: The bias voltage used to calculate local density of states to simulate scanning tunneling microscope, see details in out_ldos. When using three parameters:

  - The first parameter specifies the initial bias voltage value.
  - The second parameter defines the voltage increment (step size between consecutive bias values).
  - The third parameter determines the total number of voltage points
- **Default**: 1.0
- **Unit**: V

### ldos_line

- **Type**: Real*6 Integer(optional)
- **Description**: Specify the path of the three-dimensional space and display LDOS in the form of a two-dimensional color chart, see details in out_ldos. The first three paramenters are the direct coordinates of the start point, the next three paramenters are the direct coordinates of the end point, and the final one is the number of points along the path, whose default is 100.
- **Default**: 0.0 0.0 0.0 0.0 0.0 1.0 100

[back to top](#full-list-of-input-keywords)

## NAOs

### bessel_nao_ecut

- **Type**: String
- **Description**: "Energy cutoff" (in Ry) of spherical Bessel functions. The number of spherical Bessel functions that constitute the radial parts of NAOs is determined by sqrt(bessel_nao_ecut)*bessel_nao_rcut/.
- **Default**: ecutwfc

### bessel_nao_tolerence

- **Type**: Real
- **Description**: Tolerance when searching for the zeros of spherical Bessel functions.
- **Default**: 1.0e-12

### bessel_nao_rcut

- **Type**: Vector of Real (N values)
- **Description**: Cutoff radius (in Bohr) and the common node of spherical Bessel functions used to construct the NAOs.
- **Default**: 6.0

### bessel_nao_smooth

- **Type**: Boolean
- **Description**: If True, NAOs will be smoothed near the cutoff radius. See bessel_nao_rcut and bessel_nao_sigma for parameters.
- **Default**: True

### bessel_nao_sigma

- **Type**: Real
- **Description**: Smoothing range (in Bohr). See also bessel_nao_smooth.
- **Default**: 0.1

[back to top](#full-list-of-input-keywords)

## DeePKS

### deepks_out_labels

- **Type**: Integer
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Print labels and descriptors for DeePKS in OUT.${suffix}. The names of these files start with "deepks".
  - 0 : No output.
  - 1 : Output intermediate files needed during DeePKS training.
  - 2 : Output target labels for label preperation. The label files are named as deepks_&lt;property&gt;.npy or deepks_&lt;property&gt;.csr, where the units and formats are the same as label files &lt;property&gt;.npy or &lt;property&gt;.csr required for training, except that the first dimension (nframes) is excluded. System structrue files are also given in deepks_atom.npy and deepks_box.npy in the unit of Bohr, which means lattice_constant should be set to 1 when training.

  > Note: When deepks_out_labels equals 1, the path of a numerical descriptor (an orb file) is needed to be specified under the NUMERICAL_DESCRIPTOR tag in the STRU file. This is not needed when deepks_out_labels equals 2.
- **Default**: 0

### deepks_out_freq_elec

- **Type**: Integer
- **Availability**: *Numerical atomic orbital basis*
- **Description**: When deepks_out_freq_elec is greater than 0, print labels and descriptors for DeePKS in OUT.${suffix}/DeePKS_Labels_Elec per deepks_out_freq_elec electronic iterations, with suffix _e* to distinguish different steps. Often used with deepks_out_labels equals 1.
- **Default**: 0

### deepks_out_base

- **Type**: String
- **Availability**: *Numerical atomic orbital basis and deepks_out_freq_elec is greater than 0*
- **Description**: Print labels and descriptors calculated by base functional ( determined by deepks_out_base ) and target functional ( determined by dft_functional ) for DeePKS in per deepks_out_freq_elec electronic iterations. The SCF process, labels and descriptors output of the target functional are all consistent with those when the target functional is used alone. The only additional output under this configuration is the labels of the base functional. Often used with deepks_out_labels equals 1.
- **Default**: None

### deepks_scf

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis*
- **Description**: perform self-consistent field iteration in DeePKS method

  > Note: A trained, traced model file is needed.
- **Default**: False

### deepks_equiv

- **Type**: Boolean
- **Availability**: *Numerical atomic orbital basis*
- **Description**: whether to use equivariant version of DeePKS

  > Note: The equivariant version of DeePKS-kit is still under development, so this feature is currently only intended for internal usage.
- **Default**: False

### deepks_model

- **Type**: String
- **Availability**: *Numerical atomic orbital basis and deepks_scf is true*
- **Description**: the path of the trained, traced neural network model file generated by deepks-kit
- **Default**: None

### bessel_descriptor_lmax

- **Type**: Integer
- **Availability**: *gen_bessel calculation*
- **Description**: the maximum angular momentum of the Bessel functions generated as the projectors in DeePKS - NOte: To generate such projectors, set calculation type to gen_bessel in ABACUS. See also calculation.
- **Default**: 2

### bessel_descriptor_ecut

- **Type**: String
- **Availability**: *gen_bessel calculation*
- **Description**: energy cutoff of Bessel functions
- **Default**: same as ecutwfc
- **Unit**: Ry

### bessel_descriptor_tolerence

- **Type**: Real
- **Availability**: *gen_bessel calculation*
- **Description**: tolerance for searching the zeros of Bessel functions
- **Default**: 1.0e-12

### bessel_descriptor_rcut

- **Type**: Real
- **Availability**: *gen_bessel calculation*
- **Description**: cutoff radius of Bessel functions
- **Default**: 6.0
- **Unit**: Bohr

### bessel_descriptor_smooth

- **Type**: Boolean
- **Availability**: *gen_bessel calculation*
- **Description**: smooth the Bessel functions at radius cutoff
- **Default**: False

### bessel_descriptor_sigma

- **Type**: Real
- **Availability**: *gen_bessel calculation*
- **Description**: smooth parameter at the cutoff radius of projectors
- **Default**: 0.1
- **Unit**: Bohr

### deepks_bandgap

- **Type**: Integer
- **Availability**: *Numerical atomic orbital basis and deepks_scf is true*
- **Description**: include bandgap label for DeePKS training
  - 0: Don't include bandgap label
  - 1: Include target bandgap label (see deepks_band_range for more details)
  - 2: Include multiple bandgap label (see deepks_band_range for more details)
  - 3: Used for systems containing H atoms. Here HOMO is defined as the max occupation except H atoms and the bandgap label is the energy between HOMO and (HOMO + 1)
- **Default**: 0

### deepks_band_range

- **Type**: Integer*2
- **Availability**: *Numerical atomic orbital basis, deepks_scf is true, and deepks_bandgap is 1 or 2*
- **Description**: The first value should not be larger than the second one and the meaning differs in different cases below
  - deepks_bandgap is 1: Bandgap label is the energy between LUMO + deepks_band_range[0] and LUMO + deepks_band_range[1]. If not set, it will calculate energy between HOMO and LUMO states.
  - deepks_bandgap is 2: Bandgap labels are energies between HOMO and all states in range [LUMO + deepks_band_range[0], LUMO + deepks_band_range[1]] (Thus there are deepks_band_range[1] - deepks_band_range[0] + 1 bandgaps in total). If HOMO is included in the setting range, it will be ignored since it will always be zero and has no valuable messages (deepks_band_range[1] - deepks_band_range[0] bandgaps in this case). NOTICE: The set range can be greater than, less than, or include the value of HOMO. In the bandgap label, we always calculate the energy of the state in the set range minus the energy of HOMO state, so the bandgap can be negative if the state is lower than HOMO.
- **Default**: -1 0

### deepks_v_delta

- **Type**: Integer
- **Availability**: *Numerical atomic orbital basis*
- **Description**: Include V_delta/V_delta_R (Hamiltonian in k/real space) label for DeePKS training. When deepks_out_labels is true and deepks_v_delta &gt; 0 (k space), ABACUS will output deepks_hbase.npy, deepks_vdelta.npy and deepks_htot.npy(htot=hbase+vdelta). When deepks_out_labels is true and deepks_v_delta &lt; 0 (real space), ABACUS will output deepks_hrtot.csr, deepks_hrdelta.csr. Some more files output for different settings. NOTICE: To match the unit Normally used in DeePKS, the unit of Hamiltonian in k space is Hartree. However, currently in R space the unit is still Ry.
  - deepks_v_delta = 1: deepks_vdpre.npy, which is used to calculate V_delta during DeePKS training.
  - deepks_v_delta = 2: deepks_phialpha.npy and deepks_gevdm.npy, which can be used to calculate deepks_vdpre.npy. A recommanded method for memory saving.
  - deepks_v_delta = -1: deepks_vdrpre.npy, which is used to calculate V_delta_R during DeePKS training.
  - deepks_v_delta = -2: deepks_phialpha_r.npy and deepks_gevdm.npy, which can be used to calculate deepks_vdrpre.npy. A recommanded method for memory saving.
- **Default**: 0

### deepks_out_unittest

- **Type**: Boolean
- **Description**: generate files for constructing DeePKS unit test

  > Note: Not relevant when running actual calculations. When set to 1, ABACUS needs to be run with only 1 process.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## OFDFT: orbital free density functional theory

### of_kinetic

- **Type**: String
- **Availability**: *OFDFT*
- **Description**: Kinetic energy functional type:
  - tf: Thomas-Fermi (TF) functional
  - vw: von Weizsacker (vW) functional
  - tf+: TF + vW functional
  - wt: Wang-Teter (WT) functional
  - xwm: XWM functional
  - lkt: Luo-Karasiev-Trickey (LKT) functional
  - ml: Machine learning KEDF
  - mpn: MPN KEDF (automatically sets ml parameters)
  - cpn5: CPN5 KEDF (automatically sets ml parameters)
- **Default**: wt

### of_method

- **Type**: String
- **Availability**: *OFDFT*
- **Description**: The optimization method used in OFDFT.
  - cg1: Polak-Ribiere. Standard CG algorithm.
  - cg2: Hager-Zhang (generally faster than cg1).
  - tn: Truncated Newton algorithm.
- **Default**: tn

### of_conv

- **Type**: String
- **Availability**: *OFDFT*
- **Description**: Criterion used to check the convergence of OFDFT.
  - energy: Total energy changes less than of_tole.
  - potential: The norm of potential is less than of_tolp.
  - both: Both energy and potential must satisfy the convergence criterion.
- **Default**: energy

### of_tole

- **Type**: Real
- **Availability**: *OFDFT*
- **Description**: Tolerance of the energy change for determining the convergence.
- **Default**: 2e-6
- **Unit**: Ry

### of_tolp

- **Type**: Real
- **Availability**: *OFDFT*
- **Description**: Tolerance of potential for determining the convergence.
- **Default**: 1e-5
- **Unit**: Ry

### of_tf_weight

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=tf, tf+, wt, xwm*
- **Description**: Weight of TF KEDF (kinetic energy density functional).
- **Default**: 1.0

### of_vw_weight

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=vw, tf+, wt, lkt, xwm*
- **Description**: Weight of vW KEDF (kinetic energy density functional).
- **Default**: 1.0

### of_wt_alpha

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=wt*
- **Description**: Parameter alpha of WT KEDF (kinetic energy density functional).

### of_wt_beta

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=wt*
- **Description**: Parameter beta of WT KEDF (kinetic energy density functional).

### of_wt_rho0

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=wt*
- **Description**: The average density of system.
- **Default**: 0.0
- **Unit**: Bohr^-3

### of_hold_rho0

- **Type**: Boolean
- **Availability**: *OFDFT with of_kinetic=wt*
- **Description**: Whether to fix the average density rho0.
  - True: rho0 will be fixed even if the volume of system has changed, it will be set to True automatically if of_wt_rho0 is not zero.
  - False: rho0 will change if volume of system has changed.
- **Default**: False

### of_lkt_a

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=lkt*
- **Description**: Parameter a of LKT KEDF (kinetic energy density functional).
- **Default**: 1.3

### of_xwm_rho_ref

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=xwm*
- **Description**: Reference charge density for XWM kinetic energy functional. If set to 0, the program will use average charge density.
- **Default**: 0.0

### of_xwm_kappa

- **Type**: Real
- **Availability**: *OFDFT with of_kinetic=xwm*
- **Description**: Parameter for XWM kinetic energy functional. See PHYSICAL REVIEW B 100, 205132 (2019) for optimal values.
- **Default**: 0.0

### of_read_kernel

- **Type**: Boolean
- **Availability**: *OFDFT with of_kinetic=wt*
- **Description**: Whether to read in the kernel file.
  - True: The kernel of WT KEDF (kinetic energy density functional) will be filled from the file specified by of_kernel_file.
  - False: The kernel of WT KEDF (kinetic energy density functional) will be filled from formula.
- **Default**: False

### of_kernel_file

- **Type**: String
- **Availability**: *OFDFT with of_read_kernel=True*
- **Description**: The name of WT kernel file.
- **Default**: WTkernel.txt

### of_full_pw

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Whether to use full planewaves.
  - True: Ecut will be ignored while collecting planewaves, so that all planewaves will be used in FFT.
  - False: Only use the planewaves inside ecut, the same as KSDFT.
- **Default**: True

### of_full_pw_dim

- **Type**: Integer
- **Availability**: *OFDFT with of_full_pw = True*
- **Description**: Specify the parity of FFT dimensions.
  - 0: either odd or even.
  - 1: odd only.
  - 2: even only.

  Note: Even dimensions may cause slight errors in FFT. It should be ignorable in ofdft calculation, but it may make Cardinal B-spline interpolation unstable, so please set of_full_pw_dim = 1 if nbspline != -1.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## ML-KEDF: machine learning based kinetic energy density functional for OFDFT

### of_ml_gene_data

- **Type**: Boolean
- **Availability**: *Used only for KSDFT with plane wave basis*
- **Description**: Controls the generation of machine learning training data. When enabled, training data in .npy format will be saved in the directory OUT.${suffix}/.
- **Default**: False

### of_ml_device

- **Type**: String
- **Availability**: *OFDFT*
- **Description**: Run Neural Network on GPU or CPU.
  - cpu: CPU
  - gpu: GPU
- **Default**: cpu

### of_ml_feg

- **Type**: Integer
- **Availability**: *OFDFT*
- **Description**: The method to incorporate the Free Electron Gas (FEG) limit.
  - 0: Do not incorporate the FEG limit.
  - 1: Incorporate the FEG limit by translation.
  - 3: Incorporate the FEG limit by nonlinear transformation using softplus function.
- **Default**: 0

### of_ml_nkernel

- **Type**: Integer
- **Availability**: *OFDFT*
- **Description**: Number of kernel functions.
- **Default**: 1

### of_ml_kernel

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the type of the i-th kernel function.
  - 1: Wang-Teter kernel function.
  - 2: Modified Yukawa function, and alpha is specified by of_ml_yukawa_alpha.
  - 3: Truncated kinetic kernel (TKK), the file containing TKK is specified by of_ml_kernel_file.
- **Default**: 1

### of_ml_kernel_scaling

- **Type**: Vector of Real
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the RECIPROCAL of scaling parameter of the i-th kernel function.
- **Default**: 1.0

### of_ml_yukawa_alpha

- **Type**: Vector of Real
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the parameter alpha of i-th kernel function. ONLY used for Yukawa kernel function.
- **Default**: 1.0

### of_ml_kernel_file

- **Type**: Vector of String
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the file containing the i-th kernel function. ONLY used for TKK.
- **Default**: none

### of_ml_gamma

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Local descriptor: gamma = (rho / rho0)^(1/3).
- **Default**: False

### of_ml_p

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Semi-local descriptor: p = |nabla rho|^2 / [2 (3 pi^2)^(1/3) rho^(4/3)]^2.
- **Default**: False

### of_ml_q

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Semi-local descriptor: q = nabla^2 rho / [4 (3 pi^2)^(2/3) rho^(5/3)].
- **Default**: False

### of_ml_tanhp

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Semi-local descriptor: tanhp = tanh(chi_p * p).
- **Default**: False

### of_ml_tanhq

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: Semi-local descriptor: tanhq = tanh(chi_q * q).
- **Default**: False

### of_ml_chi_p

- **Type**: Real
- **Availability**: *OFDFT*
- **Description**: Hyperparameter chi_p: tanhp = tanh(chi_p * p).
- **Default**: 1.0

### of_ml_chi_q

- **Type**: Real
- **Availability**: *OFDFT*
- **Description**: Hyperparameter chi_q: tanhq = tanh(chi_q * q).
- **Default**: 1.0

### of_ml_gammanl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor gammanl defined by the i-th kernel function.
- **Default**: 0

### of_ml_pnl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor pnl defined by the i-th kernel function.
- **Default**: 0

### of_ml_qnl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor qnl defined by the i-th kernel function.
- **Default**: 0

### of_ml_xi

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor xi defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanhxi

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhxi defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanhxi_nl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhxi_nl defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanh_pnl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanh_pnl defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanh_qnl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanh_qnl defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanhp_nl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhp_nl defined by the i-th kernel function.
- **Default**: 0

### of_ml_tanhq_nl

- **Type**: Vector of Integer
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element controls the non-local descriptor tanhq_nl defined by the i-th kernel function.
- **Default**: 0

### of_ml_chi_xi

- **Type**: Vector of Real
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_xi of non-local descriptor tanhxi defined by the i-th kernel function.
- **Default**: 1.0

### of_ml_chi_pnl

- **Type**: Vector of Real
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_pnl of non-local descriptor tanh_pnl defined by the i-th kernel function.
- **Default**: 1.0

### of_ml_chi_qnl

- **Type**: Vector of Real
- **Availability**: *OFDFT*
- **Description**: Containing nkernel (see of_ml_nkernel) elements. The i-th element specifies the hyperparameter chi_qnl of non-local descriptor tanh_qnl defined by the i-th kernel function.
- **Default**: 1.0

### of_ml_local_test

- **Type**: Boolean
- **Availability**: *OFDFT*
- **Description**: FOR TEST. Read in the density, and output the F and Pauli potential.
- **Default**: False

### ml_exx

- **Type**: Boolean
- **Description**: Whether to use machine learning based exact exchange (ML-EXX).
- **Default**: False

[back to top](#full-list-of-input-keywords)

## TDOFDFT: time dependent orbital free density functional theory

### of_cd

- **Type**: Boolean
- **Availability**: *TDOFDFT*
- **Description**: Added the current dependent(CD) potential. (https://doi.org/10.1103/PhysRevB.98.144302)
  - True: Added the CD potential.
  - False: Not added the CD potential.
- **Default**: False

### of_mcd_alpha

- **Type**: Real
- **Availability**: *TDOFDFT*
- **Description**: The value of the parameter alpha in modified CD potential method. mCDPotential=alpha*CDPotential (proposed in paper PhysRevB.98.144302)
- **Default**: 1.0

[back to top](#full-list-of-input-keywords)

## Electric field and dipole correction

### efield_flag

- **Type**: Boolean
- **Description**: Added the electric field.
  - True: A saw-like potential simulating an electric field is added to the bare ionic potential.
  - False: Not added the electric field.
- **Default**: False

### dip_cor_flag

- **Type**: Boolean
- **Availability**: *With dip_cor_flag = True and efield_flag = True.*
- **Description**: Added a dipole correction to the bare ionic potential.
  - True: A dipole correction is also added to the bare ionic potential.
  - False: A dipole correction is not added to the bare ionic potential.

  > Note: If you do not want any electric field, the parameter efield_amp should be set to zero. This should ONLY be used in a slab geometry for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
- **Default**: False

### efield_dir

- **Type**: Integer
- **Availability**: *with efield_flag = True.*
- **Description**: The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir can set to 0, 1 or 2.
  - 0: parallel to the first reciprocal lattice vector
  - 1: parallel to the second reciprocal lattice vector
  - 2: parallel to the third reciprocal lattice vector
- **Default**: 2

### efield_pos_max

- **Type**: Real
- **Availability**: *with efield_flag = True.*
- **Description**: Position of the maximum of the saw-like potential along crystal axis efield_dir, within the unit cell, 0 &lt;= efield_pos_max &lt; 1.
- **Default**: Autoset to center of vacuum - width of vacuum / 20

### efield_pos_dec

- **Type**: Real
- **Availability**: *with efield_flag = True.*
- **Description**: Zone in the unit cell where the saw-like potential decreases, 0 &lt; efield_pos_dec &lt; 1.
- **Default**: Autoset to width of vacuum / 10

### efield_amp

- **Type**: Real
- **Availability**: *with efield_flag = True.*
- **Description**: Amplitude of the electric field. The saw-like potential increases with slope efield_amp in the region from efield_pos_max+efield_pos_dec-1) to (efield_pos_max), then decreases until (efield_pos_max+efield_pos_dec), in units of the crystal vector efield_dir.

  > Note: The change of slope of this potential must be located in the empty region, or else unphysical forces will result.
- **Default**: 0.0
- **Unit**: a.u., 1 a.u. = 51.4220632*10^10 V/m.

[back to top](#full-list-of-input-keywords)

## Gate field (compensating charge)

### gate_flag

- **Type**: Boolean
- **Description**: Controls the addition of compensating charge by a charged plate for charged cells.
  - true: A charged plate is placed at the zgate position to add compensating charge. The direction is determined by efield_dir.
  - false: No compensating charge is added.
- **Default**: false

### zgate

- **Type**: Real
- **Description**: Position of the charged plate in the unit cell
- **Default**: 0.5
- **Unit**: Unit cell size

### block

- **Type**: Boolean
- **Description**: Controls the addition of a potential barrier to prevent electron spillover.
  - true: A potential barrier is added from block_down to block_up with a height of block_height. If dip_cor_flag is set to true, efield_pos_dec is used to smoothly increase and decrease the potential barrier.
  - false: No potential barrier is added.
- **Default**: false

### block_down

- **Type**: Real
- **Description**: Lower beginning of the potential barrier
- **Default**: 0.45
- **Unit**: Unit cell size

### block_up

- **Type**: Real
- **Description**: Upper beginning of the potential barrier
- **Default**: 0.55
- **Unit**: Unit cell size

### block_height

- **Type**: Real
- **Description**: Height of the potential barrier
- **Default**: 0.1
- **Unit**: Rydberg

[back to top](#full-list-of-input-keywords)

## Exact Exchange (Common)

### exx_fock_alpha

- **Type**: Real
- **Description**: Fraction of full-ranged Fock exchange $1/r$ in range-separated hybrid functionals.
- **Default**: see hybrid_func_params

### exx_erfc_alpha

- **Type**: Real
- **Description**: Fraction of short-ranged Fock exchange $\mathrm{erfc}(\omega r)/r$ in range-separated hybrid functionals.
- **Default**: see hybrid_func_params

### exx_erfc_omega

- **Type**: Real
- **Description**: Range-separation parameter $\omega$ in the short-ranged Fock term $\mathrm{erfc}(\omega r)/r$.
- **Default**: see hybrid_func_params

### exx_separate_loop

- **Type**: Boolean
- **Description**: There are two types of iterative approaches provided by ABACUS to evaluate Fock exchange.
  - False: Start with a GGA-Loop, and then Hybrid-Loop, in which EXX Hamiltonian is updated with electronic iterations.
  - True: A two-step method is employed, i.e. in the inner iterations, density matrix is updated, while in the outer iterations, is calculated based on density matrix that converges in the inner iteration.
- **Default**: True

### exx_hybrid_step

- **Type**: Integer
- **Availability**: *exx_separate_loop==1*
- **Description**: The maximal iteration number of the outer-loop, where the Fock exchange is calculated
- **Default**: 100

### exx_mixing_beta

- **Type**: Real
- **Availability**: *exx_separate_loop==1*
- **Description**: Mixing parameter for densty matrix in each iteration of the outer-loop
- **Default**: 1.0

[back to top](#full-list-of-input-keywords)

## Exact Exchange (LCAO in PW)

### exx_fock_lambda

- **Type**: Real
- **Availability**: *basis_type==lcao_in_pw*
- **Description**: It is used to compensate for divergence points at G=0 in the evaluation of Fock exchange using lcao_in_pw method.
- **Default**: 0.3

[back to top](#full-list-of-input-keywords)

## Exact Exchange (LCAO)

### exx_pca_threshold

- **Type**: Real
- **Description**: To accelerate the evaluation of four-center integrals (), the product of atomic orbitals are expanded in the basis of auxiliary basis functions (ABF): . The size of the ABF (i.e. number of ) is reduced using principal component analysis. When a large PCA threshold is used, the number of ABF will be reduced, hence the calculation becomes faster. However, this comes at the cost of computational accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. Smaller components (less than exx_c_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_cs_inv_thr

- **Type**: Real
- **Description**: By default, the Coulomb matrix inversion required for obtaining LRI coefficients is performed using LU decomposition. However, this approach may suffer from numerical instabilities when a large set of auxiliary basis functions (ABFs) is employed. When exx_cs_inv_thr &gt; 0, the inversion is instead carried out via matrix diagonalization. Eigenvalues smaller than exx_cs_inv_thr are discarded to improve numerical stability. A relatively safe and commonly recommended value is 1e-5.
- **Default**: -1

### exx_v_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_dm_threshold

- **Type**: Real
- **Description**: The Fock exchange can be expressed as where D is the density matrix. Smaller values of the density matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_c_grad_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. is used in force. Smaller components (less than exx_c_grad_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_grad_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. is used in force. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_c_grad_r_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. is used in stress. Smaller components (less than exx_c_grad_r_threshold) of the matrix are neglected to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 1e-4.
- **Default**: 1E-4

### exx_v_grad_r_threshold

- **Type**: Real
- **Description**: See also the entry exx_pca_threshold. With the approximation , the four-center integral in Fock exchange is expressed as , where is a double-center integral. is used in force and stress. Smaller values of the V matrix can be truncated to accelerate calculation. The larger the threshold is, the faster the calculation and the lower the accuracy. A relatively safe choice of the value is 0, i.e. no truncation.
- **Default**: 1E-1

### exx_ccp_rmesh_times

- **Type**: String
- **Description**: This parameter determines how many times larger the radial mesh required for calculating Columb potential is to that of atomic orbitals. The value should be larger than 0. Reducing this value can effectively increase the speed of self-consistent calculations using hybrid functionals.

### exx_opt_orb_lmax

- **Type**: Integer
- **Availability**: *calculation==gen_opt_abfs*
- **Description**: The maximum l of the spherical Bessel functions, when the radial part of opt-ABFs are generated as linear combinations of spherical Bessel functions. A reasonable choice is 2.
- **Default**: 0

### exx_opt_orb_ecut

- **Type**: Real
- **Availability**: *calculation==gen_opt_abfs*
- **Description**: The cut-off of plane wave expansion, when the plane wave basis is used to optimize the radial ABFs. A reasonable choice is 60.
- **Default**: 0
- **Unit**: Ry

### exx_opt_orb_tolerence

- **Type**: Real
- **Availability**: *calculation==gen_opt_abfs*
- **Description**: The threshold when solving for the zeros of spherical Bessel functions. A reasonable choice is 1e-12.
- **Default**: 1E-12

### exx_real_number

- **Type**: String
- **Description**: - True: Enforce LibRI to use double data type.
  - False: Enforce LibRI to use complex data type. Setting it to True can effectively improve the speed of self-consistent calculations with hybrid functionals.
- **Default**: depends on the gamma_only option

### exx_singularity_correction

- **Type**: String
- **Description**: - spencer: see Phys. Rev. B 77, 193110 (2008).
  - revised_spencer: see Phys. Rev. Mater. 5, 013807 (2021). Set the scheme of Coulomb singularity correction.
- **Default**: default

### rpa_ccp_rmesh_times

- **Type**: Real
- **Description**: How many times larger the radial mesh required is to that of atomic orbitals in the postprocess calculation of the bare Coulomb matrix for RPA, GW, etc.
- **Default**: 10

### exx_symmetry_realspace

- **Type**: Boolean
- **Availability**: *symmetry==1 and exx calculation (dft_fuctional==hse/hf/pbe0/scan0 or rpa==True)*
- **Description**: - False: only rotate k-space density matrix D(k) from irreducible k-points to accelerate diagonalization
  - True: rotate both D(k) and Hexx(R) to accelerate both diagonalization and EXX calculation
- **Default**: True

### out_ri_cv

- **Type**: Boolean
- **Description**: Whether to output the coefficient tensor C(R) and ABFs-representation Coulomb matrix V(R) for each atom pair and cell in real space.
- **Default**: false

[back to top](#full-list-of-input-keywords)

## Exact Exchange (PW)

### exxace

- **Type**: Boolean
- **Availability**: *exx_separate_loop==True.*
- **Description**: Whether to use the ACE method (https://doi.org/10.1021/acs.jctc.6b00092) to accelerate the calculation the Fock exchange matrix. Should be set to true most of the time.
  - True: Use the ACE method to calculate the Fock exchange operator.
  - False: Use the traditional method to calculate the Fock exchange operator.
- **Default**: True

### exx_gamma_extrapolation

- **Type**: Boolean
- **Description**: Whether to use the gamma point extrapolation method to calculate the Fock exchange operator. See https://doi.org/10.1103/PhysRevB.79.205114 for details. Should be set to true most of the time.
- **Default**: True

### ecutexx

- **Type**: Real
- **Description**: The energy cutoff for EXX (Fock) exchange operator in plane wave basis calculations. Reducing ecutexx below ecutrho may significantly accelerate EXX computations. This speed improvement comes with a reduced numerical accuracy in the exchange energy calculation.
- **Default**: same as ecutrho
- **Unit**: Ry

### exx_thr_type

- **Type**: String
- **Description**: The type of threshold used to judge whether the outer loop has converged in the separate loop EXX calculation.
  - energy: use the change of exact exchange energy to judge convergence.
  - density: if the change of charge density difference between two successive outer loop iterations is seen as converged according to scf_thr, then the outer loop is seen as converged.
- **Default**: density

### exx_ene_thr

- **Type**: Real
- **Availability**: *exx_thr_type==energy*
- **Description**: The threshold for the change of exact exchange energy to judge convergence of the outer loop in the separate loop EXX calculation.
- **Default**: 1e-5
- **Unit**: Ry

[back to top](#full-list-of-input-keywords)

## Molecular dynamics

### md_type

- **Type**: String
- **Description**: Control the algorithm to integrate the equation of motion for molecular dynamics (MD), see md.md in detail.

  - fire: a MD-based relaxation algorithm, named fast inertial relaxation engine.
  - nve: NVE ensemble with velocity Verlet algorithm.
  - nvt: NVT ensemble, see md_thermostat in detail.
  - npt: Nose-Hoover style NPT ensemble, see md_pmode in detail.
  - langevin: NVT ensemble with Langevin thermostat, see md_damp in detail.
  - msst: MSST method, see msst_direction, msst_vel, msst_qmass, msst_vis, msst_tscale in detail.
- **Default**: nvt

### md_nstep

- **Type**: Integer
- **Description**: The total number of molecular dynamics steps.
- **Default**: 10

### md_dt

- **Type**: Real
- **Description**: The time step used in molecular dynamics calculations.
- **Default**: 1.0
- **Unit**: fs

### md_thermostat

- **Type**: String
- **Description**: Specify the temperature control method used in NVT ensemble.

  - nhc: Nose-Hoover chain, see md_tfreq and md_tchain in detail.
  - anderson: Anderson thermostat, see md_nraise in detail.
  - berendsen: Berendsen thermostat, see md_nraise in detail.
  - rescaling: velocity Rescaling method 1, see md_tolerance in detail.
  - rescale_v: velocity Rescaling method 2, see md_nraise in detail.
- **Default**: nhc

### md_tfirst

- **Type**: Real
- **Description**: The temperature used in molecular dynamics calculations.

  If md_tfirst is unset or less than zero, init_vel is autoset to be true. If init_vel is true, the initial temperature will be determined by the velocities read from STRU. In this case, if velocities are unspecified in STRU, the initial temperature is set to zero.

  If md_tfirst is set to a positive value and init_vel is true simultaneously, please make sure they are consistent, otherwise abacus will exit immediately.

  Note that md_tlast is only used in NVT/NPT simulations. If md_tlast is unset or less than zero, md_tlast is set to md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.
- **Default**: No default
- **Unit**: K

### md_tlast

- **Type**: Real
- **Description**: The temperature used in molecular dynamics calculations.

  If md_tfirst is unset or less than zero, init_vel is autoset to be true. If init_vel is true, the initial temperature will be determined by the velocities read from STRU. In this case, if velocities are unspecified in STRU, the initial temperature is set to zero.

  If md_tfirst is set to a positive value and init_vel is true simultaneously, please make sure they are consistent, otherwise abacus will exit immediately.

  Note that md_tlast is only used in NVT/NPT simulations. If md_tlast is unset or less than zero, md_tlast is set to md_tfirst. If md_tlast is set to be different from md_tfirst, ABACUS will automatically change the temperature from md_tfirst to md_tlast.
- **Default**: No default
- **Unit**: K

### md_prec_level

- **Type**: Integer
- **Description**: Determine the precision level of variable-cell molecular dynamics calculations.
  - 0: FFT grids do not change, only G vectors and K vectors are changed due to the change of lattice vector. This level is suitable for cases where the variation of the volume and shape is not large, and the efficiency is relatively higher.
  - 2: FFT grids change per step. This level is suitable for cases where the variation of the volume and shape is large, such as the MSST method. However, accuracy comes at the cost of efficiency.
- **Default**: 0

### md_restart

- **Type**: Boolean
- **Description**: Control whether to restart molecular dynamics calculations and time-dependent density functional theory calculations.
  - True: ABACUS will read in {md_step}, then read in the corresponding STRU_MD_suffix/STRU/ automatically. For tddft, ABACUS will also read in WFC_NAO_K${kpoint} of the last step (You need to set out_wfc_lcao=1 and out_app_flag=0 to obtain this file).
  - False: ABACUS will start molecular dynamics calculations normally from the first step.
- **Default**: False

### md_restartfreq

- **Type**: Integer
- **Description**: The output frequency of OUT.{suffix}/STRIU/, which are used to restart molecular dynamics calculations, see md_restart in detail.
- **Default**: 5

### md_dumpfreq

- **Type**: Integer
- **Description**: The output frequency of OUT.${suffix}/MD_dump in molecular dynamics calculations, which including the information of lattices and atoms.
- **Default**: 1

### dump_force

- **Type**: Boolean
- **Description**: Whether to output atomic forces into the file OUT.${suffix}/MD_dump.
- **Default**: True

### dump_vel

- **Type**: Boolean
- **Description**: Whether to output atomic velocities into the file OUT.${suffix}/MD_dump.
- **Default**: True

### dump_virial

- **Type**: Boolean
- **Description**: Whether to output lattice virials into the file OUT.${suffix}/MD_dump.
- **Default**: True

### md_seed

- **Type**: Integer
- **Description**: The random seed to initialize random numbers used in molecular dynamics calculations.
  - &lt; 0: No srand() function is called.
  - &gt;= 0: The function srand(md_seed) is called.
- **Default**: -1

### md_tfreq

- **Type**: Real
- **Description**: Control the frequency of temperature oscillations during the simulation. If it is too large, the temperature will fluctuate violently; if it is too small, the temperature will take a very long time to equilibrate with the atomic system.

  Note: It is a system-dependent empirical parameter, ranging from 1/(40*md_dt) to 1/(100*md_dt). An improper choice might lead to the failure of jobs.
- **Default**: 1/40/md_dt

### md_tchain

- **Type**: Integer
- **Description**: Number of thermostats coupled with the particles in the NVT/NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
- **Default**: 1

### md_pmode

- **Type**: String
- **Description**: Determine the precision level of variable-cell molecular dynamics calculations.
  - 0: FFT grids do not change, only G vectors and K vectors are changed due to the change of lattice vector. This level is suitable for cases where the variation of the volume and shape is not large, and the efficiency is relatively higher.
  - 2: FFT grids change per step. This level is suitable for cases where the variation of the volume and shape is large, such as the MSST method. However, accuracy comes at the cost of efficiency.
- **Default**: iso

### ref_cell_factor

- **Type**: Real
- **Description**: Construct a reference cell bigger than the initial cell. The reference cell has to be large enough so that the lattice vectors of the fluctuating cell do not exceed the reference lattice vectors during MD. Typically, 1.02 ~ 1.10 is sufficient. However, the cell fluctuations depend on the specific system and thermodynamic conditions. So users must test for a proper choice. This parameters should be used in conjunction with erf_ecut, erf_height, and erf_sigma.
- **Default**: 1.0

### md_pcouple

- **Type**: String
- **Description**: The coupled lattice vectors will scale proportionally in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
  - none: Three lattice vectors scale independently.
  - xyz: Lattice vectors x, y, and z scale proportionally.
  - xy: Lattice vectors x and y scale proportionally.
  - xz: Lattice vectors x and z scale proportionally.
  - yz: Lattice vectors y and z scale proportionally.
- **Default**: none

### md_pfirst

- **Type**: Real
- **Description**: The target pressure used in NPT ensemble simulations, the default value of md_plast is md_pfirst. If md_plast is set to be different from md_pfirst, ABACUS will automatically change the target pressure from md_pfirst to md_plast.
- **Default**: -1.0
- **Unit**: kbar

### md_plast

- **Type**: Real
- **Description**: The target pressure used in NPT ensemble simulations, the default value of md_plast is md_pfirst. If md_plast is set to be different from md_pfirst, ABACUS will automatically change the target pressure from md_pfirst to md_plast.
- **Default**: -1.0
- **Unit**: kbar

### md_pfreq

- **Type**: Real
- **Description**: The frequency of pressure oscillations during the NPT ensemble simulation. If it is too large, the pressure will fluctuate violently; if it is too small, the pressure will take a very long time to equilibrate with the atomic system.

  Note: It is a system-dependent empirical parameter. An improper choice might lead to the failure of jobs.
- **Default**: 1/400/md_dt

### md_pchain

- **Type**: Integer
- **Description**: The number of thermostats coupled with the barostat in the NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.
- **Default**: 1

### lj_rule

- **Type**: Integer
- **Description**: The Lennard-Jones potential between two atoms equals: $\sigma_k\sigma(i,j)$
- **Default**: 2

### lj_eshift

- **Type**: Boolean
- **Description**: It True, the LJ potential is shifted by a constant such that it is zero at the cut-off distance.
- **Default**: False

### lj_rcut

- **Type**: Real
- **Description**: Cut-off radius for Leonard Jones potential, beyond which the interaction will be neglected. It can be a single value, which means that all pairs of atoms types share the same cut-off radius. Otherwise, it should be a multiple-component vector, containing values, see details in lj_rule.
- **Default**: No default
- **Unit**: Angstrom

### lj_epsilon

- **Type**: Real
- **Description**: The vector representing the matrix for Leonard Jones potential. See details in lj_rule.
- **Default**: No default
- **Unit**: eV

### lj_sigma

- **Type**: Real
- **Description**: The vector representing the matrix for Leonard Jones potential. See details in lj_rule.
- **Default**: No default
- **Unit**: Angstrom

### pot_file

- **Type**: String
- **Description**: The filename of DP/NEP potential files, see md.md in detail.
- **Default**: graph.pb

### dp_rescaling

- **Type**: Real
- **Availability**: *esolver_type = dp.*
- **Description**: Rescaling factor to use a temperature-dependent DP. Energy, stress and force calculated by DP will be multiplied by this factor.
- **Default**: 1.0

### dp_fparam

- **Type**: Real
- **Availability**: *esolver_type = dp.*
- **Description**: The frame parameter for dp potential. The array size is dim_fparam, then all frames are assumed to be provided with the same fparam.
- **Default**: {}

### dp_aparam

- **Type**: Real
- **Availability**: *esolver_type = dp.*
- **Description**: The atomic parameter for dp potential. The array size can be (1) natoms x dim_aparam, then all frames are assumed to be provided with the same aparam; (2) dim_aparam, then all frames and atoms are assumed to be provided with the same aparam.
- **Default**: {}

### msst_direction

- **Type**: Integer
- **Description**: The direction of the shock wave in the MSST method.
  - 0: x direction
  - 1: y direction
  - 2: z direction
- **Default**: 2

### msst_vel

- **Type**: Real
- **Description**: The velocity of the shock wave in the MSST method.
- **Default**: 0.0
- **Unit**: Angstrom/fs

### msst_vis

- **Type**: Real
- **Description**: Artificial viscosity in the MSST method.
- **Default**: 0.0
- **Unit**: g/(mol*Angstrom*fs)

### msst_tscale

- **Type**: Real
- **Description**: The reduction percentage of the initial temperature used to compress volume in the MSST method.
- **Default**: 0.01

### msst_qmass

- **Type**: Real
- **Description**: Inertia of the extended system variable. You should set a number larger than 0.
- **Default**: No default

### md_damp

- **Type**: Real
- **Description**: The damping parameter used to add fictitious force in the Langevin method.
- **Default**: 1.0
- **Unit**: fs

### md_tolerance

- **Type**: Real
- **Description**: The temperature tolerance for velocity rescaling. Velocities are rescaled if the current and target temperature differ more than md_tolerance.
- **Default**: 100.0
- **Unit**: K

### md_nraise

- **Type**: Integer
- **Description**: - Anderson: The "collision frequency" parameter is given as 1/md_nraise.
  - Berendsen: The "rise time" parameter is given in units of the time step: tau = md_nraise*md_dt, so md_dt/tau = 1/md_nraise.
  - Rescale_v: Every md_nraise steps the current temperature is rescaled to the target temperature.
- **Default**: 1

### cal_syns

- **Type**: Boolean
- **Description**: Whether to calculate and output asynchronous overlap matrix for Hefei-NAMD interface. When enabled, calculates &lt;phi(t-1)|phi(t)&gt; by computing overlap between basis functions at atomic positions from previous time step and current time step. The overlap is calculated by shifting atom positions backward by velocity x md_dt. Output file: OUT.*/syns_nao.csr in CSR format.

  > Note: Only works with LCAO basis and molecular dynamics calculations. Requires atomic velocities. Output starts from the second MD step (istep &gt; 0).
- **Default**: False

### dmax

- **Type**: Real
- **Description**: The maximum displacement of all atoms in one step. This parameter is useful when cal_syns = True.
- **Default**: 0.01
- **Unit**: bohr

[back to top](#full-list-of-input-keywords)

## DFT+U correction

### dft_plus_u

- **Type**: Integer
- **Description**: Determines whether to calculate the plus U correction, which is especially important for correlated electrons.
  - 1: Calculate plus U correction with radius-adjustable localized projections (with parameter onsite_radius).
  - 2: Calculate plus U correction using first zeta of NAOs as projections (this is old method for testing).
  - 0: Do not calculate plus U correction.
- **Default**: 0

### dft_plus_dmft

- **Type**: Boolean
- **Availability**: *basis_type==lcao*
- **Description**: Whether to enable DFT+DMFT calculation. True: DFT+DMFT; False: standard DFT calculation.
- **Default**: False

### orbital_corr

- **Type**: Vector of Integer (n values where n is the number of atomic types)
- **Description**: Specifies which orbits need plus U correction for each atom type ( for atom type 1, 2, 3, respectively).
  - -1: The plus U correction will not be calculated for this atom.
  - 1: For p-electron orbits, the plus U correction is needed.
  - 2: For d-electron orbits, the plus U correction is needed.
  - 3: For f-electron orbits, the plus U correction is needed.
- **Default**: -1

### hubbard_u

- **Type**: Vector of Real (n values where n is the number of atomic types)
- **Description**: Specifies the Hubbard Coulomb interaction parameter U (eV) in plus U correction, which should be specified for each atom unless the Yukawa potential is used.

  > Note: Since only the simplified scheme by Duradev is implemented, the 'U' here is actually U-effective, which is given by Hubbard U minus Hund J.
- **Default**: 0.0

### yukawa_potential

- **Type**: Boolean
- **Description**: Determines whether to use the local screen Coulomb potential method to calculate the values of U and J.
  - True: hubbard_u does not need to be specified.
  - False: hubbard_u does need to be specified.
- **Default**: False

### yukawa_lambda

- **Type**: Real
- **Availability**: *DFT+U with yukawa_potential = True.*
- **Description**: The screen length of Yukawa potential. If left to default, the screen length will be calculated as an average of the entire system. It's better to stick to the default setting unless there is a very good reason.
- **Default**: Calculated on the fly.

### uramping

- **Type**: Real
- **Availability**: *DFT+U calculations with mixing_restart &gt; 0.*
- **Description**: Once uramping &gt; 0.15 eV. DFT+U calculations will start SCF with U = 0 eV, namely normal LDA/PBE calculations. Once SCF restarts when drho&lt;mixing_restart, U value will increase by uramping eV. SCF will repeat above calcuations until U values reach target defined in hubbard_u. As for uramping=1.0 eV, the recommendations of mixing_restart is around 5e-4.
- **Default**: -1.0.
- **Unit**: eV

### omc

- **Type**: Integer
- **Description**: The parameter controls the form of occupation matrix control used.
  - 0: No occupation matrix control is performed, and the onsite density matrix will be calculated from wavefunctions in each SCF step.
  - 1: The first SCF step will use an initial density matrix read from a file named initial_onsite.dm, but for later steps, the onsite density matrix will be updated.
  - 2: The same onsite density matrix from initial_onsite.dm will be used throughout the entire calculation.

  > Note: The easiest way to create initial_onsite.dm is to run a DFT+U calculation, look for a file named onsite.dm in the OUT.prefix directory, and make replacements there. The format of the file is rather straight-forward.
- **Default**: 0

### onsite_radius

- **Type**: Real
- **Availability**: *dft_plus_u is set to 1*
- **Description**: - The onsite_radius parameter facilitates modulation of the single-zeta portion of numerical atomic orbitals used for DFT+U projections.
  - The modulation algorithm applies a smooth truncation to the orbital tail followed by normalization. A representative profile is $f(r)=\frac{1}{2}\left[1+\operatorname{erf}\!\left(\frac{r_c-r}{\sigma}\right)\right]$, where $r_c$ is the cutoff radius and $\sigma=\gamma r_c$ controls smoothness.
- **Default**: 3.0
- **Unit**: Bohr

[back to top](#full-list-of-input-keywords)

## Spin-Constrained DFT

### sc_mag_switch

- **Type**: Boolean
- **Description**: Switch to control spin-constrained DFT calculation
- **Default**: False

### decay_grad_switch

- **Type**: Boolean
- **Description**: Switch to control gradient break condition in spin-constrained DFT
- **Default**: False

### sc_thr

- **Type**: Real
- **Availability**: *sc_mag_switch is true*
- **Description**: Convergence criterion of spin-constrained iteration (RMS) in uB
- **Default**: 1.0e-6
- **Unit**: uB

### nsc

- **Type**: Integer
- **Availability**: *sc_mag_switch is true*
- **Description**: Maximal number of spin-constrained iteration
- **Default**: 100

### nsc_min

- **Type**: Integer
- **Availability**: *sc_mag_switch is true*
- **Description**: Minimum number of spin-constrained iteration
- **Default**: 2

### sc_scf_nmin

- **Type**: Integer
- **Availability**: *sc_mag_switch is true*
- **Description**: Minimum number of outer scf loop before initializing lambda loop
- **Default**: 2

### alpha_trial

- **Type**: Real
- **Availability**: *sc_mag_switch is true*
- **Description**: Initial trial step size for lambda in eV/uB^2
- **Default**: 0.01
- **Unit**: eV/uB^2

### sccut

- **Type**: Real
- **Availability**: *sc_mag_switch is true*
- **Description**: Maximal step size for lambda in eV/uB
- **Default**: 3.0
- **Unit**: eV/uB

### sc_drop_thr

- **Type**: Real
- **Availability**: *sc_mag_switch is true*
- **Description**: Convergence criterion ratio of lambda iteration in Spin-constrained DFT
- **Default**: 1.0e-2

### sc_scf_thr

- **Type**: Real
- **Availability**: *sc_mag_switch is true*
- **Description**: Density error threshold for inner loop of spin-constrained SCF
- **Default**: 1.0e-4

[back to top](#full-list-of-input-keywords)

## vdW correction

### vdw_method

- **Type**: String
- **Description**: Specifies the method used for Van der Waals (VdW) correction. Available options are:
  - d2: Grimme's D2 dispersion correction method
  - d3_0: Grimme's DFT-D3(0) dispersion correction method (zero-damping)
  - d3_bj: Grimme's DFTD3(BJ) dispersion correction method (BJ-damping)
  - none: no vdW correction

  > Note: ABACUS supports automatic setting of DFT-D3 parameters for common functionals. To benefit from this feature, please specify the parameter dft_functional explicitly, otherwise the autoset procedure will crash. If not satisfied with the built-in parameters, any manual setting on vdw_s6, vdw_s8, vdw_a1 and vdw_a2 will overwrite the automatic values.
- **Default**: none

### vdw_s6

- **Type**: String
- **Availability**: *vdw_method is set to d2, d3_0, or d3_bj*
- **Description**: This scale factor is used to optimize the interaction energy deviations in van der Waals (vdW) corrected calculations. The recommended values of this parameter are dependent on the chosen vdW correction method and the DFT functional being used. For DFT-D2, the recommended values are 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). If not set, will use values of PBE functional. For DFT-D3, recommended values with different DFT functionals can be found on the here. If not set, will search in ABACUS built-in dataset based on the dft_functional keywords. User set value will overwrite the searched value.

### vdw_s8

- **Type**: String
- **Availability**: *vdw_method is set to d3_0 or d3_bj*
- **Description**: This scale factor is relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the webpage. If not set, will search in ABACUS built-in dataset based on the dft_functional keywords. User set value will overwrite the searched value.

### vdw_a1

- **Type**: String
- **Availability**: *vdw_method is set to d3_0 or d3_bj*
- **Description**: This damping function parameter is relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the webpage. If not set, will search in ABACUS built-in dataset based on the dft_functional keywords. User set value will overwrite the searched value.

### vdw_a2

- **Type**: String
- **Availability**: *vdw_method is set to d3_0 or d3_bj*
- **Description**: This damping function parameter is only relevant for D3(0) and D3(BJ) van der Waals (vdW) correction methods. The recommended values of this parameter with different DFT functionals can be found on the webpage. If not set, will search in ABACUS built-in dataset based on the dft_functional keywords. User set value will overwrite the searched value.

### vdw_d

- **Type**: Real
- **Availability**: *vdw_method is set to d2*
- **Description**: Controls the damping rate of the damping function in the DFT-D2 method.
- **Default**: 20

### vdw_abc

- **Type**: Boolean
- **Availability**: *vdw_method is set to d3_0 or d3_bj*
- **Description**: Determines whether three-body terms are calculated for DFT-D3 methods.
  - True: ABACUS will calculate the three-body term.
  - False: The three-body term is not included.
- **Default**: False

### vdw_c6_file

- **Type**: String
- **Availability**: *vdw_method is set to d2*
- **Description**: Specifies the name of the file containing parameters for each element when using the D2 method. If not set, ABACUS uses the default parameters (Jnm6/mol) stored in the program. To manually set the parameters, provide a file containing the parameters. An example is given by:

  H 0.1 Si 9.0

  Namely, each line contains the element name and the corresponding parameter.
- **Default**: default

### vdw_c6_unit

- **Type**: String
- **Availability**: *vdw_C6_file is not default*
- **Description**: Specifies the unit of the provided parameters in the D2 method. Available options are:
  - Jnm6/mol (J nm^6/mol)
  - eVA (eV Angstrom)
- **Default**: Jnm6/mol

### vdw_r0_file

- **Type**: String
- **Availability**: *vdw_method is set to d2*
- **Description**: Specifies the name of the file containing parameters for each element when using the D2 method. If not set, ABACUS uses the default parameters (Angstrom) stored in the program. To manually set the parameters, provide a file containing the parameters. An example is given by:

  Li 1.0 Cl 2.0

  Namely, each line contains the element name and the corresponding parameter.
- **Default**: default

### vdw_r0_unit

- **Type**: String
- **Availability**: *vdw_R0_file is not default*
- **Description**: Specifies the unit for the parameters in the D2 method when manually set by the user. Available options are:
  - A (Angstrom)
  - Bohr
- **Default**: A

### vdw_cutoff_type

- **Type**: String
- **Description**: Determines the method used for specifying the cutoff radius in periodic systems when applying Van der Waals correction. Available options are:
  - radius: The supercell is selected within a sphere centered at the origin with a radius defined by vdw_cutoff_radius.
  - period: The extent of the supercell is explicitly specified using the vdw_cutoff_period keyword.
- **Default**: radius

### vdw_cutoff_radius

- **Type**: String
- **Availability**: *vdw_cutoff_type is set to radius*
- **Description**: Defines the radius of the cutoff sphere when vdw_cutoff_type is set to radius. The default values depend on the chosen vdw_method.
- **Unit**: defined by vdw_radius_unit (default Bohr)

### vdw_radius_unit

- **Type**: String
- **Availability**: *vdw_cutoff_type is set to radius*
- **Description**: Specify the unit of vdw_cutoff_radius. Available options are:
  - A(Angstrom)
  - Bohr
- **Default**: Bohr

### vdw_cutoff_period

- **Type**: Integer Integer Integer
- **Availability**: *vdw_cutoff_type is set to period*
- **Description**: The three integers supplied here explicitly specify the extent of the supercell in the directions of the three basis lattice vectors.
- **Default**: 3 3 3

### vdw_cn_thr

- **Type**: Real
- **Availability**: *vdw_method is set to d3_0 or d3_bj*
- **Description**: The cutoff radius when calculating coordination numbers.
- **Default**: 40
- **Unit**: defined by vdw_cn_thr_unit (default: Bohr)

### vdw_cn_thr_unit

- **Type**: String
- **Description**: Unit of the coordination number cutoff (vdw_cn_thr). Available options are:
  - A(Angstrom)
  - Bohr
- **Default**: Bohr

[back to top](#full-list-of-input-keywords)

## Berry phase and wannier90 interface

### berry_phase

- **Type**: Boolean
- **Description**: Controls the calculation of Berry phase
  - true: Calculate Berry phase.
  - false: Do not calculate Berry phase.
- **Default**: false

### gdir

- **Type**: Integer
- **Description**: The direction of the polarization in the lattice vector for Berry phase calculation
  - 1: Calculate the polarization in the direction of the lattice vector a_1 defined in the STRU file.
  - 2: Calculate the polarization in the direction of the lattice vector a_2 defined in the STRU file.
  - 3: Calculate the polarization in the direction of the lattice vector a_3 defined in the STRU file.
- **Default**: 3

### towannier90

- **Type**: Boolean
- **Description**: Controls the generation of files for the Wannier90 code.
  - 1: Generate files for the Wannier90 code.
  - 0: Do not generate files for the Wannier90 code.
- **Default**: 0

### nnkpfile

- **Type**: String
- **Description**: The file name generated when running "wannier90 -pp ..." command
- **Default**: seedname.nnkp

### wannier_method

- **Type**: Integer
- **Description**: Only available on LCAO basis, using different methods to generate "\.mmn" file and "\.amn" file.
  - 1: Calculated using the lcao_in_pw method, the calculation accuracy can be improved by increasing ecutwfc to maintain consistency with the pw basis set results.
  - 2: The overlap between atomic orbitals is calculated using grid integration. The radial grid points are generated using the Gauss-Legendre method, while the spherical grid points are generated using the Lebedev-Laikov method.
- **Default**: 1

### wannier_spin

- **Type**: String
- **Description**: The spin direction for the Wannier function calculation when nspin is set to 2
  - up: Calculate spin up for the Wannier function.
  - down: Calculate spin down for the Wannier function.
- **Default**: up

### out_wannier_mmn

- **Type**: Boolean
- **Description**: Write the "*.mmn" file or not.
  - 0: don't write the "*.mmn" file.
  - 1: write the "*.mmn" file.
- **Default**: 1

### out_wannier_amn

- **Type**: Boolean
- **Description**: Write the "*.amn" file or not.
  - 0: don't write the "*.amn" file.
  - 1: write the "*.amn" file.
- **Default**: 1

### out_wannier_eig

- **Type**: Boolean
- **Description**: Write the "*.eig" file or not.
  - 0: don't write the "*.eig" file.
  - 1: write the "*.eig" file.
- **Default**: 1

### out_wannier_unk

- **Type**: Boolean
- **Description**: Write the "UNK.*" file or not.
  - 0: don't write the "UNK.*" file.
  - 1: write the "UNK.*" file.
- **Default**: 0

### out_wannier_wvfn_formatted

- **Type**: Boolean
- **Description**: Write the "UNK.*" file in ASCII format or binary format.
  - 0: write the "UNK.*" file in binary format.
  - 1: write the "UNK.*" file in ASCII format (text file format).
- **Default**: 1

[back to top](#full-list-of-input-keywords)

## RT-TDDFT: Real-Time Time-Dependent Density Functional Theory

### estep_per_md

- **Type**: Integer
- **Description**: The number of electronic propagation steps between two ionic steps.
- **Default**: 1

### td_dt

- **Type**: Real
- **Description**: The time step used in electronic propagation. Setting td_dt will reset the value of md_dt to td_dt * estep_per_md.
- **Default**: md_dt / estep_per_md
- **Unit**: fs

### td_edm

- **Type**: Integer
- **Description**: Method to calculate the energy-density matrix, mainly affects the calculation of force and stress.
  - 0: Using the original formula.
  - 1: Using the formula for ground state (deprecated). Note that this usually does not hold if wave function is not the eigenstate of the Hamiltonian.
- **Default**: 0

### td_print_eij

- **Type**: Real
- **Description**: Controls the printing of Hamiltonian matrix elements.
  - &lt; 0: Suppress all output.
  - &gt;= 0: Print only elements with either i or j exceeding td_print_eij.
- **Default**: -1
- **Unit**: Ry

### td_propagator

- **Type**: Integer
- **Description**: Methods of electronic propagation.
  - 0: Crank-Nicolson, based on matrix inversion.
  - 1: 4th-order Taylor expansion of exponential.
  - 2: Enforced time-reversal symmetry (ETRS).
  - 3: Crank-Nicolson, based on solving linear equation.
- **Default**: 0

### td_vext

- **Type**: Boolean
- **Description**: - True: Add a laser-material interaction (external electric field).
  - False: No external electric field.
- **Default**: False

### td_vext_dire

- **Type**: String
- **Description**: Specifies the direction(s) of the external electric field when td_vext is enabled. For example, td_vext_dire 1 2 indicates that external electric fields are applied to both the x and y directions simultaneously. Electric field parameters can also be written as strings. For example, td_gauss_phase 0 1.5707963 indicates that the Gaussian type electric fields in the x and y directions have a phase delay of pi/2.
  - 1: The external field direction is along the x-axis.
  - 2: The external field direction is along the y-axis.
  - 3: The external field direction is along the z-axis.
- **Default**: 1

### td_stype

- **Type**: Integer
- **Description**: Type of electric field in the space domain, i.e. the gauge of the electric field.
  - 0: Length gauge.
  - 1: Velocity gauge.
  - 2: Hybrid gauge. See J. Chem. Theory Comput. 2025, 21, 3335-3341 for more information.
- **Default**: 0

### td_ttype

- **Type**: String
- **Description**: Type of electric field in the time domain.
  - 0: Gaussian type function.
  - 1: Trapezoid type function.
  - 2: Trigonometric type function.
  - 3: Heaviside type function.
- **Default**: 0

### td_tstart

- **Type**: Integer
- **Description**: The initial time step when the time-dependent electric field is activated.
- **Default**: 1

### td_tend

- **Type**: Integer
- **Description**: The final time step when the time-dependent electric field is deactivated. The field remains active between td_tstart and td_tend.
- **Default**: 1000

### td_lcut1

- **Type**: Real
- **Description**: The lower bound of the interval in the length gauge RT-TDDFT, where the coordinate is the fractional coordinate.
- **Default**: 0.05

### td_lcut2

- **Type**: Real
- **Description**: The upper bound of the interval in the length gauge RT-TDDFT, where the coordinate is the fractional coordinate.
- **Default**: 0.95

### td_gauss_freq

- **Type**: String
- **Description**: Frequency of the Gaussian type electric field.
- **Default**: 22.13
- **Unit**: 1/fs

### td_gauss_phase

- **Type**: String
- **Description**: Phase of the Gaussian type electric field.
- **Default**: 0.0

### td_gauss_sigma

- **Type**: String
- **Description**: Pulse width (standard deviation) of the Gaussian type electric field.
- **Default**: 30.0
- **Unit**: fs

### td_gauss_t0

- **Type**: String
- **Description**: Step number of the time center of the Gaussian type electric field.
- **Default**: 100

### td_gauss_amp

- **Type**: String
- **Description**: Amplitude of the Gaussian type electric field.
- **Default**: 0.25
- **Unit**: V/Angstrom

### td_trape_freq

- **Type**: String
- **Description**: Frequency of the trapezoid type electric field.
- **Default**: 1.60
- **Unit**: 1/fs

### td_trape_phase

- **Type**: String
- **Description**: Phase of the trapezoid type electric field.
- **Default**: 0.0

### td_trape_t1

- **Type**: String
- **Description**: Step number of the time interval t1 of the trapezoid type electric field.
- **Default**: 1875

### td_trape_t2

- **Type**: String
- **Description**: Step number of the time interval t2 of the trapezoid type electric field.
- **Default**: 5625

### td_trape_t3

- **Type**: String
- **Description**: Step number of the time interval t3 of the trapezoid type electric field.
- **Default**: 7500

### td_trape_amp

- **Type**: String
- **Description**: Amplitude of the trapezoid type electric field.
- **Default**: 2.74
- **Unit**: V/Angstrom

### td_trigo_freq1

- **Type**: String
- **Description**: Frequency 1 of the trigonometric type electric field.
- **Default**: 1.164656
- **Unit**: 1/fs

### td_trigo_freq2

- **Type**: String
- **Description**: Frequency 2 of the trigonometric type electric field.
- **Default**: 0.029116
- **Unit**: 1/fs

### td_trigo_phase1

- **Type**: String
- **Description**: Phase 1 of the trigonometric type electric field.
- **Default**: 0.0

### td_trigo_phase2

- **Type**: String
- **Description**: Phase 2 of the trigonometric type electric field.
- **Default**: 0.0

### td_trigo_amp

- **Type**: String
- **Description**: Amplitude of the trigonometric type electric field.
- **Default**: 2.74
- **Unit**: V/Angstrom

### td_heavi_t0

- **Type**: String
- **Description**: Step number of the switch time of the Heaviside type electric field.
- **Default**: 100

### td_heavi_amp

- **Type**: String
- **Description**: Amplitude of the Heaviside type electric field.
- **Default**: 1.0
- **Unit**: V/Angstrom

### init_vecpot_file

- **Type**: Boolean
- **Description**: Initialize vector potential through file or not.
  - True: Initialize vector potential from file At.dat (unit: a.u.). It consists of four columns, representing the step number and vector potential on each direction.
  - False: Calculate vector potential by integrating the electric field.
- **Default**: False

### ocp

- **Type**: Boolean
- **Description**: - True: Fixes the band occupations based on the values specified in ocp_set.
  - False: Does not fix the band occupations.
- **Default**: False

### ocp_set

- **Type**: String
- **Description**: If ocp is set to 1, ocp_set must be provided as a string specifying the occupation numbers for each band across all k-points. The format follows a space-separated pattern, where occupations are assigned sequentially to bands for each k-point. A shorthand notation Nx can be used to repeat a value x for N bands.
  - Example:
  1 10*1 0 1 represents occupations for 13 bands, where the 12th band is fully unoccupied (0), and all others are occupied (1).
  - For a system with multiple k-points, the occupations must be specified for all k-points, following their order in the output file kpoints (may lead to fractional occupations).
  - Incorrect specification of ocp_set could lead to inconsistencies in electron counting, causing the calculation to terminate with an error.
- **Default**: None

### out_dipole

- **Type**: Boolean
- **Description**: - True: Output electric dipole moment.
  - False: Do not output electric dipole moment.
- **Default**: False

### out_current

- **Type**: Integer
- **Description**: - 0: Do not output current.
  - 1: Output current using the two-center integral, faster.
  - 2: Output current using the matrix commutation, more precise.
- **Default**: 0

### out_current_k

- **Type**: Boolean
- **Description**: - True: Output current for each k-points separately.
  - False: Output current in total.
- **Default**: False

### out_efield

- **Type**: Boolean
- **Description**: Whether to output the electric field data to files. When enabled, writes real-time electric field values (unit: V/A) into files named efield_[num].txt, where [num] is the sequential index of the electric field ranges from 0 to N-1 for N configured fields. It is noteworthy that the field type sequence follows td_ttype, while the direction sequence follows td_vext_dire.
  - True: Output electric field.
  - False: Do not output electric field.
- **Default**: False

### out_vecpot

- **Type**: Boolean
- **Description**: Output vector potential or not (unit: a.u.).
  - True: Output vector potential into file At.dat.
  - False: Do not output vector potential.
- **Default**: False

[back to top](#full-list-of-input-keywords)

## Variables useful for debugging

### nurse

- **Type**: Integer
- **Description**: Debugging flag for developers
- **Default**: 0

### t_in_h

- **Type**: Boolean
- **Description**: Specify whether to include kinetic term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vl_in_h

- **Type**: Boolean
- **Description**: Specify whether to include local pseudopotential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vnl_in_h

- **Type**: Boolean
- **Description**: Specify whether to include non-local pseudopotential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vh_in_h

- **Type**: Boolean
- **Description**: Specify whether to include Hartree potential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### vion_in_h

- **Type**: Boolean
- **Description**: Specify whether to include local ionic potential term in obtaining the Hamiltonian matrix.
  - 0: No.
  - 1: Yes.
- **Default**: 1

### test_force

- **Type**: Boolean
- **Description**: Specify whether to output the detailed components in forces.
  - 0: No.
  - 1: Yes.
- **Default**: 0

### test_stress

- **Type**: Boolean
- **Description**: Specify whether to output the detailed components in stress.
  - 0: No.
  - 1: Yes.
- **Default**: 0

### test_skip_ewald

- **Type**: Boolean
- **Description**: Specify whether to skip the calculation of the ewald energy.
  - 0: No.
  - 1: Yes.
- **Default**: 0

[back to top](#full-list-of-input-keywords)

## Electronic conductivities

### cal_cond

- **Type**: Boolean
- **Availability**: *basis_type = pw*
- **Description**: Whether to calculate electronic conductivities.
- **Default**: False

### cond_che_thr

- **Type**: Real
- **Availability**: *esolver_type = sdft*
- **Description**: Control the error of Chebyshev expansions for conductivities.
- **Default**: 1e-8

### cond_dw

- **Type**: Real
- **Availability**: *basis_type = pw*
- **Description**: Frequency interval () for frequency-dependent conductivities.
- **Default**: 0.1
- **Unit**: eV

### cond_wcut

- **Type**: Real
- **Availability**: *basis_type = pw*
- **Description**: Cutoff frequency for frequency-dependent conductivities.
- **Default**: 10.0
- **Unit**: eV

### cond_dt

- **Type**: Real
- **Availability**: *basis_type = pw*
- **Description**: Time interval () to integrate Onsager coefficients.
- **Default**: 0.02
- **Unit**: a.u.

### cond_dtbatch

- **Type**: Integer
- **Availability**: *esolver_type = sdft*
- **Description**: exp(iH\dt\cond_dtbatch) is expanded with Chebyshev expansion to calculate conductivities. It is faster but costs more memory.
  - If cond_dtbatch = 0: Autoset this parameter to make expansion orders larger than 100.
- **Default**: 0

### cond_smear

- **Type**: Integer
- **Description**: Smearing method for conductivities
  - 1: Gaussian smearing
  - 2: Lorentzian smearing
- **Default**: 1

### cond_fwhm

- **Type**: Real
- **Availability**: *basis_type = pw*
- **Description**: FWHM for conductivities. For Gaussian smearing, ; for Lorentzian smearing, .
- **Default**: 0.4
- **Unit**: eV

### cond_nonlocal

- **Type**: Boolean
- **Availability**: *basis_type = pw*
- **Description**: Whether to consider nonlocal potential correction when calculating velocity matrix .
  - True: .
  - False: .
- **Default**: True

[back to top](#full-list-of-input-keywords)

## Implicit solvation model

### imp_sol

- **Type**: Boolean
- **Description**: Calculate implicit solvation correction
- **Default**: False

### eb_k

- **Type**: Real
- **Availability**: *imp_sol is true.*
- **Description**: The relative permittivity of the bulk solvent, 80 for water
- **Default**: 80

### tau

- **Type**: Real
- **Description**: The effective surface tension parameter that describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent which are not captured by the electrostatic terms
- **Default**: 1.0798e-05

### sigma_k

- **Type**: Real
- **Description**: The width of the diffuse cavity that is implicitly determined by the electronic structure of the solute
- **Default**: 0.6

### nc_k

- **Type**: Real
- **Description**: The value of the electron density at which the dielectric cavity forms
- **Default**: 0.00037

[back to top](#full-list-of-input-keywords)

## Quasiatomic Orbital (QO) analysis

### qo_switch

- **Type**: Boolean
- **Description**: Whether to let ABACUS output QO analysis required files
- **Default**: False

### qo_basis

- **Type**: String
- **Description**: Type of QO basis function:
  - hydrogen: hydrogen-like basis
  - pswfc: read basis from pseudopotential
  - szv: single-zeta valence basis
- **Default**: szv

### qo_strategy

- **Type**: Vector of String (1 or n values where n is the number of atomic types)
- **Description**: Strategy to generate radial orbitals for QO analysis. For hydrogen: energy-valence, for pswfc and szv: all
- **Default**: for hydrogen: energy-valence, for pswfc and szv: all

### qo_screening_coeff

- **Type**: Vector of Real (n values where n is the number of atomic types; 1 value allowed for qo_basis=pswfc)
- **Description**: The screening coefficient for each atom type to rescale the shape of radial orbitals
- **Default**: 0.1
- **Unit**: Bohr^-1

### qo_thr

- **Type**: Real
- **Description**: The convergence threshold determining the cutoff of generated orbital. Lower threshold will yield orbital with larger cutoff radius.
- **Default**: 1.0e-6

[back to top](#full-list-of-input-keywords)

## PEXSI

### pexsi_npole

- **Type**: Integer
- **Description**: The number of poles used in the pole expansion method, should be a even number.
- **Default**: 40

### pexsi_inertia

- **Type**: Boolean
- **Description**: Whether inertia counting is used at the very beginning.
- **Default**: True

### pexsi_nmax

- **Type**: Integer
- **Description**: Maximum number of PEXSI iterations after each inertia counting procedure.
- **Default**: 80

### pexsi_comm

- **Type**: Boolean
- **Description**: Whether to construct PSelInv communication pattern.
- **Default**: True

### pexsi_storage

- **Type**: Boolean
- **Description**: Whether to use symmetric storage space used by the Selected Inversion algorithm for symmetric matrices.
- **Default**: True

### pexsi_ordering

- **Type**: Integer
- **Description**: Ordering strategy for factorization and selected inversion. 0: Parallel ordering using ParMETIS, 1: Sequential ordering using METIS, 2: Multiple minimum degree ordering
- **Default**: 0

### pexsi_row_ordering

- **Type**: Integer
- **Description**: Row permutation strategy for factorization and selected inversion, 0: No row permutation, 1: Make the diagonal entry of the matrix larger than the off-diagonal entries.
- **Default**: 1

### pexsi_nproc

- **Type**: Integer
- **Description**: Number of processors for PARMETIS. Only used if pexsi_ordering == 0.
- **Default**: 1

### pexsi_symm

- **Type**: Boolean
- **Description**: Whether the matrix is symmetric.
- **Default**: True

### pexsi_trans

- **Type**: Boolean
- **Description**: Whether to factorize the transpose of the matrix.
- **Default**: False

### pexsi_method

- **Type**: Integer
- **Description**: The pole expansion method to be used. 1 for Cauchy Contour Integral method, 2 for Moussa optimized method.
- **Default**: 1

### pexsi_nproc_pole

- **Type**: Integer
- **Description**: The point parallelizaion of PEXSI. Recommend two points parallelization.
- **Default**: 1

### pexsi_temp

- **Type**: Real
- **Description**: Temperature in Fermi-Dirac distribution, in Ry, should have the same effect as the smearing sigma when smearing method is set to Fermi-Dirac.
- **Default**: 0.015

### pexsi_gap

- **Type**: Real
- **Description**: Spectral gap, this can be set to be 0 in most cases.
- **Default**: 0

### pexsi_delta_e

- **Type**: Real
- **Description**: Upper bound for the spectral radius of S^{-1}H.
- **Default**: 20

### pexsi_mu_lower

- **Type**: Real
- **Description**: Initial guess of lower bound for mu.
- **Default**: -10

### pexsi_mu_upper

- **Type**: Real
- **Description**: Initial guess of upper bound for mu.
- **Default**: 10

### pexsi_mu

- **Type**: Real
- **Description**: Initial guess for mu (for the solver).
- **Default**: 0

### pexsi_mu_thr

- **Type**: Real
- **Description**: Stopping criterion in terms of the chemical potential for the inertia counting procedure.
- **Default**: 0.05

### pexsi_mu_expand

- **Type**: Real
- **Description**: If the chemical potential is not in the initial interval, the interval is expanded by this value.
- **Default**: 0.3

### pexsi_mu_guard

- **Type**: Real
- **Description**: Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.
- **Default**: 0.2

### pexsi_elec_thr

- **Type**: Real
- **Description**: Stopping criterion of the PEXSI iteration in terms of the number of electrons compared to numElectronExact.
- **Default**: 0.001

### pexsi_zero_thr

- **Type**: Real
- **Description**: if the absolute value of CCS matrix element is less than this value, it will be considered as zero.
- **Default**: 1e-10

[back to top](#full-list-of-input-keywords)

## Linear Response TDDFT

### ri_hartree_benchmark

- **Type**: String
- **Description**: Whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)
- **Default**: none

### aims_nbasis

- **Type**: A number(ntype) of Integers
- **Availability**: *ri_hartree_benchmark = aims*
- **Description**: Atomic basis set size for each atom type (with the same order as in STRU) in FHI-aims.
- **Default**: {} (empty list, where ABACUS use its own basis set size)

[back to top](#full-list-of-input-keywords)

## Linear Response TDDFT (Under Development Feature)

### xc_kernel

- **Type**: String
- **Description**: The exchange-correlation kernel used in the calculation. Currently supported: RPA, LDA, PBE, HSE, HF.
- **Default**: LDA

### lr_init_xc_kernel

- **Type**: Vector of String (&gt;=1 values)
- **Description**: The method to initalize the xc kernel.
  - "default": Calculate xc kernel from the ground-state charge density.
  - "file": Read the xc kernel on grid from the provided files. The following words should be the paths of ".cube" files, where the first 1 (nspin==1) or 3 (nspin==2, namely spin-aa, spin-ab and spin-bb) will be read in. The parameter xc_kernel will be invalid. Now only LDA-type kernel is supported as the potential will be calculated by directly multiplying the transition density.
  - "from_charge_file": Calculate fxc from the charge density read from the provided files. The following words should be the paths of ".cube" files, where the first nspin files will be read in.
- **Default**: "default"

### lr_solver

- **Type**: String
- **Description**: The method to solve the Casida equation in LR-TDDFT under Tamm-Dancoff approximation (TDA).
  - dav/dav_subspace/cg: Construct and diagonalize the Hamiltonian matrix iteratively with Davidson/Non-ortho-Davidson/CG algorithm.
  - lapack: Construct the full matrix and directly diagonalize with LAPACK.
  - spectrum: Calculate absorption spectrum only without solving Casida equation.
- **Default**: dav

### lr_thr

- **Type**: Real
- **Description**: The convergence threshold of iterative diagonalization solver for LR-TDDFT. It is a pure-math number with the same meaning as pw_diag_thr, but since the Casida equation is a one-shot eigenvalue problem, it is also the convergence threshold of LR-TDDFT.
- **Default**: 1e-2

### nocc

- **Type**: Integer
- **Description**: The number of occupied orbitals (up to HOMO) used in the LR-TDDFT calculation.
  - Note: If the value is illegal ( &gt; nelec/2 or &lt;= 0), it will be autoset to nelec/2.
- **Default**: nband

### nvirt

- **Type**: Integer
- **Description**: The number of virtual orbitals (starting from LUMO) used in the LR-TDDFT calculation.
- **Default**: 1

### lr_nstates

- **Type**: Integer
- **Description**: The number of 2-particle states to be solved.
- **Default**: 0

### lr_unrestricted

- **Type**: Boolean
- **Description**: Whether to use unrestricted construction for LR-TDDFT (the matrix size will be doubled).
  - True: Always use unrestricted LR-TDDFT.
  - False: Use unrestricted LR-TDDFT only when the system is open-shell.
- **Default**: False

### abs_wavelen_range

- **Type**: Real Real
- **Description**: The range of the wavelength for the absorption spectrum calculation.
- **Default**: 0.0 0.0
- **Unit**: nm

### out_wfc_lr

- **Type**: Boolean
- **Description**: Whether to output the eigenstates (excitation energy) and eigenvectors (excitation amplitude) of the LR-TDDFT calculation. The output files are OUT.{suffix}/Excitation_Amplitude_${processor_rank}.dat.
- **Default**: False

### abs_gauge

- **Type**: String
- **Description**: Whether to use length or velocity gauge to calculate the absorption spectrum in LR-TDDFT.
- **Default**: length

### abs_broadening

- **Type**: Real
- **Description**: The broadening factor for the absorption spectrum calculation.
- **Default**: 0.01

[back to top](#full-list-of-input-keywords)

## Reduced Density Matrix Functional Theory

### rdmft

- **Type**: Boolean
- **Description**: Whether to perform rdmft calculation (reduced density matrix funcional theory)
- **Default**: false

### rdmft_power_alpha

- **Type**: Real
- **Description**: The alpha parameter of power-functional(or other exx-type/hybrid functionals) which used in RDMFT, g(occ_number) = occ_number^alpha
- **Default**: 0.656

[back to top](#full-list-of-input-keywords)

#include "module_io/read_input.h"

#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of input.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Default()
 *     - read empty INPUT file
 *   - Read()
 *     - read input parameters from input files
 *   - Check()
 *     - check_mode = true
 */

#define private public
#include "module_io/input.h"

class InputTest : public testing::Test
{
  protected:
};

TEST_F(InputTest, Default)
{
    ModuleIO::ReadInput readinput(0);
    readinput.check_ntype_flag = false;
    Parameter param;
    readinput.readin_parameters(param, "./support/empty_INPUT");
    const std::string suffix = param.get().suffix;
    EXPECT_EQ(suffix, "ABACUS");
    EXPECT_EQ(param.get().stru_file, "STRU");
    EXPECT_EQ(param.get().kpoint_file, "KPT");
    EXPECT_EQ(param.get().pseudo_dir, "");
    EXPECT_EQ(param.get().orbital_dir, "");
    EXPECT_EQ(param.get().read_file_dir, "auto");
    EXPECT_EQ(param.get().wannier_card, "none");
    EXPECT_EQ(param.get().latname, "none");
    EXPECT_EQ(param.get().calculation, "scf");
    EXPECT_EQ(param.get().esolver_type, "ksdft");
    EXPECT_DOUBLE_EQ(param.get().pseudo_rcut, 15.0);
    EXPECT_FALSE(param.get().pseudo_mesh);
    EXPECT_EQ(param.get().ntype, 0);
    EXPECT_EQ(param.get().nbands, 0);
    EXPECT_EQ(param.get().nbands_sto, 256);
    EXPECT_EQ(param.get().nbands_istate, 5);
    EXPECT_EQ(param.get().pw_seed, 1);
    EXPECT_EQ(param.get().emin_sto, 0.0);
    EXPECT_EQ(param.get().emax_sto, 0.0);
    EXPECT_EQ(param.get().nche_sto, 100);
    EXPECT_EQ(param.get().seed_sto, 0);
    EXPECT_EQ(param.get().initsto_ecut, 0.0);
    EXPECT_EQ(param.get().bndpar, 1);
    EXPECT_EQ(param.get().kpar, 1);
    EXPECT_EQ(param.get().initsto_freq, 0);
    EXPECT_EQ(param.get().method_sto, 2);
    EXPECT_EQ(param.get().npart_sto, 1);
    EXPECT_FALSE(param.get().cal_cond);
    EXPECT_EQ(param.get().dos_nche, 100);
    EXPECT_DOUBLE_EQ(param.get().cond_che_thr, 1e-8);
    EXPECT_DOUBLE_EQ(param.get().cond_dw, 0.1);
    EXPECT_DOUBLE_EQ(param.get().cond_wcut, 10);
    EXPECT_EQ(param.get().cond_dt, 0.02);
    EXPECT_EQ(param.get().cond_dtbatch, 0);
    EXPECT_EQ(param.get().cond_smear, 1);
    EXPECT_DOUBLE_EQ(param.get().cond_fwhm, 0.4);
    EXPECT_TRUE(param.get().cond_nonlocal);
    EXPECT_FALSE(param.get().berry_phase);
    EXPECT_EQ(param.get().gdir, 3);
    EXPECT_FALSE(param.get().towannier90);
    EXPECT_EQ(param.get().nnkpfile, "seedname.nnkp");
    EXPECT_EQ(param.get().wannier_spin, "up");
    EXPECT_EQ(param.get().wannier_method, 1);
    EXPECT_TRUE(param.get().out_wannier_amn);
    EXPECT_TRUE(param.get().out_wannier_mmn);
    EXPECT_FALSE(param.get().out_wannier_unk);
    EXPECT_TRUE(param.get().out_wannier_eig);
    EXPECT_TRUE(param.get().out_wannier_wvfn_formatted);
    EXPECT_DOUBLE_EQ(param.get().kspacing[0], 0.0);
    EXPECT_DOUBLE_EQ(param.get().kspacing[1], 0.0);
    EXPECT_DOUBLE_EQ(param.get().kspacing[2], 0.0);
    EXPECT_DOUBLE_EQ(param.get().min_dist_coef, 0.2);
    EXPECT_EQ(param.get().dft_functional, "default");
    EXPECT_DOUBLE_EQ(param.get().xc_temperature, 0.0);
    EXPECT_EQ(param.get().nspin, 1);
    EXPECT_DOUBLE_EQ(param.get().nelec, 0.0);
    EXPECT_EQ(param.get().lmaxmax, 2);
    EXPECT_EQ(param.get().basis_type, "pw");
    EXPECT_EQ(param.get().ks_solver, "cg");
    EXPECT_DOUBLE_EQ(param.get().search_radius, -1.0);
    EXPECT_TRUE(param.get().search_pbc);
    EXPECT_EQ(param.get().symmetry, "1");
    EXPECT_FALSE(param.get().init_vel);
    EXPECT_DOUBLE_EQ(param.get().ref_cell_factor, 1.0);
    EXPECT_DOUBLE_EQ(param.get().symmetry_prec, 1.0e-6);
    EXPECT_TRUE(param.get().symmetry_autoclose);
    EXPECT_EQ(param.get().cal_force, 0);
    EXPECT_DOUBLE_EQ(param.get().force_thr, 1.0e-3);
    EXPECT_DOUBLE_EQ(param.get().force_thr_ev2, 0);
    EXPECT_DOUBLE_EQ(param.get().stress_thr, 0.5);
    EXPECT_DOUBLE_EQ(param.get().press1, 0.0);
    EXPECT_DOUBLE_EQ(param.get().press2, 0.0);
    EXPECT_DOUBLE_EQ(param.get().press3, 0.0);
    EXPECT_FALSE(param.get().cal_stress);
    EXPECT_EQ(param.get().fixed_axes, "None");
    EXPECT_FALSE(param.get().fixed_ibrav);
    EXPECT_FALSE(param.get().fixed_atoms);
    EXPECT_EQ(param.get().relax_method, "cg");
    EXPECT_DOUBLE_EQ(param.get().relax_cg_thr, 0.5);
    EXPECT_EQ(param.get().out_level, "ie");
    EXPECT_FALSE(param.get().sup.out_md_control);
    EXPECT_TRUE(param.get().relax_new);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w1, 0.01);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w2, 0.5);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmax, 0.8);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmin, 1e-5);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_init, 0.5);
    EXPECT_DOUBLE_EQ(param.get().relax_scale_force, 0.5);
    EXPECT_EQ(param.get().nbspline, -1);
    EXPECT_FALSE(param.get().gamma_only);
    EXPECT_FALSE(param.get().sup.gamma_only_local);
    EXPECT_DOUBLE_EQ(param.get().ecutwfc, 50.0);
    EXPECT_DOUBLE_EQ(param.get().erf_ecut, 0.0);
    EXPECT_DOUBLE_EQ(param.get().erf_height, 0.0);
    EXPECT_DOUBLE_EQ(param.get().erf_sigma, 0.1);
    EXPECT_EQ(param.get().fft_mode, 0);
    EXPECT_EQ(param.get().nx, 0);
    EXPECT_EQ(param.get().ny, 0);
    EXPECT_EQ(param.get().nz, 0);
    EXPECT_EQ(param.get().bx, 1);
    EXPECT_EQ(param.get().by, 1);
    EXPECT_EQ(param.get().bz, 1);
    EXPECT_EQ(param.get().ndx, 0);
    EXPECT_EQ(param.get().ndy, 0);
    EXPECT_EQ(param.get().ndz, 0);
    EXPECT_EQ(param.get().diago_proc, 1);
    EXPECT_EQ(param.get().pw_diag_nmax, 50);
    EXPECT_EQ(param.get().diago_cg_prec, 1);
    EXPECT_EQ(param.get().pw_diag_ndim, 4);
    EXPECT_DOUBLE_EQ(param.get().pw_diag_thr, 1.0e-2);
    EXPECT_EQ(param.get().nb2d, 0);
    EXPECT_EQ(param.get().nurse, 0);
    EXPECT_EQ(param.get().colour, 0);
    EXPECT_EQ(param.get().t_in_h, 1);
    EXPECT_EQ(param.get().vl_in_h, 1);
    EXPECT_EQ(param.get().vnl_in_h, 1);
    EXPECT_EQ(param.get().vh_in_h, 1);
    EXPECT_EQ(param.get().vion_in_h, 1);
    EXPECT_EQ(param.get().test_force, 0);
    EXPECT_EQ(param.get().test_stress, 0);
    EXPECT_DOUBLE_EQ(param.get().scf_thr, 1e-9);
    EXPECT_EQ(param.get().scf_thr_type, 1);
    EXPECT_EQ(param.get().scf_nmax, 100);
    EXPECT_EQ(param.get().relax_nmax, 1);
    EXPECT_EQ(param.get().out_stru, 0);
    EXPECT_EQ(param.get().smearing_method, "gauss");
    EXPECT_DOUBLE_EQ(param.get().smearing_sigma, 0.015);
    EXPECT_EQ(param.get().mixing_mode, "broyden");
    EXPECT_DOUBLE_EQ(param.get().mixing_beta, 0.8);
    EXPECT_EQ(param.get().mixing_ndim, 8);
    EXPECT_DOUBLE_EQ(param.get().mixing_gg0, 1.00);
    EXPECT_EQ(param.get().init_wfc, "atomic");
    EXPECT_EQ(param.get().mem_saver, 0);
    EXPECT_EQ(param.get().printe, 100);
    EXPECT_EQ(param.get().init_chg, "atomic");
    EXPECT_EQ(param.get().chg_extrap, "atomic");
    EXPECT_EQ(param.get().out_freq_elec, 0);
    EXPECT_EQ(param.get().out_freq_ion, 0);
    EXPECT_EQ(param.get().out_chg, 0);
    EXPECT_EQ(param.get().out_dm, 0);
    EXPECT_EQ(param.get().out_dm1, 0);
    EXPECT_EQ(param.get().deepks_out_labels, 0);
    EXPECT_EQ(param.get().deepks_scf, 0);
    EXPECT_EQ(param.get().deepks_equiv, 0);
    EXPECT_EQ(param.get().deepks_bandgap, 0);
    EXPECT_EQ(param.get().deepks_out_unittest, 0);
    EXPECT_EQ(param.get().out_pot, 0);
    EXPECT_EQ(param.get().out_wfc_pw, 0);
    EXPECT_EQ(param.get().out_wfc_r, 0);
    EXPECT_EQ(param.get().out_dos, 0);
    EXPECT_EQ(param.get().out_band[0], 0);
    EXPECT_EQ(param.get().out_band[1], 8);
    EXPECT_EQ(param.get().out_proj_band, 0);
    EXPECT_EQ(param.get().out_mat_hs[0], 0);
    EXPECT_EQ(param.get().out_mat_hs[1], 8);
    EXPECT_EQ(param.get().out_mat_hs2, 0);
    EXPECT_EQ(param.get().out_mat_xc, 0);
    EXPECT_EQ(param.get().out_interval, 1);
    EXPECT_EQ(param.get().out_app_flag, 1);
    EXPECT_EQ(param.get().out_mat_r, 0);
    EXPECT_EQ(param.get().out_wfc_lcao, 0);
    EXPECT_FALSE(param.get().out_alllog);
    EXPECT_DOUBLE_EQ(param.get().dos_emin_ev, -15);
    EXPECT_DOUBLE_EQ(param.get().dos_emax_ev, 15);
    EXPECT_DOUBLE_EQ(param.get().dos_edelta_ev, 0.01);
    EXPECT_DOUBLE_EQ(param.get().dos_scale, 0.01);
    EXPECT_DOUBLE_EQ(param.get().dos_sigma, 0.07);
    EXPECT_FALSE(param.get().out_element_info);
    EXPECT_DOUBLE_EQ(param.get().lcao_ecut, 0);
    EXPECT_DOUBLE_EQ(param.get().lcao_dk, 0.01);
    EXPECT_DOUBLE_EQ(param.get().lcao_dr, 0.01);
    EXPECT_DOUBLE_EQ(param.get().lcao_rmax, 30);
    EXPECT_TRUE(param.get().bessel_nao_smooth);
    EXPECT_DOUBLE_EQ(param.get().bessel_nao_sigma, 0.1);
    EXPECT_EQ(std::stod( param.get().bessel_nao_ecut), 50);
    EXPECT_DOUBLE_EQ(param.get().sup.bessel_nao_rcut, 6.0);
    EXPECT_DOUBLE_EQ(param.get().bessel_nao_tolerence, 1E-12);
    EXPECT_EQ(param.get().bessel_descriptor_lmax, 2);
    EXPECT_TRUE(param.get().bessel_descriptor_smooth);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_sigma, 0.1);
    EXPECT_EQ(std::stod( param.get().bessel_descriptor_ecut), 50);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_rcut, 6.0);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_tolerence, 1E-12);

    EXPECT_FALSE(param.get().efield_flag);
    EXPECT_FALSE(param.get().dip_cor_flag);
    EXPECT_EQ(param.get().efield_dir, 2);
    EXPECT_DOUBLE_EQ(param.get().efield_pos_max, -1.0);
    EXPECT_DOUBLE_EQ(param.get().efield_pos_dec, -1.0);
    EXPECT_DOUBLE_EQ(param.get().efield_amp, 0.0);
    EXPECT_FALSE(param.get().gate_flag);
    EXPECT_DOUBLE_EQ(param.get().zgate, 0.5);
    EXPECT_FALSE(param.get().relax);
    EXPECT_FALSE(param.get().block);
    EXPECT_DOUBLE_EQ(param.get().block_down, 0.45);
    EXPECT_DOUBLE_EQ(param.get().block_up, 0.55);
    EXPECT_DOUBLE_EQ(param.get().block_height, 0.1);
    EXPECT_EQ(param.get().vdw_method, "none");
    EXPECT_EQ(param.get().vdw_s6, "default");
    EXPECT_EQ(param.get().vdw_s8, "default");
    EXPECT_EQ(param.get().vdw_a1, "default");
    EXPECT_EQ(param.get().vdw_a2, "default");
    EXPECT_DOUBLE_EQ(param.get().vdw_d, 20);
    EXPECT_FALSE(param.get().vdw_abc);
    EXPECT_EQ(param.get().vdw_cutoff_radius, "default");
    EXPECT_EQ(param.get().vdw_radius_unit, "Bohr");
    EXPECT_DOUBLE_EQ(param.get().vdw_cn_thr, 40.0);
    EXPECT_EQ(param.get().vdw_cn_thr_unit, "Bohr");
    EXPECT_EQ(param.get().vdw_C6_file, "default");
    EXPECT_EQ(param.get().vdw_C6_unit, "Jnm6/mol");
    EXPECT_EQ(param.get().vdw_R0_file, "default");
    EXPECT_EQ(param.get().vdw_R0_unit, "A");
    EXPECT_EQ(param.get().vdw_cutoff_type, "radius");
    EXPECT_EQ(param.get().vdw_cutoff_period[0], 3);
    EXPECT_EQ(param.get().vdw_cutoff_period[1], 3);
    EXPECT_EQ(param.get().vdw_cutoff_period[2], 3);
    EXPECT_EQ(param.get().exx_hybrid_alpha, "0");
    EXPECT_EQ(param.get().exx_real_number, "0");
    EXPECT_DOUBLE_EQ(param.get().exx_hse_omega, 0.11);
    EXPECT_TRUE(param.get().exx_separate_loop);
    EXPECT_EQ(param.get().exx_hybrid_step, 100);
    EXPECT_DOUBLE_EQ(param.get().exx_lambda, 0.3);
    EXPECT_DOUBLE_EQ(param.get().exx_mixing_beta, 1.0);
    EXPECT_DOUBLE_EQ(param.get().exx_pca_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(param.get().exx_c_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(param.get().exx_v_threshold, 1E-1);
    EXPECT_DOUBLE_EQ(param.get().exx_dm_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(param.get().exx_schwarz_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_threshold, 1E-7);
    EXPECT_DOUBLE_EQ(param.get().exx_c_grad_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(param.get().exx_v_grad_threshold, 1E-1);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_force_threshold, 1E-7);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_stress_threshold, 1E-7);
    EXPECT_EQ(param.get().exx_ccp_rmesh_times, "1");
    EXPECT_DOUBLE_EQ(param.get().rpa_ccp_rmesh_times, 10.0);
    EXPECT_EQ(param.get().exx_distribute_type, "htime");
    EXPECT_EQ(param.get().exx_opt_orb_lmax, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_ecut, 0.0);
    EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_tolerence, 0.0);
    EXPECT_FALSE(param.get().noncolin);
    EXPECT_FALSE(param.get().lspinorb);
    EXPECT_DOUBLE_EQ(param.get().soc_lambda, 1.0);
    EXPECT_DOUBLE_EQ(param.get().td_force_dt, 0.02);
    EXPECT_FALSE(param.get().td_vext);
    EXPECT_EQ(param.get().td_vext_dire[0], 1);
    EXPECT_EQ(param.get().propagator, 0);
    EXPECT_EQ(param.get().td_stype, 0);
    EXPECT_EQ(param.get().td_ttype, "0");
    EXPECT_EQ(param.get().td_tstart, 1);
    EXPECT_EQ(param.get().td_tend, 1000);
    EXPECT_EQ(param.get().td_lcut1, 0.05);
    EXPECT_EQ(param.get().td_lcut2, 0.95);
    EXPECT_EQ(param.get().td_gauss_amp, 0.25);
    EXPECT_EQ(param.get().td_gauss_freq, 22.13);
    EXPECT_EQ(param.get().td_gauss_phase, 0.0);
    EXPECT_EQ(param.get().td_gauss_t0, 100.0);
    EXPECT_EQ(param.get().td_gauss_sigma, 30.0);
    EXPECT_EQ(param.get().td_trape_amp, 2.74);
    EXPECT_EQ(param.get().td_trape_freq, 1.60);
    EXPECT_EQ(param.get().td_trape_phase, 0.0);
    EXPECT_EQ(param.get().td_trape_t1, 1875);
    EXPECT_EQ(param.get().td_trape_t2, 5625);
    EXPECT_EQ(param.get().td_trape_t3, 7500);
    EXPECT_EQ(param.get().td_trigo_freq1, 1.164656);
    EXPECT_EQ(param.get().td_trigo_freq2, 0.029116);
    EXPECT_EQ(param.get().td_trigo_phase1, 0.0);
    EXPECT_EQ(param.get().td_trigo_phase2, 0.0);
    EXPECT_EQ(param.get().td_trigo_amp, 2.74);
    EXPECT_EQ(param.get().td_heavi_t0, 100);
    EXPECT_EQ(param.get().td_heavi_amp, 1.0);
    EXPECT_EQ(param.get().out_dipole, 0);
    EXPECT_EQ(param.get().out_efield, 0);
    EXPECT_EQ(param.get().td_print_eij, -1.0);
    EXPECT_EQ(param.get().td_edm, 0);
    EXPECT_DOUBLE_EQ(param.get().cell_factor, 1.2);
    EXPECT_EQ(param.get().out_mul, 0);
    EXPECT_FALSE(param.get().restart_save);
    EXPECT_FALSE(param.get().restart_load);
    EXPECT_FALSE(param.get().test_skip_ewald);
    EXPECT_EQ(param.get().dft_plus_u, 0);
    EXPECT_FALSE(param.get().yukawa_potential);
    EXPECT_DOUBLE_EQ(param.get().yukawa_lambda, -1.0);
    EXPECT_EQ(param.get().onsite_radius, 0.0);
    EXPECT_EQ(param.get().omc, 0);
    EXPECT_FALSE(param.get().dft_plus_dmft);
    EXPECT_FALSE(param.get().rpa);
    EXPECT_EQ(param.get().imp_sol, 0);
    EXPECT_DOUBLE_EQ(param.get().eb_k, 80.0);
    EXPECT_DOUBLE_EQ(param.get().tau, 1.0798 * 1e-5);
    EXPECT_DOUBLE_EQ(param.get().sigma_k, 0.6);
    EXPECT_DOUBLE_EQ(param.get().nc_k, 0.00037);
    EXPECT_EQ(param.get().of_kinetic, "wt");
    EXPECT_EQ(param.get().of_method, "tn");
    EXPECT_EQ(param.get().of_conv, "energy");
    EXPECT_DOUBLE_EQ(param.get().of_tole, 1e-6);
    EXPECT_DOUBLE_EQ(param.get().of_tolp, 1e-5);
    EXPECT_DOUBLE_EQ(param.get().of_tf_weight, 1.);
    EXPECT_DOUBLE_EQ(param.get().of_vw_weight, 1.);
    EXPECT_DOUBLE_EQ(param.get().of_wt_alpha, 5. / 6.);
    EXPECT_DOUBLE_EQ(param.get().of_wt_beta, 5. / 6.);
    EXPECT_DOUBLE_EQ(param.get().of_wt_rho0, 0.);
    EXPECT_FALSE(param.get().of_hold_rho0);
    EXPECT_DOUBLE_EQ(param.get().of_lkt_a, 1.3);
    EXPECT_TRUE(param.get().of_full_pw);
    EXPECT_EQ(param.get().of_full_pw_dim, 0);
    EXPECT_FALSE(param.get().of_read_kernel);
    EXPECT_EQ(param.get().of_kernel_file, "WTkernel.txt");
    EXPECT_EQ(param.get().device, "cpu");
    EXPECT_DOUBLE_EQ(param.get().ecutrho, 0.0);
    EXPECT_EQ(param.get().sup.ncx, 0);
    EXPECT_EQ(param.get().sup.ncy, 0);
    EXPECT_EQ(param.get().sup.ncz, 0);
    EXPECT_NEAR(param.get().mdp.lj_epsilon, 0.01032, 1e-7);
    EXPECT_NEAR(param.get().mdp.lj_rcut, 8.5, 1e-7);
    EXPECT_NEAR(param.get().mdp.lj_sigma, 3.405, 1e-7);
    EXPECT_EQ(param.get().mdp.md_damp, 1);
    EXPECT_EQ(param.get().mdp.md_dt, 1);
    EXPECT_EQ(param.get().mdp.md_dumpfreq, 1);
    EXPECT_EQ(param.get().mdp.md_nraise, 1);
    EXPECT_EQ(param.get().cal_syns, 0);
    EXPECT_EQ(param.get().dmax, 0.01);
    EXPECT_EQ(param.get().mdp.md_nstep, 10);
    EXPECT_EQ(param.get().mdp.md_pchain, 1);
    EXPECT_EQ(param.get().mdp.md_pcouple, "none");
    EXPECT_DOUBLE_EQ(param.get().mdp.md_pfirst, -1);
    EXPECT_DOUBLE_EQ(param.get().mdp.md_pfreq, 0);
    EXPECT_DOUBLE_EQ(param.get().mdp.md_plast, -1);
    EXPECT_EQ(param.get().mdp.md_pmode, "iso");
    EXPECT_EQ(param.get().mdp.md_restart, 0);
    EXPECT_EQ(param.get().mdp.md_restartfreq, 5);
    EXPECT_EQ(param.get().mdp.md_seed, -1);
    EXPECT_EQ(param.get().mdp.md_prec_level, 0);
    EXPECT_EQ(param.get().mdp.md_tchain, 1);
    EXPECT_DOUBLE_EQ(param.get().mdp.md_tfirst, -1);
    EXPECT_DOUBLE_EQ(param.get().mdp.md_tfreq, 0);
    EXPECT_EQ(param.get().mdp.md_thermostat, "nhc");
    EXPECT_DOUBLE_EQ(param.get().mdp.md_tlast, -1);
    EXPECT_DOUBLE_EQ(param.get().mdp.md_tolerance, 100);
    EXPECT_EQ(param.get().mdp.md_type, "nvt");
    EXPECT_EQ(param.get().mdp.msst_direction, 2);
    EXPECT_DOUBLE_EQ(param.get().mdp.msst_qmass, -1);
    EXPECT_DOUBLE_EQ(param.get().mdp.msst_tscale, 0.01);
    EXPECT_DOUBLE_EQ(param.get().mdp.msst_vel, 0);
    EXPECT_DOUBLE_EQ(param.get().mdp.msst_vis, 0);
    EXPECT_EQ(param.get().mdp.pot_file, "graph.pb");
    EXPECT_TRUE(param.get().mdp.dump_force);
    EXPECT_TRUE(param.get().mdp.dump_vel);
    EXPECT_TRUE(param.get().mdp.dump_virial);
    EXPECT_EQ(param.get().sc_mag_switch, 0);
    EXPECT_FALSE(param.get().decay_grad_switch);
    EXPECT_DOUBLE_EQ(param.get().sc_thr, 1e-6);
    EXPECT_EQ(param.get().nsc, 100);
    EXPECT_EQ(param.get().nsc_min, 2);
    EXPECT_EQ(param.get().sc_scf_nmin, 2);
    EXPECT_DOUBLE_EQ(param.get().alpha_trial, 0.01);
    EXPECT_DOUBLE_EQ(param.get().sccut, 3.0);
    EXPECT_EQ(param.get().sc_file, "none");
}

// TEST_F(InputTest, Read)
// {
// 	std::string input_file = "./support/param.get()";
// 	param.get().Read(input_file);
// 	EXPECT_EQ(param.get().suffix,"autotest");
// 	EXPECT_EQ(param.get().stru_file,"./support/STRU");
// 	EXPECT_EQ(param.get().kpoint_file,"KPT");
// 	EXPECT_EQ(param.get().pseudo_dir,"../../PP_ORB/");
// 	EXPECT_EQ(param.get().orbital_dir,"../../PP_ORB/");
// 	EXPECT_EQ(param.get().read_file_dir,"auto");
// 	EXPECT_EQ(param.get().wannier_card,"none");
// 	EXPECT_EQ(param.get().latname,"none");
// 	EXPECT_EQ(param.get().calculation,"scf");
// 	EXPECT_EQ(param.get().esolver_type,"ksdft");
// 	EXPECT_DOUBLE_EQ(param.get().pseudo_rcut,15.0);
// 	EXPECT_FALSE(param.get().pseudo_mesh);
// 	EXPECT_EQ(param.get().ntype,1);
// 	EXPECT_EQ(param.get().nbands,8);
// 	EXPECT_EQ(param.get().nbands_sto,256);
// 	EXPECT_EQ(param.get().nbands_istate,5);
// 	EXPECT_EQ(param.get().pw_seed,1);
// 	EXPECT_EQ(param.get().emin_sto,0.0);
// 	EXPECT_EQ(param.get().emax_sto,0.0);
// 	EXPECT_EQ(param.get().nche_sto,100);
//         EXPECT_EQ(param.get().seed_sto,0);
// 		EXPECT_EQ(param.get().initsto_ecut,0.0);
//         EXPECT_EQ(param.get().bndpar,1);
//         EXPECT_EQ(param.get().kpar,1);
//         EXPECT_EQ(param.get().initsto_freq,0);
//         EXPECT_EQ(param.get().method_sto,3);
//         EXPECT_EQ(param.get().npart_sto,1);
//         EXPECT_FALSE(param.get().cal_cond);
//         EXPECT_EQ(param.get().dos_nche,100);
//         EXPECT_DOUBLE_EQ(param.get().cond_che_thr,1e-8);
//         EXPECT_DOUBLE_EQ(param.get().cond_dw,0.1);
//         EXPECT_DOUBLE_EQ(param.get().cond_wcut,10);
//         EXPECT_EQ(param.get().cond_dt,0.07);
// 		EXPECT_EQ(param.get().cond_dtbatch,2);
//         EXPECT_DOUBLE_EQ(param.get().cond_fwhm,0.3);
//         EXPECT_TRUE(param.get().cond_nonlocal);
//         EXPECT_FALSE(param.get().berry_phase);
//         EXPECT_EQ(param.get().gdir,3);
//         EXPECT_FALSE(param.get().towannier90);
//         EXPECT_EQ(param.get().nnkpfile,"seedname.nnkp");
//         EXPECT_EQ(param.get().wannier_spin,"up");
//         EXPECT_EQ(param.get().wannier_method,1);
// 		EXPECT_TRUE(param.get().out_wannier_amn);
// 		EXPECT_TRUE(param.get().out_wannier_mmn);
// 		EXPECT_TRUE(param.get().out_wannier_unk);
// 		EXPECT_TRUE(param.get().out_wannier_eig);
//         EXPECT_TRUE(param.get().out_wannier_wvfn_formatted);
//         EXPECT_DOUBLE_EQ(param.get().kspacing[0], 0.0);
//         EXPECT_DOUBLE_EQ(param.get().kspacing[1],0.0);
//         EXPECT_DOUBLE_EQ(param.get().kspacing[2],0.0);
//         EXPECT_DOUBLE_EQ(param.get().min_dist_coef,0.2);
//         EXPECT_EQ(param.get().dft_functional,"hse");
//         EXPECT_DOUBLE_EQ(param.get().xc_temperature,0.0);
//         EXPECT_EQ(param.get().nspin,1);
//         EXPECT_DOUBLE_EQ(param.get().nelec,0.0);
//         EXPECT_EQ(param.get().lmaxmax,2);
//         EXPECT_EQ(param.get().basis_type,"lcao");
//         EXPECT_EQ(param.get().ks_solver,"genelpa");
//         EXPECT_DOUBLE_EQ(param.get().search_radius,-1.0);
//         EXPECT_TRUE(param.get().search_pbc);
//         EXPECT_EQ(param.get().symmetry,"1");
//         EXPECT_FALSE(param.get().init_vel);
//         EXPECT_DOUBLE_EQ(param.get().symmetry_prec, 1.0e-6);
//         EXPECT_TRUE(param.get().symmetry_autoclose);
//         EXPECT_EQ(param.get().cal_force, 0);
//         EXPECT_NEAR(param.get().force_thr,1.0e-3,1.0e-7);
//         EXPECT_DOUBLE_EQ(param.get().force_thr_ev2,0);
//         EXPECT_DOUBLE_EQ(param.get().stress_thr,1.0e-2);
//         EXPECT_DOUBLE_EQ(param.get().press1,0.0);
//         EXPECT_DOUBLE_EQ(param.get().press2,0.0);
//         EXPECT_DOUBLE_EQ(param.get().press3,0.0);
//         EXPECT_FALSE(param.get().cal_stress);
//         EXPECT_EQ(param.get().fixed_axes,"None");
//         EXPECT_FALSE(param.get().fixed_ibrav);
//         EXPECT_FALSE(param.get().fixed_atoms);
//         EXPECT_EQ(param.get().relax_method,"cg");
//         EXPECT_DOUBLE_EQ(param.get().relax_cg_thr,0.5);
//         EXPECT_EQ(param.get().out_level,"ie");
//         EXPECT_TRUE(param.get().out_md_control);
//         EXPECT_TRUE(param.get().relax_new);
//         EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w1,0.01);
//         EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w2,0.5);
//         EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmax,0.8);
//         EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmin,1e-5);
//         EXPECT_DOUBLE_EQ(param.get().relax_bfgs_init,0.5);
//         EXPECT_DOUBLE_EQ(param.get().relax_scale_force,0.5);
//         EXPECT_EQ(param.get().nbspline,-1);
//         EXPECT_TRUE(param.get().gamma_only);
//         EXPECT_TRUE(param.get().gamma_only_local);
//         EXPECT_DOUBLE_EQ(param.get().ecutwfc,20.0);
//         EXPECT_DOUBLE_EQ(param.get().erf_ecut, 20.0);
//         EXPECT_DOUBLE_EQ(param.get().erf_height, 20.0);
//         EXPECT_DOUBLE_EQ(param.get().erf_sigma, 4.0);
//         EXPECT_DOUBLE_EQ(param.get().ecutrho, 0.0);
// 		EXPECT_EQ(param.get().fft_mode,0);
//         EXPECT_EQ(param.get().ncx,0);
//         EXPECT_EQ(param.get().ncy,0);
//         EXPECT_EQ(param.get().ncz,0);
//         EXPECT_EQ(param.get().nx,0);
//         EXPECT_EQ(param.get().ny,0);
//         EXPECT_EQ(param.get().nz,0);
//         EXPECT_EQ(param.get().bx,2);
//         EXPECT_EQ(param.get().by,2);
//         EXPECT_EQ(param.get().bz,2);
//         EXPECT_EQ(param.get().ndx, 0);
//         EXPECT_EQ(param.get().ndy, 0);
//         EXPECT_EQ(param.get().ndz, 0);
//         EXPECT_EQ(param.get().diago_proc,4);
//         EXPECT_EQ(param.get().pw_diag_nmax,50);
//         EXPECT_EQ(param.get().diago_cg_prec,1);
//         EXPECT_EQ(param.get().pw_diag_ndim,4);
//         EXPECT_DOUBLE_EQ(param.get().pw_diag_thr,1.0e-2);
//         EXPECT_EQ(param.get().nb2d,0);
//         EXPECT_EQ(param.get().nurse,0);
//         EXPECT_EQ(param.get().colour,0);
//         EXPECT_EQ(param.get().t_in_h,1);
//         EXPECT_EQ(param.get().vl_in_h,1);
//         EXPECT_EQ(param.get().vnl_in_h,1);
//         EXPECT_EQ(param.get().vh_in_h,1);
//         EXPECT_EQ(param.get().vion_in_h,1);
//         EXPECT_EQ(param.get().test_force,0);
//         EXPECT_EQ(param.get().test_stress,0);
//         EXPECT_NEAR(param.get().scf_thr,1.0e-8,1.0e-15);
//         EXPECT_EQ(param.get().scf_nmax,50);
//         EXPECT_EQ(param.get().relax_nmax,1);
//         EXPECT_EQ(param.get().out_stru,0);
//         EXPECT_EQ(param.get().occupations,"smearing");
//         EXPECT_EQ(param.get().smearing_method,"gauss");
//         EXPECT_DOUBLE_EQ(param.get().smearing_sigma,0.002);
//         EXPECT_EQ(param.get().mixing_mode,"broyden");
//         EXPECT_DOUBLE_EQ(param.get().mixing_beta,0.7);
//         EXPECT_EQ(param.get().mixing_ndim,8);
//         EXPECT_DOUBLE_EQ(param.get().mixing_gg0,0.00);
//         EXPECT_EQ(param.get().init_wfc,"atomic");
//         EXPECT_EQ(param.get().mem_saver,0);
//         EXPECT_EQ(param.get().printe,100);
//         EXPECT_EQ(param.get().init_chg,"atomic");
//         EXPECT_EQ(param.get().chg_extrap,"atomic");
//         EXPECT_EQ(param.get().out_freq_elec,0);
//         EXPECT_EQ(param.get().out_freq_ion,0);
//         EXPECT_EQ(param.get().out_chg,0);
//         EXPECT_EQ(param.get().out_dm,0);
//         EXPECT_EQ(param.get().out_dm1,0);
//         EXPECT_EQ(param.get().deepks_out_labels,0);
//         EXPECT_EQ(param.get().deepks_scf,0);
// 		EXPECT_EQ(param.get().deepks_equiv,0);
//         EXPECT_EQ(param.get().deepks_bandgap,0);
//         EXPECT_EQ(param.get().deepks_out_unittest,0);
//         EXPECT_EQ(param.get().out_pot,2);
//         EXPECT_EQ(param.get().out_wfc_pw,0);
//         EXPECT_EQ(param.get().out_wfc_r,0);
//         EXPECT_EQ(param.get().out_dos,0);
//         EXPECT_EQ(param.get().out_band[0],0);
// 		EXPECT_EQ(param.get().out_band[1],8);
//         EXPECT_EQ(param.get().out_proj_band,0);
//         EXPECT_EQ(param.get().out_mat_hs[0],0);
// 		EXPECT_EQ(param.get().out_mat_hs[1],8);
//         EXPECT_EQ(param.get().out_mat_hs2,0);
//         EXPECT_EQ(param.get().out_mat_xc, 0);
//         EXPECT_EQ(param.get().out_interval,1);
//         EXPECT_EQ(param.get().out_app_flag,0);
//         EXPECT_EQ(param.get().out_mat_r,0);
//         EXPECT_FALSE(param.get().out_wfc_lcao);
//         EXPECT_FALSE(param.get().out_alllog);
//         EXPECT_DOUBLE_EQ(param.get().dos_emin_ev,-15);
//         EXPECT_DOUBLE_EQ(param.get().dos_emax_ev,15);
//         EXPECT_DOUBLE_EQ(param.get().dos_edelta_ev,0.01);
//         EXPECT_DOUBLE_EQ(param.get().dos_scale,0.01);
//         EXPECT_DOUBLE_EQ(param.get().dos_sigma,0.07);
//         EXPECT_FALSE(param.get().out_element_info);
//         EXPECT_DOUBLE_EQ(param.get().lcao_ecut,20);
//         EXPECT_DOUBLE_EQ(param.get().lcao_dk,0.01);
//         EXPECT_DOUBLE_EQ(param.get().lcao_dr,0.01);
//         EXPECT_DOUBLE_EQ(param.get().lcao_rmax,30);
// 		EXPECT_TRUE(param.get().bessel_nao_smooth);
// 		EXPECT_DOUBLE_EQ(param.get().bessel_nao_sigma, 0.1);
// 		EXPECT_EQ(param.get().bessel_nao_ecut, "default");
// 		EXPECT_DOUBLE_EQ(param.get().bessel_nao_rcut, 6.0);
// 		EXPECT_DOUBLE_EQ(param.get().bessel_nao_tolerence, 1E-12);
// 		EXPECT_EQ(param.get().bessel_descriptor_lmax, 2);
// 		EXPECT_TRUE(param.get().bessel_descriptor_smooth);
// 		EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_sigma, 0.1);
// 		EXPECT_EQ(param.get().bessel_descriptor_ecut, "default");
// 		EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_rcut, 6.0);
// 		EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_tolerence, 1E-12);
//         EXPECT_FALSE(param.get().efield_flag);
//         EXPECT_FALSE(param.get().dip_cor_flag);
//         EXPECT_EQ(param.get().efield_dir,2);
//         EXPECT_DOUBLE_EQ(param.get().efield_pos_max,0.5);
//         EXPECT_DOUBLE_EQ(param.get().efield_pos_dec,0.1);
//         EXPECT_DOUBLE_EQ(param.get().efield_amp ,0.0);
//         EXPECT_FALSE(param.get().gate_flag);
//         EXPECT_DOUBLE_EQ(param.get().zgate,0.5);
//         EXPECT_FALSE(param.get().relax);
//         EXPECT_FALSE(param.get().block);
//         EXPECT_DOUBLE_EQ(param.get().block_down,0.45);
//         EXPECT_DOUBLE_EQ(param.get().block_up,0.55);
//         EXPECT_DOUBLE_EQ(param.get().block_height,0.1);
//         EXPECT_EQ(param.get().vdw_method,"d2");
//         EXPECT_EQ(param.get().vdw_s6,"default");
//         EXPECT_EQ(param.get().vdw_s8,"default");
//         EXPECT_EQ(param.get().vdw_a1,"default");
//         EXPECT_EQ(param.get().vdw_a2,"default");
//         EXPECT_DOUBLE_EQ(param.get().vdw_d,20);
//         EXPECT_FALSE(param.get().vdw_abc);
//         EXPECT_EQ(param.get().vdw_cutoff_radius,"default");
//         EXPECT_EQ(param.get().vdw_radius_unit,"Bohr");
//         EXPECT_DOUBLE_EQ(param.get().vdw_cn_thr,40.0);
//         EXPECT_EQ(param.get().vdw_cn_thr_unit,"Bohr");
//         EXPECT_EQ(param.get().vdw_C6_file,"default");
//         EXPECT_EQ(param.get().vdw_C6_unit,"Jnm6/mol");
//         EXPECT_EQ(param.get().vdw_R0_file,"default");
//         EXPECT_EQ(param.get().vdw_R0_unit,"A");
//         EXPECT_EQ(param.get().vdw_cutoff_type,"radius");
//         EXPECT_EQ(param.get().vdw_cutoff_period[0],3);
//         EXPECT_EQ(param.get().vdw_cutoff_period[1],3);
//         EXPECT_EQ(param.get().vdw_cutoff_period[2],3);
//         EXPECT_EQ(param.get().exx_hybrid_alpha,"default");
//         EXPECT_EQ(param.get().exx_real_number,"default");
//         EXPECT_DOUBLE_EQ(param.get().exx_hse_omega,0.11);
//         EXPECT_TRUE(param.get().exx_separate_loop);
//         EXPECT_EQ(param.get().exx_hybrid_step,100);
//         EXPECT_DOUBLE_EQ(param.get().exx_lambda,0.3);
// 		EXPECT_DOUBLE_EQ(param.get().exx_mixing_beta,1.0);
//         EXPECT_DOUBLE_EQ(param.get().exx_pca_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_c_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_v_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_dm_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_schwarz_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_cauchy_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_c_grad_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_v_grad_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_cauchy_force_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_cauchy_stress_threshold,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_ccp_threshold,1E-8);
//         EXPECT_EQ(param.get().exx_ccp_rmesh_times, "default");
//         EXPECT_DOUBLE_EQ(param.get().rpa_ccp_rmesh_times, 10.0);
//         EXPECT_EQ(param.get().exx_distribute_type, "htime");
//         EXPECT_EQ(param.get().exx_opt_orb_lmax,0);
//         EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_ecut,0.0);
//         EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_tolerence,0.0);
//         EXPECT_FALSE(param.get().noncolin);
//         EXPECT_FALSE(param.get().lspinorb);
//         EXPECT_DOUBLE_EQ(param.get().soc_lambda,1.0);
//         EXPECT_EQ(param.get().input_error,0);
//         EXPECT_DOUBLE_EQ(param.get().td_force_dt,0.02);
//         EXPECT_EQ(param.get().td_vext,0);
//         // EXPECT_EQ(param.get().td_vext_dire,"1");
// 		EXPECT_EQ(param.get().propagator,0);
// 		EXPECT_EQ(param.get().td_stype,0);
// 		// EXPECT_EQ(param.get().td_ttype,"0");
// 		EXPECT_EQ(param.get().td_tstart,1);
// 		EXPECT_EQ(param.get().td_tend,1000);
// 		EXPECT_EQ(param.get().td_lcut1,0.05);
// 		EXPECT_EQ(param.get().td_lcut2,0.95);
// 		// EXPECT_EQ(param.get().td_gauss_amp,"0.25");
// 		// EXPECT_EQ(param.get().td_gauss_freq,"22.13");
// 		// EXPECT_EQ(param.get().td_gauss_phase,"0.0");
// 		// EXPECT_EQ(param.get().td_gauss_t0,"100.0");
// 		// EXPECT_EQ(param.get().td_gauss_sigma,"30.0");
// 		// EXPECT_EQ(param.get().td_trape_amp,"2.74");
// 		// EXPECT_EQ(param.get().td_trape_freq,"1.60");
// 		// EXPECT_EQ(param.get().td_trape_phase,"0.0");
// 		// EXPECT_EQ(param.get().td_trape_t1,"1875");
// 		// EXPECT_EQ(param.get().td_trape_t2,"5625");
// 		// EXPECT_EQ(param.get().td_trape_t3,"7500");
// 		// EXPECT_EQ(param.get().td_trigo_freq1,"1.164656");
// 		// EXPECT_EQ(param.get().td_trigo_freq2,"0.029116");
// 		// EXPECT_EQ(param.get().td_trigo_phase1,"0.0");
// 		// EXPECT_EQ(param.get().td_trigo_phase2,"0.0");
// 		// EXPECT_EQ(param.get().td_trigo_amp,"2.74");
// 		// EXPECT_EQ(param.get().td_heavi_t0,"100");
// 		// EXPECT_EQ(param.get().td_heavi_amp,"1.0");
//         EXPECT_EQ(param.get().out_dipole,0);
// 		EXPECT_EQ(param.get().out_efield,0);
// 		EXPECT_EQ(param.get().td_print_eij,-1.0);
// 		EXPECT_EQ(param.get().td_edm,0);
//         EXPECT_DOUBLE_EQ(param.get().cell_factor,1.2);
//         EXPECT_EQ(param.get().out_mul,0);
//         EXPECT_FALSE(param.get().restart_save);
//         EXPECT_FALSE(param.get().restart_load);
//         EXPECT_FALSE(param.get().test_skip_ewald);
//         EXPECT_EQ(param.get().dft_plus_u, 0);
//         EXPECT_FALSE(param.get().yukawa_potential);
//         EXPECT_DOUBLE_EQ(param.get().yukawa_lambda,-1.0);
// 		EXPECT_EQ(param.get().onsite_radius, 0.0);
//         EXPECT_EQ(param.get().omc,0);
//         EXPECT_FALSE(param.get().dft_plus_dmft);
//         EXPECT_FALSE(param.get().rpa);
//         EXPECT_EQ(param.get().imp_sol,0);
//         EXPECT_DOUBLE_EQ(param.get().eb_k,80.0);
//         EXPECT_DOUBLE_EQ(param.get().tau,1.0798 * 1e-5);
//         EXPECT_DOUBLE_EQ(param.get().sigma_k,0.6);
//         EXPECT_DOUBLE_EQ(param.get().nc_k,0.00037);
//         EXPECT_EQ(param.get().of_kinetic,"vw");
//         EXPECT_EQ(param.get().of_method,"tn");
//         EXPECT_EQ(param.get().of_conv,"energy");
//         EXPECT_DOUBLE_EQ(param.get().of_tole,1e-6);
//         EXPECT_DOUBLE_EQ(param.get().of_tolp,1e-5);
//         EXPECT_DOUBLE_EQ(param.get().of_tf_weight,1.);
//         EXPECT_DOUBLE_EQ(param.get().of_vw_weight,1.);
//         EXPECT_DOUBLE_EQ(param.get().of_wt_alpha,0.833333);
//         EXPECT_DOUBLE_EQ(param.get().of_wt_beta,0.833333);
//         EXPECT_DOUBLE_EQ(param.get().of_wt_rho0,1.);
//         EXPECT_FALSE(param.get().of_hold_rho0);
//         EXPECT_DOUBLE_EQ(param.get().of_lkt_a, 1.3);
//         EXPECT_FALSE(param.get().of_full_pw);
//         EXPECT_EQ(param.get().of_full_pw_dim,0);
//         EXPECT_FALSE(param.get().of_read_kernel);
//         EXPECT_EQ(param.get().of_kernel_file,"WTkernel.txt");
//         EXPECT_EQ(param.get().device, "cpu");
//         EXPECT_EQ(param.get().ncx,0);
//         EXPECT_EQ(param.get().ncy,0);
//         EXPECT_EQ(param.get().ncz,0);
//         //EXPECT_NEAR(param.get().force_thr_ev,0.0257112,1e-8);
//         EXPECT_DOUBLE_EQ(param.get().hubbard_u[0],0);
// 	EXPECT_EQ(param.get().orbital_corr[0],-1);
// 	EXPECT_NEAR(param.get().mdp.lj_epsilon,0.01032,1e-7);
// 	EXPECT_NEAR(param.get().mdp.lj_rcut,8.5,1e-7);
// 	EXPECT_NEAR(param.get().mdp.lj_sigma,3.405,1e-7);
// 	EXPECT_EQ(param.get().mdp.md_damp,1);
// 	EXPECT_EQ(param.get().mdp.md_dt,1);
// 	EXPECT_EQ(param.get().mdp.md_dumpfreq,1);
// 	EXPECT_EQ(param.get().mdp.md_nraise,1);
// 	EXPECT_EQ(param.get().cal_syns,0);
// 	EXPECT_EQ(param.get().dmax,0.01);
// 	EXPECT_EQ(param.get().mdp.md_nstep,10);
// 	EXPECT_EQ(param.get().mdp.md_pchain,1);
// 	EXPECT_EQ(param.get().mdp.md_pcouple,"none");
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_pfirst,-1);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_pfreq,0);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_plast,-1);
// 	EXPECT_EQ(param.get().mdp.md_pmode,"iso");
// 	EXPECT_EQ(param.get().mdp.md_restart,0);
// 	EXPECT_EQ(param.get().mdp.md_restartfreq,5);
// 	EXPECT_EQ(param.get().mdp.md_seed,-1);
//     EXPECT_EQ(param.get().mdp.md_prec_level, 2);
//     EXPECT_DOUBLE_EQ(param.get().ref_cell_factor,1.2);
// 	EXPECT_EQ(param.get().mdp.md_tchain,1);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_tfirst,-1);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_tfreq,0);
// 	EXPECT_EQ(param.get().mdp.md_thermostat,"nhc");
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_tlast,-1);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_tolerance,100);
// 	EXPECT_EQ(param.get().mdp.md_type,"nvt");
// 	EXPECT_EQ(param.get().mdp.msst_direction,2);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.msst_qmass,-1);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.msst_tscale,0.01);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.msst_vel,0);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.msst_vis,0);
// 	EXPECT_EQ(param.get().mdp.pot_file,"graph.pb");
//     EXPECT_FALSE(param.get().mdp.dump_force);
//     EXPECT_FALSE(param.get().mdp.dump_vel);
//     EXPECT_FALSE(param.get().mdp.dump_virial);
//     EXPECT_EQ(param.get().sc_mag_switch, 0);
//     EXPECT_TRUE(param.get().decay_grad_switch);
//     EXPECT_DOUBLE_EQ(param.get().sc_thr, 1e-4);
//     EXPECT_EQ(param.get().nsc, 50);
// 	EXPECT_EQ(param.get().nsc_min, 4);
// 	EXPECT_EQ(param.get().sc_scf_nmin, 4);
//     EXPECT_DOUBLE_EQ(param.get().alpha_trial, 0.02);
// 	EXPECT_DOUBLE_EQ(param.get().sccut, 4.0);
//     EXPECT_EQ(param.get().sc_file, "sc.json");
// }

// TEST_F(InputTest, Default_2)
// {
// 	//==================================================
// 	// prepare default parameters for the 1st calling
// 	EXPECT_EQ(param.get().vdw_method,"d2");
// 	EXPECT_EQ(param.get().vdw_s6,"default");
//         EXPECT_EQ(param.get().vdw_s8,"default");
//         EXPECT_EQ(param.get().vdw_a1,"default");
//         EXPECT_EQ(param.get().vdw_a2,"default");
//         EXPECT_EQ(param.get().vdw_cutoff_radius,"default");
// 	EXPECT_NE(param.get().esolver_type,"sdft");
// 	EXPECT_NE(param.get().method_sto,1);
// 	EXPECT_NE(param.get().method_sto,2);
// 	EXPECT_NE(param.get().of_wt_rho0,0.0);
// 	EXPECT_EQ(param.get().exx_hybrid_alpha,"default");
// 	EXPECT_EQ(param.get().dft_functional,"hse");
// 	EXPECT_EQ(param.get().exx_real_number,"default");
// 	EXPECT_TRUE(param.get().gamma_only);
//         EXPECT_EQ(param.get().exx_ccp_rmesh_times,"default");
// 	EXPECT_EQ(param.get().diago_proc,4);
// 	EXPECT_EQ(GlobalV::NPROC,1);
// 	EXPECT_EQ(param.get().calculation,"scf");
// 	EXPECT_EQ(param.get().basis_type,"lcao");
// 	param.get().ks_solver = "default";
// 	param.get().lcao_ecut = 0;
// 	param.get().scf_thr = -1.0;
// 	param.get().scf_thr_type = -1;
//     EXPECT_DOUBLE_EQ(param.get().ecutwfc, 20.0);
//     EXPECT_DOUBLE_EQ(param.get().erf_ecut, 20.0);
//     EXPECT_DOUBLE_EQ(param.get().erf_height, 20.0);
//     EXPECT_DOUBLE_EQ(param.get().erf_sigma, 4.0);
//     param.get().nbndsto_str = "all";
//     param.get().nx = param.get().ny = param.get().nz = 4;
//     param.get().ndx = param.get().ndy = param.get().ndz = 0;
//     // the 1st calling
//     param.get().Default_2();
//     // ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().ndx, 4);
//     EXPECT_EQ(param.get().ndy, 4);
//     EXPECT_EQ(param.get().ndz, 4);
//     EXPECT_FALSE(GlobalV::double_grid);
//     EXPECT_DOUBLE_EQ(param.get().ecutrho, 80.0);
//     EXPECT_EQ(param.get().vdw_s6, "0.75");
//     EXPECT_EQ(param.get().vdw_cutoff_radius, "56.6918");
//     EXPECT_EQ(param.get().bndpar,1);
// 	EXPECT_EQ(param.get().method_sto,2);
// 	EXPECT_TRUE(param.get().of_hold_rho0);
//         EXPECT_EQ(param.get().of_full_pw_dim,0);
// 	EXPECT_EQ(param.get().exx_hybrid_alpha,"0.25");
// 	EXPECT_EQ(param.get().exx_real_number,"1");
//         EXPECT_EQ(param.get().exx_ccp_rmesh_times,"1.5");
// 	EXPECT_EQ(param.get().diago_proc,1);
// 	EXPECT_EQ(param.get().mem_saver,0);
// 	EXPECT_EQ(param.get().relax_nmax,1);
// 	EXPECT_DOUBLE_EQ(param.get().scf_thr,1.0e-7);
// 	EXPECT_EQ(param.get().scf_thr_type,2);
// #ifdef __ELPA
// 	EXPECT_EQ(param.get().ks_solver,"genelpa");
// #else
// 	EXPECT_EQ(param.get().ks_solver,"scalapack_gvx");
// #endif
//         EXPECT_DOUBLE_EQ(param.get().lcao_ecut,20.0);
// 	EXPECT_EQ(param.get().nbands_sto, 0);
// 	//==================================================
// 	// prepare default parameters for the 2nd calling
// 	param.get().vdw_method = "d3_0";
// 	param.get().vdw_s6 = "default";
// 	param.get().vdw_s8 = "default";
// 	param.get().vdw_a1 = "default";
// 	param.get().vdw_a2 = "default";
// 	param.get().vdw_cutoff_radius = "default";
// 	param.get().exx_hybrid_alpha = "default";
// 	param.get().dft_functional = "hf";
// 	param.get().exx_real_number = "default";
// 	param.get().gamma_only = 0;
//         param.get().exx_ccp_rmesh_times = "default";
// 	param.get().diago_proc = 0;
// 	param.get().calculation = "relax";
//     param.get().chg_extrap = "default";
//     param.get().relax_nmax = 0;
//     param.get().basis_type = "pw";
//     param.get().ks_solver = "default";
//     param.get().gamma_only_local = 1;
// 	param.get().scf_thr = -1.0;
// 	param.get().scf_thr_type = -1;
//     param.get().nbndsto_str = "0";
//     param.get().esolver_type = "sdft";
//     param.get().nx = param.get().ny = param.get().nz = 0;
//     param.get().ndx = param.get().ndy = param.get().ndz = 4;
//     // the 2nd calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().nx, 4);
//     EXPECT_EQ(param.get().ny, 4);
//     EXPECT_EQ(param.get().nz, 4);
//     EXPECT_FALSE(GlobalV::double_grid);
//     EXPECT_EQ(param.get().chg_extrap, "first-order");
//     EXPECT_EQ(param.get().vdw_s6, "1.0");
//     EXPECT_EQ(param.get().vdw_s8, "0.722");
//     EXPECT_EQ(param.get().vdw_a1, "1.217");
//     EXPECT_EQ(param.get().vdw_a2,"1.0");
// 	EXPECT_EQ(param.get().vdw_cutoff_radius,"95");
// 	EXPECT_EQ(param.get().exx_hybrid_alpha,"1");
// 	EXPECT_EQ(param.get().exx_real_number,"0");
//         EXPECT_EQ(param.get().exx_ccp_rmesh_times,"5");
// 	EXPECT_EQ(param.get().diago_proc,1);
// 	EXPECT_EQ(param.get().mem_saver,0);
// 	EXPECT_EQ(param.get().cal_force,1);
// 	EXPECT_EQ(param.get().relax_nmax,50);
// 	EXPECT_EQ(param.get().ks_solver,"cg");
// 	EXPECT_EQ(param.get().gamma_only_local,0);
// 	EXPECT_EQ(param.get().bx,1);
// 	EXPECT_EQ(param.get().by,1);
// 	EXPECT_EQ(param.get().bz,1);
// 	EXPECT_DOUBLE_EQ(param.get().scf_thr,1.0e-9);
// 	EXPECT_EQ(param.get().scf_thr_type,1);
// 	EXPECT_EQ(param.get().esolver_type, "ksdft");
// 	//==================================================
// 	// prepare default parameters for the 3rd calling
// 	param.get().vdw_method = "d3_bj";
// 	param.get().vdw_s6 = "default";
// 	param.get().vdw_s8 = "default";
// 	param.get().vdw_a1 = "default";
// 	param.get().vdw_a2 = "default";
// 	param.get().vdw_cutoff_radius = "default";
// 	param.get().calculation = "get_S";
//     param.get().chg_extrap = "default";
//     param.get().basis_type = "pw";
//     param.get().pw_diag_thr = 1.0e-2;
//     param.get().cal_force = 1;
// 	param.get().init_chg = "atomic";
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "cg";
// 	GlobalV::NPROC = 8;
// 	param.get().diago_proc = 1;
//     param.get().nx = param.get().ny = param.get().nz = 4;
//     param.get().ndx = param.get().ndy = param.get().ndz = 6;
//     // the 3rd calling
//     param.get().Default_2();
//     // ^^^^^^^^^^^^^^
//     EXPECT_TRUE(GlobalV::double_grid);
//     EXPECT_EQ(param.get().chg_extrap, "atomic");
//     EXPECT_EQ(param.get().vdw_s6, "1.0");
//     EXPECT_EQ(param.get().vdw_s8, "0.7875");
//     EXPECT_EQ(param.get().vdw_a1,"0.4289");
// 	EXPECT_EQ(param.get().vdw_a2,"4.4407");
// 	EXPECT_EQ(param.get().vdw_cutoff_radius,"95");
// 	EXPECT_EQ(GlobalV::CALCULATION,"nscf");
// 	EXPECT_EQ(param.get().relax_nmax,1);
// 	EXPECT_EQ(param.get().out_stru,0);
// 	EXPECT_DOUBLE_EQ(param.get().pw_diag_thr,1.0e-5);
// 	EXPECT_FALSE(param.get().cal_force);
// 	EXPECT_EQ(param.get().init_chg,"file");
// 	EXPECT_EQ(param.get().diago_proc,8);
// 	//==================================================
// 	// prepare default parameters for the 4th calling
// 	param.get().calculation = "get_pchg";
//     param.get().chg_extrap = "default";
//     param.get().symmetry = "default";
//     param.get().ecutwfc = 10;
//     param.get().ecutrho = 100;
//     param.get().nx = param.get().ny = param.get().nz = 0;
//     param.get().ndx = param.get().ndy = param.get().ndz = 0;
//     GlobalV::double_grid = false;
//     // the 4th calling
//     param.get().Default_2();
//     // ^^^^^^^^^^^^^^
//     EXPECT_TRUE(GlobalV::double_grid);
//     EXPECT_EQ(GlobalV::CALCULATION, "get_pchg");
//     EXPECT_EQ(param.get().relax_nmax, 1);
//     EXPECT_EQ(param.get().out_stru, 0);
//     EXPECT_EQ(param.get().symmetry, "0");
// 	EXPECT_EQ(param.get().out_band[0],0);
// 	EXPECT_EQ(param.get().out_band[1],8);
// 	EXPECT_EQ(param.get().out_proj_band,0);
// 	EXPECT_EQ(param.get().cal_force,0);
// 	EXPECT_EQ(param.get().init_wfc,"file");
// 	EXPECT_EQ(param.get().init_chg,"atomic");
// 	EXPECT_EQ(param.get().chg_extrap,"atomic");
// 	EXPECT_EQ(param.get().out_chg,1);
// 	EXPECT_EQ(param.get().out_dm,0);
// 	EXPECT_EQ(param.get().out_dm1,0);
// 	EXPECT_EQ(param.get().out_pot,0);
// 	//==================================================
// 	// prepare default parameters for the 5th calling
// 	param.get().calculation = "get_wf";
//     param.get().symmetry = "default";
//     param.get().chg_extrap = "default";
//     // the 5th calling
//     param.get().Default_2();
//     // ^^^^^^^^^^^^^^
// 	EXPECT_EQ(GlobalV::CALCULATION,"get_wf");
//     EXPECT_EQ(param.get().relax_nmax, 1);
//     EXPECT_EQ(param.get().symmetry, "0");
//     EXPECT_EQ(param.get().out_stru, 0);
// 	EXPECT_EQ(param.get().out_band[0],0);
// 	EXPECT_EQ(param.get().out_band[1],8);
// 	EXPECT_EQ(param.get().out_proj_band,0);
// 	EXPECT_EQ(param.get().cal_force,0);
// 	EXPECT_EQ(param.get().init_wfc,"file");
// 	EXPECT_EQ(param.get().init_chg,"atomic");
// 	EXPECT_EQ(param.get().chg_extrap,"atomic");
// 	EXPECT_EQ(param.get().out_chg,1);
// 	EXPECT_EQ(param.get().out_dm,0);
// 	EXPECT_EQ(param.get().out_dm1,0);
// 	EXPECT_EQ(param.get().out_pot,0);
// 	//==================================================
// 	// prepare default parameters for the 6th calling
// 	param.get().calculation = "md";
//     param.get().chg_extrap = "default";
//     param.get().mdp.md_nstep = 0;
// 	param.get().out_md_control = 0;
// 	param.get().mdp.md_tlast = -1.0;
// 	param.get().mdp.md_plast = -1.0;
// 	param.get().mdp.md_tfreq = 0;
// 	param.get().mdp.md_pfreq = 0;
// 	param.get().mdp.md_restart = 1;
// 	param.get().mdp.md_type = "npt";
// 	param.get().mdp.md_pmode = "iso";
// 	// the 6th calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
// 	EXPECT_EQ(GlobalV::CALCULATION,"md");
//     EXPECT_EQ(param.get().chg_extrap, "second-order");
//     EXPECT_EQ(param.get().symmetry, "0");
//     EXPECT_EQ(param.get().cal_force, 1);
//     EXPECT_EQ(param.get().mdp.md_nstep,50);
//     EXPECT_EQ(param.get().out_level, "m");
//     EXPECT_DOUBLE_EQ(param.get().mdp.md_plast, param.get().mdp.md_pfirst);
//     EXPECT_DOUBLE_EQ(param.get().mdp.md_tfreq,1.0/40/param.get().mdp.md_dt);
// 	EXPECT_DOUBLE_EQ(param.get().mdp.md_pfreq,1.0/400/param.get().mdp.md_dt);
// 	EXPECT_EQ(param.get().init_vel,1);
// 	EXPECT_EQ(param.get().cal_stress,1);
// 	//==================================================
// 	// prepare default parameters for the 7th calling
// 	param.get().calculation = "cell-relax";
//     param.get().chg_extrap = "default";
//     param.get().relax_nmax = 0;
//     // the 7th calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().chg_extrap, "first-order");
//     EXPECT_EQ(param.get().cal_force, 1);
//     EXPECT_EQ(param.get().cal_stress,1);
// 	EXPECT_EQ(param.get().relax_nmax,50);
// 	//==================================================
// 	// prepare default parameters for the 8th calling
// 	param.get().calculation = "test_memory";
//     param.get().chg_extrap = "default";
//     // the 8th calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().chg_extrap, "atomic");
//     EXPECT_EQ(param.get().relax_nmax,1);
// 	//==================================================
// 	// prepare default parameters for the 9th calling
// 	param.get().calculation = "test_neighbour";
//     param.get().chg_extrap = "default";
//     // the 9th calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().chg_extrap, "atomic");
//     EXPECT_EQ(param.get().relax_nmax,1);
// 	//==================================================
// 	// prepare default parameters for the 10th calling
// 	param.get().calculation = "gen_bessel";
//     param.get().chg_extrap = "default";
//     // the 10th calling
// 	param.get().Default_2();
// 	// ^^^^^^^^^^^^^^
//     EXPECT_EQ(param.get().chg_extrap, "atomic");
//     EXPECT_EQ(param.get().relax_nmax,1);
// 	//==================================================
// 	remove("param.get()");
// 	remove("STRU");
// }

// TEST_F(InputTest, Check)
// {
//     param.get().ecutwfc = 20.0;
//     param.get().ecutrho = 10;
//     testing::internal::CaptureStdout();
//     EXPECT_EXIT(param.get().Check(), ::testing::ExitedWithCode(0), "");
//     output = testing::internal::GetCapturedStdout();
//     EXPECT_THAT(output, testing::HasSubstr("ecutrho/ecutwfc must >= 4"));
//     param.get().ecutrho = 100.0;
//     param.get().nx = param.get().ny = param.get().nz = 10;
//     param.get().ndx = param.get().ndy = param.get().ndz = 8;
//     testing::internal::CaptureStdout();
//     EXPECT_EXIT(param.get().Check(), ::testing::ExitedWithCode(0), "");
//     output = testing::internal::GetCapturedStdout();
//     EXPECT_THAT(output, testing::HasSubstr("smooth grids is denser than dense grids"));
//     param.get().ndx = param.get().ndy = param.get().ndz = 11;
//     //
//     param.get().nbands = -1;
//     testing::internal::CaptureStdout();
//     EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("NBANDS must >= 0"));
// 	param.get().nbands = 2;
// 	//
// 	param.get().nb2d = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nb2d must > 0"));
// 	param.get().nb2d = 1;
// 	//
// 	param.get().ntype = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("ntype must > 0"));
// 	param.get().ntype = 1;
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().diago_proc = 2;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("please don't set diago_proc with lcao base"));
// 	param.get().diago_proc = 1;
// 	//
// 	param.get().kspacing[0] = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("kspacing must > 0"));
// 	param.get().kspacing[0] = param.get().kspacing[1] = param.get().kspacing[2] = 0.8;
// 	//
// 	param.get().nelec = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nelec < 0 is not allowed !"));
// 	param.get().nelec = 100;
// 	//
// 	param.get().efield_flag = 0;
// 	param.get().dip_cor_flag = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("dipole correction is not active if efield_flag=false !"));
// 	param.get().dip_cor_flag = 0;
// 	//
// 	param.get().efield_flag = 1;
// 	param.get().gate_flag = 1;
// 	param.get().dip_cor_flag = 0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("gate field cannot be used with efield if dip_cor_flag=false !"));
// 	param.get().gate_flag = 0;
// 	//
// 	param.get().calculation = "nscf";
// 	param.get().out_dos = 3;
// 	param.get().symmetry = "1";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("symmetry can't be used for out_dos==3(Fermi Surface Plotting) by now."));
// 	param.get().symmetry = "0";
// 	param.get().out_dos = 0;
// 	//
// 	param.get().calculation = "get_pchg";
// 	param.get().basis_type = "pw";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculate = get_pchg is only availble for LCAO"));
// 	param.get().basis_type = "lcao";
// 	//
// 	param.get().calculation = "get_wf";
// 	param.get().basis_type = "pw";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculate = get_wf is only availble for LCAO"));
// 	param.get().basis_type = "lcao";
// 	//
// 	param.get().calculation = "md";
// 	param.get().mdp.md_dt = -1.0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("time interval of MD calculation should be set!"));
// 	param.get().mdp.md_dt = 1.0;
// 	//
// 	param.get().mdp.md_type = "msst";
// 	param.get().mdp.msst_qmass = -1.0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("msst_qmass must be greater than 0!"));
// 	param.get().mdp.msst_qmass = 1.0;
// 	//
// 	param.get().esolver_type = "dp";
// 	param.get().mdp.pot_file = "graph.pb";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("Can not find DP model !"));
// 	param.get().esolver_type = "ksdft";
// 	//
// 	param.get().calculation = "gen_bessel";
// 	param.get().basis_type = "lcao";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("to generate descriptors, please use pw basis"));
// 	param.get().basis_type = "pw";
// 	//
// 	param.get().calculation = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("check 'calculation' !"));
// 	param.get().calculation = "scf";
// 	//
// 	param.get().init_chg = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("wrong 'init_chg',not 'atomic', 'file',please check"));
// 	param.get().init_chg = "atomic";
// 	//
// 	param.get().gamma_only_local = 0;
// 	param.get().out_dm = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("out_dm with k-point algorithm is not implemented yet."));
// 	param.get().out_dm = 0;
// 	//
// 	param.get().gamma_only_local = 1;
// 	param.get().out_dm1 = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("out_dm1 is only for multi-k"));
// 	param.get().gamma_only_local = 0;
// 	param.get().out_dm1 = 0;
// 	//
// 	param.get().nbands = 100001;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nbnd >100000, out of range"));
// 	param.get().nbands = 100;
// 	//
// 	param.get().nelec = 2*param.get().nbands + 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nelec > 2*nbnd , bands not enough!"));
// 	param.get().nelec = param.get().nbands;
// 	//
// 	param.get().nspin = 3;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nspin does not equal to 1, 2, or 4!"));
// 	param.get().nspin = 1;
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "genelpa";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("genelpa can not be used with plane wave basis."));
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "scalapack_gvx";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("scalapack_gvx can not be used with plane wave basis."));
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "lapack";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("lapack can not be used with plane wave basis."));
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("please check the ks_solver parameter!"));
// 	param.get().ks_solver = "cg";
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().gamma_only = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("gamma_only not implemented for plane wave now."));
// 	param.get().gamma_only = 0;
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().out_proj_band = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("out_proj_band not implemented for plane wave now."));
// 	param.get().out_proj_band = 0;
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().out_dos = 3;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("Fermi Surface Plotting not implemented for plane wave now."));
// 	param.get().out_dos = 0;
// 	//
// 	param.get().basis_type = "pw";
// 	param.get().sc_mag_switch = 3;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("Non-colliner Spin-constrained DFT not implemented for plane wave now."));
// 	param.get().sc_mag_switch = 0;
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "cg";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("not ready for cg method in lcao ."));
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "genelpa";
// #ifndef __MPI
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("genelpa can not be used for series version."));
// #endif
// #ifndef __ELPA
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("Can not use genelpa if abacus is not compiled with ELPA. Please change
// ks_solver to scalapack_gvx.")); #endif
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "scalapack_gvx";
// #ifndef __MPI
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("scalapack_gvx can not be used for series version."));
// #endif
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "lapack";
// #ifdef __MPI
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("ks_solver=lapack is not an option for parallel version of ABACUS (try
// genelpa)")); #endif
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "cusolver";
// #ifndef __MPI
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("Cusolver can not be used for series version."));
// #endif
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("please check the ks_solver parameter!"));
// 	param.get().ks_solver = "genelpa";
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().kpar = 2;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("kpar > 1 has not been supported for lcao calculation."));
// 	param.get().kpar = 1;
// 	//
// 	param.get().basis_type = "lcao";
// 	param.get().out_wfc_lcao = 3;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("out_wfc_lcao must be 0, 1, or 2"));
// 	param.get().out_wfc_lcao = 0;
// 	//
// 	param.get().basis_type = "lcao_in_pw";
// 	param.get().ks_solver = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("LCAO in plane wave can only done with lapack."));
// 	param.get().ks_solver = "default";
// 	//
// 	param.get().basis_type = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("please check the basis_type parameter!"));
// 	param.get().basis_type = "pw";
// 	//
// 	param.get().relax_method = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("relax_method can only be sd, cg, bfgs or cg_bfgs."));
// 	param.get().relax_method = "cg";
// 	//
// 	param.get().bx = 11; param.get().by = 1; param.get().bz = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("bx, or by, or bz is larger than 10!"));
// 	param.get().bx = 1;
// 	//
// 	param.get().vdw_method = "d2";
// 	param.get().vdw_C6_unit = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_C6_unit must be Jnm6/mol or eVA6"));
// 	param.get().vdw_C6_unit = "eVA6";
// 	//
// 	param.get().vdw_R0_unit = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_R0_unit must be A or Bohr"));
// 	param.get().vdw_R0_unit = "A";
// 	//
// 	param.get().vdw_cutoff_type = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_type must be radius or period"));
// 	param.get().vdw_cutoff_type = "radius";
// 	//
// 	param.get().vdw_cutoff_period.x = 0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_period <= 0 is not allowd"));
// 	param.get().vdw_cutoff_period.x = 3;
// 	//
// 	param.get().vdw_cutoff_radius = "0";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_radius <= 0 is not allowd"));
// 	param.get().vdw_cutoff_radius = "1.0";
// 	//
// 	param.get().vdw_radius_unit = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_radius_unit must be A or Bohr"));
// 	param.get().vdw_radius_unit = "A";
// 	//
// 	param.get().vdw_cn_thr = 0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_cn_thr <= 0 is not allowd"));
// 	param.get().vdw_cn_thr = 1.0;
// 	//
// 	param.get().vdw_cn_thr_unit = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("vdw_cn_thr_unit must be A or Bohr"));
// 	param.get().vdw_cn_thr_unit = "A";
// 	//
// 	param.get().dft_functional = "scan0";
// 	param.get().exx_hybrid_alpha = "1.25";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("must 0 <= exx_hybrid_alpha <= 1"));
// 	param.get().exx_hybrid_alpha = "0.25";
// 	//
// 	param.get().exx_hybrid_step = 0;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("must exx_hybrid_step > 0"));
// 	param.get().exx_hybrid_step = 1;
// 	//
// 	param.get().exx_ccp_rmesh_times = "-1";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("must exx_ccp_rmesh_times >= 1"));
// 	param.get().exx_ccp_rmesh_times = "1.5";
//     //
//     param.get().rpa = true;
//     param.get().rpa_ccp_rmesh_times = -1;
//     testing::internal::CaptureStdout();
//     EXPECT_EXIT(param.get().Check(), ::testing::ExitedWithCode(0), "");
//     output = testing::internal::GetCapturedStdout();
//     EXPECT_THAT(output, testing::HasSubstr("must rpa_ccp_rmesh_times >= 1"));
//     param.get().rpa_ccp_rmesh_times = 10.0;
// 	//
// 	param.get().exx_distribute_type = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("exx_distribute_type must be htime or kmeans2 or kmeans1"));
// 	param.get().exx_distribute_type = "htime";
// 	//
// 	param.get().dft_functional = "opt_orb";
// 	param.get().exx_opt_orb_lmax = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_lmax must >=0"));
// 	param.get().exx_opt_orb_lmax = 0;
// 	//
// 	param.get().exx_opt_orb_ecut = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_ecut must >=0"));
// 	param.get().exx_opt_orb_ecut = 0;
// 	//
// 	param.get().exx_opt_orb_tolerence = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_tolerence must >=0"));
// 	param.get().exx_opt_orb_tolerence = 0;
// 	//
// 	param.get().berry_phase = 1;
// 	param.get().basis_type = "lcao_in_pw";
// 	param.get().ks_solver = "lapack";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set basis_type = pw or lcao"));
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "cg";
// 	//
// 	param.get().calculation = "scf";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set calculation = nscf"));
// 	param.get().calculation = "nscf";
// 	//
// 	param.get().gdir = 4;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set gdir = 1 or 2 or 3"));
// 	param.get().gdir = 3;
// 	param.get().berry_phase = 0;
// 	//
// 	param.get().towannier90 = 1;
// 	// due to the repair of lcao_in_pw, original warning has been deprecated, 2023/12/23, ykhuang
// 	// param.get().basis_type = "lcao_in_pw";
// 	// param.get().ks_solver = "lapack";
// 	// testing::internal::CaptureStdout();
// 	// EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	// output = testing::internal::GetCapturedStdout();
// 	// EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set basis_type = pw or lcao"));
// 	param.get().basis_type = "pw";
// 	param.get().ks_solver = "cg";
// 	//
// 	param.get().calculation = "scf";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set calculation = nscf"));
// 	param.get().calculation = "nscf";
// 	//
// 	param.get().nspin = 2;
// 	param.get().wannier_spin = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set wannier_spin = up or down"));
// 	param.get().wannier_spin = "up";
// 	param.get().towannier90 = 0;
// 	//
// 	param.get().read_file_dir = "arbitrary";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("please set right files directory for reading in."));
// 	param.get().read_file_dir = "auto";
// 	/*
// 	// Start to check deltaspin parameters
// 	param.get().sc_mag_switch = 1;
// 	param.get().sc_file = "none";
// 	param.get().basis_type = "lcao";
// 	param.get().ks_solver = "genelpa";
// 	// warning 1 of Deltaspin
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("sc_file (json format) must be set when sc_mag_switch > 0"));
// 	// warning 2 of Deltaspin
// 	param.get().sc_file = "sc.json";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("sc_file does not exist"));
// 	param.get().sc_file = "./support/sc.json";
// 	// warning 3 of Deltaspin
// 	param.get().nspin = 1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nspin must be 2 or 4 when sc_mag_switch > 0"));
// 	param.get().nspin = 4;
// 	// warning 4 of Deltaspin
// 	param.get().calculation = "nscf";
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("calculation must be scf when sc_mag_switch > 0"));
// 	param.get().calculation = "scf";
// 	// warning 5 of Deltaspin
// 	param.get().sc_thr = -1;
// 		testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("sc_thr must > 0"));
// 	param.get().sc_thr = 1e-6;
// 	// warning 6 of Deltaspin
// 	param.get().nsc = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nsc must > 0"));
// 	param.get().nsc = 100;
// 	// warning 7 of Deltaspin
// 	param.get().nsc_min = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nsc_min must > 0"));
// 	param.get().nsc_min = 2;
// 	// warning 8 of Deltapsin
//     param.get().alpha_trial = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("alpha_trial must > 0"));
// 	param.get().alpha_trial = 0.01;
// 	// warning 9 of Deltapsin
//     param.get().sccut = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("sccut must > 0"));
// 	param.get().sccut = 3.0;
// 	// warning 10 of Deltaspin
// 	param.get().sc_scf_nmin = -1;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("sc_scf_nmin must >= 2"));
// 	param.get().sc_scf_nmin = 2;
// 	// warning 10 of Deltaspin
// 	param.get().nupdown = 4;
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr("nupdown should not be set when sc_mag_switch > 0"));
// 	param.get().nupdown = 0;
//     // restore to default values
//     param.get().nspin = 1;
// 	param.get().sc_file = "none";
// 	param.get().sc_mag_switch = 0;
// 	param.get().ks_solver = "default";
// 	param.get().basis_type = "pw";
// 	// End of checking Deltaspin parameters
// 	*/

// 	/*
// 	testing::internal::CaptureStdout();
// 	EXPECT_EXIT(param.get().Check(),::testing::ExitedWithCode(0), "");
// 	output = testing::internal::GetCapturedStdout();
// 	EXPECT_THAT(output,testing::HasSubstr(""));
// 	*/
// }

// bool strcmp_inbuilt(const std::string& str1, const std::string& str2)
// {
// 	if(str1.size() != str2.size())
// 		return false;
// 	for(int i=0; i<str1.size(); i++)
// 	{
// 		if(str1[i] != str2[i])
// 			return false;
// 	}
// 	return true;
// }

// TEST_F(InputTest, ReadValue2stdvector)
// {
// 	std::string input_file = "./support/INPUT_list";
// 	std::ifstream ifs(input_file);
// 	std::string word;
// 	std::vector<int> value;
// 	while(!ifs.eof())
// 	{
// 		ifs >> word;
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case0"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 1);
// 			EXPECT_EQ(value[0], 7);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case1"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 1);
// 			EXPECT_EQ(value[0], 7);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case2"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 1);
// 			EXPECT_EQ(value[0], 7);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case3"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 1);
// 			EXPECT_EQ(value[0], 7);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case4"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 1);
// 			EXPECT_EQ(value[0], 7);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case5"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 4);
// 			EXPECT_EQ(value[0], 7);
// 			EXPECT_EQ(value[1], 8);
// 			EXPECT_EQ(value[2], 9);
// 			EXPECT_EQ(value[3], 10);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case6"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 4);
// 			EXPECT_EQ(value[0], 7);
// 			EXPECT_EQ(value[1], 8);
// 			EXPECT_EQ(value[2], 9);
// 			EXPECT_EQ(value[3], 10);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case7"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 4);
// 			EXPECT_EQ(value[0], 7);
// 			EXPECT_EQ(value[1], 8);
// 			EXPECT_EQ(value[2], 9);
// 			EXPECT_EQ(value[3], 10);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case8"))
// 		{
// 			value.clear(); value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, value);
// 			EXPECT_EQ(value.size(), 4);
// 			EXPECT_EQ(value[0], 7);
// 			EXPECT_EQ(value[1], 8);
// 			EXPECT_EQ(value[2], 9);
// 			EXPECT_EQ(value[3], 10);
// 		}
// 		std::vector<std::string> str_value;
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case9"))
// 		{
// 			str_value.clear(); str_value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, str_value);
// 			EXPECT_EQ(str_value.size(), 1);
// 			EXPECT_EQ(str_value[0], "string1");
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case10"))
// 		{
// 			str_value.clear(); str_value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, str_value);
// 			EXPECT_EQ(str_value.size(), 4);
// 			EXPECT_EQ(str_value[0], "string1");
// 			EXPECT_EQ(str_value[1], "string2");
// 			EXPECT_EQ(str_value[2], "string3");
// 			EXPECT_EQ(str_value[3], "string4");
// 		}
// 		std::vector<double> double_value;
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case11"))
// 		{
// 			double_value.clear(); double_value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, double_value);
// 			EXPECT_EQ(double_value.size(), 1);
// 			EXPECT_EQ(double_value[0], 1.23456789);
// 		}
// 		if(strcmp_inbuilt(word, "bessel_nao_rcut_case12"))
// 		{
// 			double_value.clear(); double_value.shrink_to_fit();
// 			param.get().read_value2stdvector(ifs, double_value);
// 			EXPECT_EQ(double_value.size(), 4);
// 			EXPECT_EQ(double_value[0], -1.23456789);
// 			EXPECT_EQ(double_value[1], 2.3456789);
// 			EXPECT_EQ(double_value[2], -3.456789);
// 			EXPECT_EQ(double_value[3], 4.56789);
// 		}
// 	}
// }
// #undef private

// class ReadKSpacingTest : public ::testing::Test {
// protected:
//     void SetUp() override
// 	{
//         // create a temporary file for testing
//         char tmpname[] = "tmpfile.tmp";
//         int fd = mkstemp(tmpname);
//         tmpfile = tmpname;
//         std::ofstream ofs(tmpfile);
//         ofs << "1.0"; // valid input
//         ofs.close();
//     }

//     void TearDown() override {
//         close(fd);
//         unlink(tmpfile.c_str());
//     }

//     std::string tmpfile;
//     int fd;
// };

// TEST_F(ReadKSpacingTest, ValidInputOneValue) {
//     std::ifstream ifs(tmpfile);
//     EXPECT_NO_THROW(param.get().read_kspacing(ifs));
//     EXPECT_EQ(param.get().kspacing[0], 1.0);
//     EXPECT_EQ(param.get().kspacing[1], 1.0);
//     EXPECT_EQ(param.get().kspacing[2], 1.0);
// }

// TEST_F(ReadKSpacingTest, ValidInputThreeValue) {
// 	std::ofstream ofs(tmpfile);
//     ofs << "1.0 2.0 3.0"; // invalid input
//     ofs.close();

//     std::ifstream ifs(tmpfile);
//     EXPECT_NO_THROW(param.get().read_kspacing(ifs));
//     EXPECT_EQ(param.get().kspacing[0], 1.0);
//     EXPECT_EQ(param.get().kspacing[1], 2.0);
//     EXPECT_EQ(param.get().kspacing[2], 3.0);
// }

// TEST_F(ReadKSpacingTest, InvalidInput) {
//     std::ofstream ofs(tmpfile);
//     ofs << "1.0 2.0"; // invalid input
//     ofs.close();

//     std::ifstream ifs(tmpfile);
// 	testing::internal::CaptureStdout();
// 	param.get().read_kspacing(ifs);
// 	std::string output;
// 	output = testing::internal::GetCapturedStdout();
//     EXPECT_TRUE(ifs.fail());
// }

#include "module_base/tool_quit.h"
#include "module_io/read_input.h"
#include "module_parameter/parameter.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
// #ifdef __MPI
#include "module_base/parallel_global.h"
#include "module_basis/module_pw/test/test_tool.h"
#include "mpi.h"
// #endif
/************************************************
 *  unit test of read_input_test.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ParaRead:
 *     - read INPUT file and STRU file
 *   - Check:
 *    - check_mode = true
 */

class InputParaTest : public testing::Test
{
  protected:
};

// #ifdef __MPI
TEST_F(InputParaTest, ParaRead)
{
    ModuleIO::ReadInput readinput(GlobalV::MY_RANK);
    Parameter param;
    readinput.read_parameters(param, "./support/INPUT");
    EXPECT_EQ(param.get().suffix, "autotest");
    EXPECT_EQ(param.get().stru_file, "./support/STRU");
    EXPECT_EQ(param.get().kpoint_file, "KPT");
    EXPECT_EQ(param.get().pseudo_dir, "../../PP_ORB/");
    EXPECT_EQ(param.get().orbital_dir, "../../PP_ORB/");
    EXPECT_EQ(param.get().read_file_dir, "OUT.autotest/");
    EXPECT_EQ(param.get().wannier_card, "none");
    EXPECT_EQ(param.get().latname, "none");
    EXPECT_EQ(param.get().calculation, "scf");
    EXPECT_EQ(param.get().esolver_type, "ksdft");
    EXPECT_DOUBLE_EQ(param.get().pseudo_rcut, 15.0);
    EXPECT_FALSE(param.get().pseudo_mesh);
    EXPECT_EQ(param.get().ntype, 1);
    EXPECT_EQ(param.get().nbands, 8);
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
    EXPECT_EQ(param.get().cond_dt, 0.07);
    EXPECT_EQ(param.get().cond_dtbatch, 2);
    EXPECT_DOUBLE_EQ(param.get().cond_fwhm, 0.3);
    EXPECT_TRUE(param.get().cond_nonlocal);
    EXPECT_FALSE(param.get().berry_phase);
    EXPECT_EQ(param.get().gdir, 3);
    EXPECT_FALSE(param.get().towannier90);
    EXPECT_EQ(param.get().nnkpfile, "seedname.nnkp");
    EXPECT_EQ(param.get().wannier_spin, "up");
    EXPECT_EQ(param.get().wannier_method, 1);
    EXPECT_TRUE(param.get().out_wannier_amn);
    EXPECT_TRUE(param.get().out_wannier_mmn);
    EXPECT_TRUE(param.get().out_wannier_unk);
    EXPECT_TRUE(param.get().out_wannier_eig);
    EXPECT_TRUE(param.get().out_wannier_wvfn_formatted);
    EXPECT_DOUBLE_EQ(param.get().kspacing[0], 0.0);
    EXPECT_DOUBLE_EQ(param.get().kspacing[1], 0.0);
    EXPECT_DOUBLE_EQ(param.get().kspacing[2], 0.0);
    EXPECT_DOUBLE_EQ(param.get().min_dist_coef, 0.2);
    EXPECT_EQ(param.get().dft_functional, "hse");
    EXPECT_DOUBLE_EQ(param.get().xc_temperature, 0.0);
    EXPECT_EQ(param.get().nspin, 1);
    EXPECT_DOUBLE_EQ(param.get().nelec, 0.0);
    EXPECT_EQ(param.get().lmaxmax, 2);
    EXPECT_EQ(param.get().basis_type, "lcao");
    EXPECT_EQ(param.get().ks_solver, "genelpa");
    EXPECT_DOUBLE_EQ(param.get().search_radius, -1.0);
    EXPECT_TRUE(param.get().search_pbc);
    EXPECT_EQ(param.get().symmetry, "1");
    EXPECT_FALSE(param.get().init_vel);
    EXPECT_DOUBLE_EQ(param.get().symmetry_prec, 1.0e-6);
    EXPECT_TRUE(param.get().symmetry_autoclose);
    EXPECT_EQ(param.get().cal_force, 0);
    EXPECT_NEAR(param.get().force_thr, 1.0e-3, 1.0e-7);
    EXPECT_DOUBLE_EQ(param.get().force_thr_ev2, 0);
    EXPECT_DOUBLE_EQ(param.get().stress_thr, 1.0e-2);
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
    EXPECT_TRUE(param.get().sup.out_md_control);
    EXPECT_TRUE(param.get().relax_new);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w1, 0.01);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_w2, 0.5);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmax, 0.8);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_rmin, 1e-5);
    EXPECT_DOUBLE_EQ(param.get().relax_bfgs_init, 0.5);
    EXPECT_DOUBLE_EQ(param.get().relax_scale_force, 0.5);
    EXPECT_EQ(param.get().nbspline, -1);
    EXPECT_TRUE(param.get().gamma_only);
    EXPECT_TRUE(param.get().sup.gamma_only_local);
    EXPECT_DOUBLE_EQ(param.get().ecutwfc, 20.0);
    EXPECT_DOUBLE_EQ(param.get().erf_ecut, 20.0);
    EXPECT_DOUBLE_EQ(param.get().erf_height, 20.0);
    EXPECT_DOUBLE_EQ(param.get().erf_sigma, 4.0);
    EXPECT_DOUBLE_EQ(param.get().ecutrho, 80);
    EXPECT_EQ(param.get().fft_mode, 0);
    EXPECT_EQ(param.get().sup.ncx, 0);
    EXPECT_EQ(param.get().sup.ncy, 0);
    EXPECT_EQ(param.get().sup.ncz, 0);
    EXPECT_EQ(param.get().nx, 0);
    EXPECT_EQ(param.get().ny, 0);
    EXPECT_EQ(param.get().nz, 0);
    EXPECT_EQ(param.get().bx, 2);
    EXPECT_EQ(param.get().by, 2);
    EXPECT_EQ(param.get().bz, 2);
    EXPECT_EQ(param.get().ndx, 0);
    EXPECT_EQ(param.get().ndy, 0);
    EXPECT_EQ(param.get().ndz, 0);
    EXPECT_EQ(param.get().diago_proc, std::min(GlobalV::NPROC, 4));
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
    EXPECT_NEAR(param.get().scf_thr, 1.0e-8, 1.0e-15);
    EXPECT_EQ(param.get().scf_nmax, 50);
    EXPECT_EQ(param.get().relax_nmax, 1);
    EXPECT_EQ(param.get().out_stru, 0);
    EXPECT_EQ(param.get().smearing_method, "gauss");
    EXPECT_DOUBLE_EQ(param.get().smearing_sigma, 0.002);
    EXPECT_EQ(param.get().mixing_mode, "broyden");
    EXPECT_DOUBLE_EQ(param.get().mixing_beta, 0.7);
    EXPECT_EQ(param.get().mixing_ndim, 8);
    EXPECT_DOUBLE_EQ(param.get().mixing_gg0, 0.00);
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
    EXPECT_EQ(param.get().out_pot, 2);
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
    EXPECT_EQ(param.get().out_app_flag, 0);
    EXPECT_EQ(param.get().out_mat_r, 0);
    EXPECT_FALSE(param.get().out_wfc_lcao);
    EXPECT_FALSE(param.get().out_alllog);
    EXPECT_DOUBLE_EQ(param.get().dos_emin_ev, -15);
    EXPECT_DOUBLE_EQ(param.get().dos_emax_ev, 15);
    EXPECT_DOUBLE_EQ(param.get().dos_edelta_ev, 0.01);
    EXPECT_DOUBLE_EQ(param.get().dos_scale, 0.01);
    EXPECT_DOUBLE_EQ(param.get().dos_sigma, 0.07);
    EXPECT_FALSE(param.get().out_element_info);
    EXPECT_DOUBLE_EQ(param.get().lcao_ecut, 20);
    EXPECT_DOUBLE_EQ(param.get().lcao_dk, 0.01);
    EXPECT_DOUBLE_EQ(param.get().lcao_dr, 0.01);
    EXPECT_DOUBLE_EQ(param.get().lcao_rmax, 30);
    EXPECT_TRUE(param.get().bessel_nao_smooth);
    EXPECT_DOUBLE_EQ(param.get().bessel_nao_sigma, 0.1);
    EXPECT_EQ(std::stod(param.get().bessel_nao_ecut), 20);
    EXPECT_DOUBLE_EQ(param.get().sup.bessel_nao_rcut, 6.0);
    EXPECT_DOUBLE_EQ(param.get().bessel_nao_tolerence, 1E-12);
    EXPECT_EQ(param.get().bessel_descriptor_lmax, 2);
    EXPECT_TRUE(param.get().bessel_descriptor_smooth);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_sigma, 0.1);
    EXPECT_EQ(std::stod(param.get().bessel_descriptor_ecut), 20);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_rcut, 6.0);
    EXPECT_DOUBLE_EQ(param.get().bessel_descriptor_tolerence, 1E-12);
    EXPECT_FALSE(param.get().efield_flag);
    EXPECT_FALSE(param.get().dip_cor_flag);
    EXPECT_EQ(param.get().efield_dir, 2);
    EXPECT_DOUBLE_EQ(param.get().efield_pos_max, 0.5);
    EXPECT_DOUBLE_EQ(param.get().efield_pos_dec, 0.1);
    EXPECT_DOUBLE_EQ(param.get().efield_amp, 0.0);
    EXPECT_FALSE(param.get().gate_flag);
    EXPECT_DOUBLE_EQ(param.get().zgate, 0.5);
    EXPECT_FALSE(param.get().relax);
    EXPECT_FALSE(param.get().block);
    EXPECT_DOUBLE_EQ(param.get().block_down, 0.45);
    EXPECT_DOUBLE_EQ(param.get().block_up, 0.55);
    EXPECT_DOUBLE_EQ(param.get().block_height, 0.1);
    EXPECT_EQ(param.get().vdw_method, "d2");
    EXPECT_EQ(std::stod(param.get().vdw_s6), 0.75);
    EXPECT_EQ(param.get().vdw_s8, "default");
    EXPECT_EQ(param.get().vdw_a1, "default");
    EXPECT_EQ(param.get().vdw_a2, "default");
    EXPECT_DOUBLE_EQ(param.get().vdw_d, 20);
    EXPECT_FALSE(param.get().vdw_abc);
    EXPECT_EQ(std::stod(param.get().vdw_cutoff_radius), 56.6918);
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
    EXPECT_EQ(std::stod(param.get().exx_hybrid_alpha), 0.25);
    EXPECT_EQ(param.get().exx_real_number, "1");
    EXPECT_DOUBLE_EQ(param.get().exx_hse_omega, 0.11);
    EXPECT_TRUE(param.get().exx_separate_loop);
    EXPECT_EQ(param.get().exx_hybrid_step, 100);
    EXPECT_DOUBLE_EQ(param.get().exx_lambda, 0.3);
    EXPECT_DOUBLE_EQ(param.get().exx_mixing_beta, 1.0);
    EXPECT_DOUBLE_EQ(param.get().exx_pca_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_c_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_v_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_dm_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_schwarz_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_c_grad_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_v_grad_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_force_threshold, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_cauchy_stress_threshold, 0);
    EXPECT_EQ(param.get().exx_ccp_rmesh_times, "1.5");
    EXPECT_DOUBLE_EQ(param.get().rpa_ccp_rmesh_times, 10.0);
    EXPECT_EQ(param.get().exx_distribute_type, "htime");
    EXPECT_EQ(param.get().exx_opt_orb_lmax, 0);
    EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_ecut, 0.0);
    EXPECT_DOUBLE_EQ(param.get().exx_opt_orb_tolerence, 0.0);
    EXPECT_FALSE(param.get().noncolin);
    EXPECT_FALSE(param.get().lspinorb);
    EXPECT_DOUBLE_EQ(param.get().soc_lambda, 1.0);
    EXPECT_DOUBLE_EQ(param.get().td_force_dt, 0.02);
    EXPECT_EQ(param.get().td_vext, 0);
    EXPECT_EQ(param.get().td_vext_dire[0], 1);
    EXPECT_EQ(param.get().propagator, 0);
    EXPECT_EQ(param.get().td_stype, 0);
    EXPECT_EQ(param.get().td_ttype, "0");
    EXPECT_EQ(param.get().td_tstart, 1);
    EXPECT_EQ(param.get().td_tend, 1000);
    EXPECT_EQ(param.get().td_lcut1, 0.05);
    EXPECT_EQ(param.get().td_lcut2, 0.95);
    EXPECT_EQ(param.get().td_gauss_amp, "0.25");
    EXPECT_EQ(param.get().td_gauss_freq, "22.13");
    EXPECT_EQ(param.get().td_gauss_phase, "0.0");
    EXPECT_EQ(param.get().td_gauss_t0, "100.0");
    EXPECT_EQ(param.get().td_gauss_sigma, "30.0");
    EXPECT_EQ(param.get().td_trape_amp, "2.74");
    EXPECT_EQ(param.get().td_trape_freq, "1.60");
    EXPECT_EQ(param.get().td_trape_phase, "0.0");
    EXPECT_EQ(param.get().td_trape_t1, "1875");
    EXPECT_EQ(param.get().td_trape_t2, "5625");
    EXPECT_EQ(param.get().td_trape_t3, "7500");
    EXPECT_EQ(param.get().td_trigo_freq1, "1.164656");
    EXPECT_EQ(param.get().td_trigo_freq2, "0.029116");
    EXPECT_EQ(param.get().td_trigo_phase1, "0.0");
    EXPECT_EQ(param.get().td_trigo_phase2, "0.0");
    EXPECT_EQ(param.get().td_trigo_amp, "2.74");
    EXPECT_EQ(param.get().td_heavi_t0, "100");
    EXPECT_EQ(param.get().td_heavi_amp, "1.0");

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
    EXPECT_EQ(param.get().of_kinetic, "vw");
    EXPECT_EQ(param.get().of_method, "tn");
    EXPECT_EQ(param.get().of_conv, "energy");
    EXPECT_DOUBLE_EQ(param.get().of_tole, 1e-6);
    EXPECT_DOUBLE_EQ(param.get().of_tolp, 1e-5);
    EXPECT_DOUBLE_EQ(param.get().of_tf_weight, 1.);
    EXPECT_DOUBLE_EQ(param.get().of_vw_weight, 1.);
    EXPECT_DOUBLE_EQ(param.get().of_wt_alpha, 0.833333);
    EXPECT_DOUBLE_EQ(param.get().of_wt_beta, 0.833333);
    EXPECT_DOUBLE_EQ(param.get().of_wt_rho0, 1.);
    EXPECT_TRUE(param.get().of_hold_rho0);
    EXPECT_DOUBLE_EQ(param.get().of_lkt_a, 1.3);
    EXPECT_FALSE(param.get().of_full_pw);
    EXPECT_EQ(param.get().of_full_pw_dim, 0);
    EXPECT_FALSE(param.get().of_read_kernel);
    EXPECT_EQ(param.get().of_kernel_file, "WTkernel.txt");
    EXPECT_EQ(param.get().device, "cpu");
    EXPECT_EQ(param.get().sup.ncx, 0);
    EXPECT_EQ(param.get().sup.ncy, 0);
    EXPECT_EQ(param.get().sup.ncz, 0);
    EXPECT_NEAR(param.get().force_thr_ev, 0.025711245953622324, 1e-8);
    EXPECT_DOUBLE_EQ(param.get().sup.hubbard_u[0], 0);
    EXPECT_EQ(param.get().orbital_corr[0], -1);
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
    EXPECT_DOUBLE_EQ(param.get().ref_cell_factor, 1.2);
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
    EXPECT_FALSE(param.get().mdp.dump_force);
    EXPECT_FALSE(param.get().mdp.dump_vel);
    EXPECT_FALSE(param.get().mdp.dump_virial);
    EXPECT_EQ(param.get().sc_mag_switch, 0);
    EXPECT_TRUE(param.get().decay_grad_switch);
    EXPECT_DOUBLE_EQ(param.get().sc_thr, 1e-4);
    EXPECT_EQ(param.get().nsc, 50);
    EXPECT_EQ(param.get().nsc_min, 4);
    EXPECT_EQ(param.get().sc_scf_nmin, 4);
    EXPECT_DOUBLE_EQ(param.get().alpha_trial, 0.02);
    EXPECT_DOUBLE_EQ(param.get().sccut, 4.0);
    EXPECT_EQ(param.get().sc_file, "sc.json");
}

TEST_F(InputParaTest, Check)
{
    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream emptyfile("./empty_INPUT");
        emptyfile << "INPUT_PARAMETERS                \n";
        emptyfile << "stru_file    ./support/STRU     \n";
        emptyfile.close();
    }
    ModuleIO::ReadInput::check_mode = true;
    ModuleIO::ReadInput readinput(GlobalV::MY_RANK);
    Parameter param;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(readinput.read_parameters(param, "./empty_INPUT"), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("INPUT parameters have been successfully checked!"));
    if (GlobalV::MY_RANK == 0)
    {
        EXPECT_TRUE(std::remove("./empty_INPUT") == 0);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
// #endif
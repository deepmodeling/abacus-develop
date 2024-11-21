#include <cstdio>
#include <fstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/tool_quit.h"
#include "module_io/read_input.h"
#include "module_parameter/parameter.h"
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
    EXPECT_EQ(param.inp.suffix, "autotest");
    EXPECT_EQ(param.inp.stru_file, "./support/STRU");
    EXPECT_EQ(param.inp.kpoint_file, "KPT");
    EXPECT_EQ(param.inp.pseudo_dir, "../../PP_ORB/");
    EXPECT_EQ(param.inp.orbital_dir, "../../PP_ORB/");
    EXPECT_EQ(param.inp.read_file_dir, "OUT.autotest/");
    EXPECT_EQ(param.inp.wannier_card, "none");
    EXPECT_EQ(param.inp.latname, "none");
    EXPECT_EQ(param.inp.calculation, "scf");
    EXPECT_EQ(param.inp.esolver_type, "ksdft");
    EXPECT_DOUBLE_EQ(param.inp.pseudo_rcut, 15.0);
    EXPECT_FALSE(param.inp.pseudo_mesh);
    EXPECT_EQ(param.inp.ntype, 1);
    EXPECT_EQ(param.inp.nbands, 8);
    EXPECT_EQ(param.inp.nbands_sto, 256);
    EXPECT_EQ(param.inp.nbands_istate, 5);
    EXPECT_EQ(param.inp.bands_to_print.size(), 0);
    EXPECT_EQ(param.inp.out_pchg.size(), 0);
    EXPECT_EQ(param.inp.out_wfc_norm.size(), 0);
    EXPECT_EQ(param.inp.out_wfc_re_im.size(), 0);
    EXPECT_FALSE(param.inp.if_separate_k);
    EXPECT_EQ(param.inp.pw_seed, 1);
    EXPECT_EQ(param.inp.emin_sto, 0.0);
    EXPECT_EQ(param.inp.emax_sto, 0.0);
    EXPECT_EQ(param.inp.nche_sto, 100);
    EXPECT_EQ(param.inp.seed_sto, 0);
    EXPECT_EQ(param.inp.initsto_ecut, 0.0);
    EXPECT_EQ(param.inp.bndpar, 1);
    EXPECT_EQ(param.inp.kpar, 1);
    EXPECT_EQ(param.inp.initsto_freq, 0);
    EXPECT_EQ(param.inp.method_sto, 2);
    EXPECT_EQ(param.inp.npart_sto, 1);
    EXPECT_FALSE(param.inp.cal_cond);
    EXPECT_EQ(param.inp.dos_nche, 100);
    EXPECT_DOUBLE_EQ(param.inp.cond_che_thr, 1e-8);
    EXPECT_DOUBLE_EQ(param.inp.cond_dw, 0.1);
    EXPECT_DOUBLE_EQ(param.inp.cond_wcut, 10);
    EXPECT_EQ(param.inp.cond_dt, 0.07);
    EXPECT_EQ(param.inp.cond_dtbatch, 2);
    EXPECT_DOUBLE_EQ(param.inp.cond_fwhm, 0.3);
    EXPECT_TRUE(param.inp.cond_nonlocal);
    EXPECT_FALSE(param.inp.berry_phase);
    EXPECT_EQ(param.inp.ocp_kb.size(), 2);
    EXPECT_EQ(param.inp.ocp_kb[0], 1);
    EXPECT_EQ(param.inp.ocp_kb[1], 1);
    EXPECT_EQ(param.inp.gdir, 3);
    EXPECT_FALSE(param.inp.towannier90);
    EXPECT_EQ(param.inp.nnkpfile, "seedname.nnkp");
    EXPECT_EQ(param.inp.wannier_spin, "up");
    EXPECT_EQ(param.inp.wannier_method, 1);
    EXPECT_TRUE(param.inp.out_wannier_amn);
    EXPECT_TRUE(param.inp.out_wannier_mmn);
    EXPECT_TRUE(param.inp.out_wannier_unk);
    EXPECT_TRUE(param.inp.out_wannier_eig);
    EXPECT_TRUE(param.inp.out_wannier_wvfn_formatted);
    EXPECT_DOUBLE_EQ(param.inp.kspacing[0], 0.0);
    EXPECT_DOUBLE_EQ(param.inp.kspacing[1], 0.0);
    EXPECT_DOUBLE_EQ(param.inp.kspacing[2], 0.0);
    EXPECT_DOUBLE_EQ(param.inp.min_dist_coef, 0.2);
    EXPECT_EQ(param.inp.dft_functional, "hse");
    EXPECT_DOUBLE_EQ(param.inp.xc_temperature, 0.0);
    EXPECT_EQ(param.inp.nspin, 1);
    EXPECT_DOUBLE_EQ(param.inp.nelec, 0.0);
    EXPECT_EQ(param.inp.lmaxmax, 2);
    EXPECT_EQ(param.inp.basis_type, "lcao");
    EXPECT_EQ(param.inp.ks_solver, "genelpa");
    EXPECT_DOUBLE_EQ(param.inp.search_radius, -1.0);
    EXPECT_TRUE(param.inp.search_pbc);
    EXPECT_EQ(param.inp.symmetry, "1");
    EXPECT_FALSE(param.inp.init_vel);
    EXPECT_DOUBLE_EQ(param.inp.symmetry_prec, 1.0e-6);
    EXPECT_TRUE(param.inp.symmetry_autoclose);
    EXPECT_EQ(param.inp.cal_force, 0);
    EXPECT_NEAR(param.inp.force_thr, 1.0e-3, 1.0e-7);
    EXPECT_DOUBLE_EQ(param.inp.force_thr_ev2, 0);
    EXPECT_DOUBLE_EQ(param.inp.stress_thr, 1.0e-2);
    EXPECT_DOUBLE_EQ(param.inp.press1, 0.0);
    EXPECT_DOUBLE_EQ(param.inp.press2, 0.0);
    EXPECT_DOUBLE_EQ(param.inp.press3, 0.0);
    EXPECT_FALSE(param.inp.cal_stress);
    EXPECT_EQ(param.inp.fixed_axes, "None");
    EXPECT_FALSE(param.inp.fixed_ibrav);
    EXPECT_FALSE(param.inp.fixed_atoms);
    EXPECT_EQ(param.inp.relax_method, "cg");
    EXPECT_DOUBLE_EQ(param.inp.relax_cg_thr, 0.5);
    EXPECT_EQ(param.inp.out_level, "ie");
    EXPECT_TRUE(param.globalv.out_md_control);
    EXPECT_TRUE(param.inp.relax_new);
    EXPECT_DOUBLE_EQ(param.inp.relax_bfgs_w1, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.relax_bfgs_w2, 0.5);
    EXPECT_DOUBLE_EQ(param.inp.relax_bfgs_rmax, 0.8);
    EXPECT_DOUBLE_EQ(param.inp.relax_bfgs_rmin, 1e-5);
    EXPECT_DOUBLE_EQ(param.inp.relax_bfgs_init, 0.5);
    EXPECT_DOUBLE_EQ(param.inp.relax_scale_force, 0.5);
    EXPECT_EQ(param.inp.nbspline, -1);
    EXPECT_FALSE(param.globalv.gamma_only_pw);
    EXPECT_TRUE(param.globalv.gamma_only_local);
    EXPECT_DOUBLE_EQ(param.inp.ecutwfc, 20.0);
    EXPECT_DOUBLE_EQ(param.inp.erf_ecut, 20.0);
    EXPECT_DOUBLE_EQ(param.inp.erf_height, 20.0);
    EXPECT_DOUBLE_EQ(param.inp.erf_sigma, 4.0);
    EXPECT_DOUBLE_EQ(param.inp.ecutrho, 80);
    EXPECT_EQ(param.inp.fft_mode, 0);
    EXPECT_EQ(param.globalv.ncx, 0);
    EXPECT_EQ(param.globalv.ncy, 0);
    EXPECT_EQ(param.globalv.ncz, 0);
    EXPECT_EQ(param.inp.nx, 0);
    EXPECT_EQ(param.inp.ny, 0);
    EXPECT_EQ(param.inp.nz, 0);
    EXPECT_EQ(param.inp.bx, 2);
    EXPECT_EQ(param.inp.by, 2);
    EXPECT_EQ(param.inp.bz, 2);
    EXPECT_EQ(param.inp.ndx, 0);
    EXPECT_EQ(param.inp.ndy, 0);
    EXPECT_EQ(param.inp.ndz, 0);
    EXPECT_EQ(param.inp.diago_proc, std::min(GlobalV::NPROC, 4));
    EXPECT_EQ(param.inp.pw_diag_nmax, 50);
    EXPECT_EQ(param.inp.diago_cg_prec, 1);
    EXPECT_EQ(param.inp.pw_diag_ndim, 4);
    EXPECT_DOUBLE_EQ(param.inp.pw_diag_thr, 1.0e-2);
    EXPECT_EQ(param.inp.nb2d, 0);
    EXPECT_EQ(param.inp.nurse, 0);
    EXPECT_EQ(param.inp.t_in_h, 1);
    EXPECT_EQ(param.inp.vl_in_h, 1);
    EXPECT_EQ(param.inp.vnl_in_h, 1);
    EXPECT_EQ(param.inp.vh_in_h, 1);
    EXPECT_EQ(param.inp.vion_in_h, 1);
    EXPECT_EQ(PARAM.inp.test_force, 0);
    EXPECT_EQ(param.inp.test_stress, 0);
    EXPECT_NEAR(param.inp.scf_thr, 1.0e-8, 1.0e-15);
    EXPECT_EQ(param.inp.scf_os_stop, 1);
    EXPECT_NEAR(param.inp.scf_os_thr, -0.02, 1.0e-15);
    EXPECT_EQ(param.inp.scf_os_ndim, 10);
    EXPECT_NEAR(param.inp.scf_ene_thr, 1.0e-6, 1.0e-15);
    EXPECT_EQ(param.inp.scf_nmax, 50);
    EXPECT_EQ(param.inp.relax_nmax, 1);
    EXPECT_EQ(param.inp.out_stru, 0);
    EXPECT_EQ(param.inp.smearing_method, "gauss");
    EXPECT_DOUBLE_EQ(param.inp.smearing_sigma, 0.002);
    EXPECT_EQ(param.inp.mixing_mode, "broyden");
    EXPECT_DOUBLE_EQ(param.inp.mixing_beta, 0.7);
    EXPECT_EQ(param.inp.mixing_ndim, 8);
    EXPECT_DOUBLE_EQ(param.inp.mixing_gg0, 0.00);
    EXPECT_EQ(param.inp.init_wfc, "atomic");
    EXPECT_EQ(param.inp.mem_saver, 0);
    EXPECT_EQ(param.inp.printe, 100);
    EXPECT_EQ(param.inp.init_chg, "atomic");
    EXPECT_EQ(param.inp.chg_extrap, "atomic");
    EXPECT_EQ(param.inp.out_freq_elec, 0);
    EXPECT_EQ(param.inp.out_freq_ion, 0);
    EXPECT_EQ(param.inp.out_chg[0], 0);
    EXPECT_EQ(param.inp.out_chg[1], 3);
    EXPECT_EQ(param.inp.out_elf[0], 0);
    EXPECT_EQ(param.inp.out_elf[1], 3);
    EXPECT_EQ(param.inp.out_dm, 0);
    EXPECT_EQ(param.inp.out_dm1, 0);
    EXPECT_EQ(param.inp.deepks_out_labels, 0);
    EXPECT_EQ(param.inp.deepks_scf, 0);
    EXPECT_EQ(param.inp.deepks_equiv, 0);
    EXPECT_EQ(param.inp.deepks_bandgap, 0);
    EXPECT_EQ(param.inp.deepks_out_unittest, 0);
    EXPECT_EQ(param.inp.out_pot, 2);
    EXPECT_EQ(param.inp.out_wfc_pw, 0);
    EXPECT_EQ(param.inp.out_wfc_r, 0);
    EXPECT_EQ(param.inp.out_dos, 0);
    EXPECT_EQ(param.inp.out_band[0], 0);
    EXPECT_EQ(param.inp.out_band[1], 8);
    EXPECT_EQ(param.inp.out_proj_band, 0);
    EXPECT_EQ(param.inp.out_mat_hs[0], 0);
    EXPECT_EQ(param.inp.out_mat_hs[1], 8);
    EXPECT_EQ(param.inp.out_mat_hs2, 0);
    EXPECT_FALSE(param.inp.out_mat_xc);
    EXPECT_FALSE(param.inp.out_eband_terms);
    EXPECT_EQ(param.inp.out_interval, 1);
    EXPECT_EQ(param.inp.out_app_flag, 0);
    EXPECT_EQ(param.inp.out_mat_r, 0);
    EXPECT_FALSE(param.inp.out_wfc_lcao);
    EXPECT_FALSE(param.inp.out_alllog);
    EXPECT_DOUBLE_EQ(param.inp.dos_emin_ev, -15);
    EXPECT_DOUBLE_EQ(param.inp.dos_emax_ev, 15);
    EXPECT_DOUBLE_EQ(param.inp.dos_edelta_ev, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.dos_scale, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.dos_sigma, 0.07);
    EXPECT_FALSE(param.inp.out_element_info);
    EXPECT_DOUBLE_EQ(param.inp.lcao_ecut, 20);
    EXPECT_DOUBLE_EQ(param.inp.lcao_dk, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.lcao_dr, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.lcao_rmax, 30);
    EXPECT_TRUE(param.inp.bessel_nao_smooth);
    EXPECT_DOUBLE_EQ(param.inp.bessel_nao_sigma, 0.1);
    EXPECT_EQ(std::stod(param.inp.bessel_nao_ecut), 20);
    EXPECT_DOUBLE_EQ(param.inp.bessel_nao_rcuts[0], 6.0);
    EXPECT_DOUBLE_EQ(param.inp.bessel_nao_tolerence, 1E-12);
    EXPECT_EQ(param.inp.bessel_descriptor_lmax, 2);
    EXPECT_TRUE(param.inp.bessel_descriptor_smooth);
    EXPECT_DOUBLE_EQ(param.inp.bessel_descriptor_sigma, 0.1);
    EXPECT_EQ(std::stod(param.inp.bessel_descriptor_ecut), 20);
    EXPECT_DOUBLE_EQ(param.inp.bessel_descriptor_rcut, 6.0);
    EXPECT_DOUBLE_EQ(param.inp.bessel_descriptor_tolerence, 1E-12);
    EXPECT_FALSE(param.inp.efield_flag);
    EXPECT_FALSE(param.inp.dip_cor_flag);
    EXPECT_EQ(param.inp.efield_dir, 2);
    EXPECT_DOUBLE_EQ(param.inp.efield_pos_max, 0.5);
    EXPECT_DOUBLE_EQ(param.inp.efield_pos_dec, 0.1);
    EXPECT_DOUBLE_EQ(param.inp.efield_amp, 0.0);
    EXPECT_FALSE(param.inp.gate_flag);
    EXPECT_DOUBLE_EQ(param.inp.zgate, 0.5);
    EXPECT_FALSE(param.inp.relax);
    EXPECT_FALSE(param.inp.block);
    EXPECT_DOUBLE_EQ(param.inp.block_down, 0.45);
    EXPECT_DOUBLE_EQ(param.inp.block_up, 0.55);
    EXPECT_DOUBLE_EQ(param.inp.block_height, 0.1);
    EXPECT_EQ(param.inp.vdw_method, "d2");
    EXPECT_EQ(std::stod(param.inp.vdw_s6), 0.75);
    EXPECT_EQ(param.inp.vdw_s8, "default");
    EXPECT_EQ(param.inp.vdw_a1, "default");
    EXPECT_EQ(param.inp.vdw_a2, "default");
    EXPECT_DOUBLE_EQ(param.inp.vdw_d, 20);
    EXPECT_FALSE(param.inp.vdw_abc);
    EXPECT_EQ(std::stod(param.inp.vdw_cutoff_radius), 56.6918);
    EXPECT_EQ(param.inp.vdw_radius_unit, "Bohr");
    EXPECT_DOUBLE_EQ(param.inp.vdw_cn_thr, 40.0);
    EXPECT_EQ(param.inp.vdw_cn_thr_unit, "Bohr");
    EXPECT_EQ(param.inp.vdw_C6_file, "default");
    EXPECT_EQ(param.inp.vdw_C6_unit, "Jnm6/mol");
    EXPECT_EQ(param.inp.vdw_R0_file, "default");
    EXPECT_EQ(param.inp.vdw_R0_unit, "A");
    EXPECT_EQ(param.inp.vdw_cutoff_type, "radius");
    EXPECT_EQ(param.inp.vdw_cutoff_period[0], 3);
    EXPECT_EQ(param.inp.vdw_cutoff_period[1], 3);
    EXPECT_EQ(param.inp.vdw_cutoff_period[2], 3);
    EXPECT_EQ(std::stod(param.inp.exx_hybrid_alpha), 0.25);
    EXPECT_EQ(param.inp.exx_real_number, "1");
    EXPECT_DOUBLE_EQ(param.inp.exx_hse_omega, 0.11);
    EXPECT_TRUE(param.inp.exx_separate_loop);
    EXPECT_EQ(param.inp.exx_hybrid_step, 100);
    EXPECT_DOUBLE_EQ(param.inp.exx_lambda, 0.3);
    EXPECT_DOUBLE_EQ(param.inp.exx_mixing_beta, 1.0);
    EXPECT_DOUBLE_EQ(param.inp.exx_pca_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_c_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_v_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_dm_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_schwarz_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_cauchy_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_c_grad_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_v_grad_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_cauchy_force_threshold, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_cauchy_stress_threshold, 0);
    EXPECT_EQ(param.inp.exx_ccp_rmesh_times, "1.5");
    EXPECT_DOUBLE_EQ(param.inp.rpa_ccp_rmesh_times, 10.0);
    EXPECT_EQ(param.inp.exx_distribute_type, "htime");
    EXPECT_EQ(param.inp.exx_opt_orb_lmax, 0);
    EXPECT_DOUBLE_EQ(param.inp.exx_opt_orb_ecut, 0.0);
    EXPECT_DOUBLE_EQ(param.inp.exx_opt_orb_tolerence, 0.0);
    EXPECT_FALSE(param.inp.noncolin);
    EXPECT_FALSE(param.inp.lspinorb);
    EXPECT_DOUBLE_EQ(param.inp.soc_lambda, 1.0);
    EXPECT_DOUBLE_EQ(param.inp.td_force_dt, 0.02);
    EXPECT_EQ(param.inp.td_vext, 0);
    EXPECT_EQ(param.inp.td_vext_dire, "1");
    EXPECT_EQ(param.inp.propagator, 0);
    EXPECT_EQ(param.inp.td_stype, 0);
    EXPECT_EQ(param.inp.td_ttype, "0");
    EXPECT_EQ(param.inp.td_tstart, 1);
    EXPECT_EQ(param.inp.td_tend, 1000);
    EXPECT_EQ(param.inp.td_lcut1, 0.05);
    EXPECT_EQ(param.inp.td_lcut2, 0.95);
    EXPECT_EQ(param.inp.td_gauss_amp, "0.25");
    EXPECT_EQ(param.inp.td_gauss_freq, "22.13");
    EXPECT_EQ(param.inp.td_gauss_phase, "0.0");
    EXPECT_EQ(param.inp.td_gauss_t0, "100.0");
    EXPECT_EQ(param.inp.td_gauss_sigma, "30.0");
    EXPECT_EQ(param.inp.td_trape_amp, "2.74");
    EXPECT_EQ(param.inp.td_trape_freq, "1.60");
    EXPECT_EQ(param.inp.td_trape_phase, "0.0");
    EXPECT_EQ(param.inp.td_trape_t1, "1875");
    EXPECT_EQ(param.inp.td_trape_t2, "5625");
    EXPECT_EQ(param.inp.td_trape_t3, "7500");
    EXPECT_EQ(param.inp.td_trigo_freq1, "1.164656");
    EXPECT_EQ(param.inp.td_trigo_freq2, "0.029116");
    EXPECT_EQ(param.inp.td_trigo_phase1, "0.0");
    EXPECT_EQ(param.inp.td_trigo_phase2, "0.0");
    EXPECT_EQ(param.inp.td_trigo_amp, "2.74");
    EXPECT_EQ(param.inp.td_heavi_t0, "100");
    EXPECT_EQ(param.inp.td_heavi_amp, "1.0");

    EXPECT_EQ(param.inp.out_dipole, 0);
    EXPECT_EQ(param.inp.out_efield, 0);
    EXPECT_EQ(param.inp.td_print_eij, -1.0);
    EXPECT_EQ(param.inp.td_edm, 0);
    EXPECT_DOUBLE_EQ(param.inp.cell_factor, 1.2);
    EXPECT_EQ(param.inp.out_mul, 0);
    EXPECT_FALSE(param.inp.restart_save);
    EXPECT_FALSE(param.inp.restart_load);
    EXPECT_FALSE(param.inp.test_skip_ewald);
    EXPECT_EQ(param.inp.dft_plus_u, 0);
    EXPECT_FALSE(param.inp.yukawa_potential);
    EXPECT_DOUBLE_EQ(param.inp.yukawa_lambda, -1.0);
    EXPECT_EQ(param.inp.onsite_radius, 0.0);
    EXPECT_EQ(param.inp.omc, 0);
    EXPECT_FALSE(param.inp.dft_plus_dmft);
    EXPECT_FALSE(param.inp.rpa);
    EXPECT_EQ(param.inp.imp_sol, 0);
    EXPECT_DOUBLE_EQ(param.inp.eb_k, 80.0);
    EXPECT_DOUBLE_EQ(param.inp.tau, 1.0798 * 1e-5);
    EXPECT_DOUBLE_EQ(param.inp.sigma_k, 0.6);
    EXPECT_DOUBLE_EQ(param.inp.nc_k, 0.00037);
    EXPECT_EQ(param.inp.of_kinetic, "vw");
    EXPECT_EQ(param.inp.of_method, "tn");
    EXPECT_EQ(param.inp.of_conv, "energy");
    EXPECT_DOUBLE_EQ(param.inp.of_tole, 1e-6);
    EXPECT_DOUBLE_EQ(param.inp.of_tolp, 1e-5);
    EXPECT_DOUBLE_EQ(param.inp.of_tf_weight, 1.);
    EXPECT_DOUBLE_EQ(param.inp.of_vw_weight, 1.);
    EXPECT_DOUBLE_EQ(param.inp.of_wt_alpha, 0.833333);
    EXPECT_DOUBLE_EQ(param.inp.of_wt_beta, 0.833333);
    EXPECT_DOUBLE_EQ(param.inp.of_wt_rho0, 1.);
    EXPECT_TRUE(param.inp.of_hold_rho0);
    EXPECT_DOUBLE_EQ(param.inp.of_lkt_a, 1.3);
    EXPECT_FALSE(param.inp.of_full_pw);
    EXPECT_EQ(param.inp.of_full_pw_dim, 0);
    EXPECT_FALSE(param.inp.of_read_kernel);
    EXPECT_EQ(param.inp.of_kernel_file, "WTkernel.txt");
    EXPECT_EQ(param.inp.device, "cpu");
    EXPECT_NEAR(param.inp.force_thr_ev, 0.025711245953622324, 1e-8);
    EXPECT_DOUBLE_EQ(param.globalv.hubbard_u[0], 0);
    EXPECT_EQ(param.inp.orbital_corr[0], -1);
    EXPECT_EQ(param.inp.mdp.lj_rule, 2);
    EXPECT_FALSE(param.inp.mdp.lj_eshift);
    EXPECT_NEAR(param.inp.mdp.lj_epsilon[0], 0.01032, 1e-7);
    EXPECT_NEAR(param.inp.mdp.lj_rcut[0], 8.5, 1e-7);
    EXPECT_NEAR(param.inp.mdp.lj_sigma[0], 3.405, 1e-7);
    EXPECT_EQ(param.inp.mdp.md_damp, 1);
    EXPECT_EQ(param.inp.mdp.md_dt, 1);
    EXPECT_EQ(param.inp.mdp.md_dumpfreq, 1);
    EXPECT_EQ(param.inp.mdp.md_nraise, 1);
    EXPECT_EQ(param.inp.cal_syns, 0);
    EXPECT_EQ(param.inp.dmax, 0.01);
    EXPECT_EQ(param.inp.mdp.md_nstep, 10);
    EXPECT_EQ(param.inp.mdp.md_pchain, 1);
    EXPECT_EQ(param.inp.mdp.md_pcouple, "xyz");
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_pfirst, -1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_pfreq, 0);
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_plast, -1);
    EXPECT_EQ(param.inp.mdp.md_pmode, "iso");
    EXPECT_EQ(param.inp.mdp.md_restart, 0);
    EXPECT_EQ(param.inp.mdp.md_restartfreq, 5);
    EXPECT_EQ(param.inp.mdp.md_seed, -1);
    EXPECT_EQ(param.inp.mdp.md_prec_level, 0);
    EXPECT_DOUBLE_EQ(param.inp.ref_cell_factor, 1.2);
    EXPECT_EQ(param.inp.mdp.md_tchain, 1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_tfirst, -1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_tfreq, 0);
    EXPECT_EQ(param.inp.mdp.md_thermostat, "nhc");
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_tlast, -1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.md_tolerance, 100);
    EXPECT_EQ(param.inp.mdp.md_type, "nvt");
    EXPECT_EQ(param.inp.mdp.msst_direction, 2);
    EXPECT_DOUBLE_EQ(param.inp.mdp.msst_qmass, 1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.msst_tscale, 0.01);
    EXPECT_DOUBLE_EQ(param.inp.mdp.msst_vel, 0);
    EXPECT_DOUBLE_EQ(param.inp.mdp.msst_vis, 0);
    EXPECT_EQ(param.inp.mdp.pot_file, "graph.pb");
    EXPECT_DOUBLE_EQ(param.inp.mdp.dp_rescaling, 1.0);
    EXPECT_EQ(param.inp.mdp.dp_fparam.size(), 2);
    EXPECT_EQ(param.inp.mdp.dp_aparam.size(), 2);
    EXPECT_DOUBLE_EQ(param.inp.mdp.dp_fparam[0], 1.0);
    EXPECT_DOUBLE_EQ(param.inp.mdp.dp_fparam[1], 1.1);
    EXPECT_DOUBLE_EQ(param.inp.mdp.dp_aparam[0], 1.0);
    EXPECT_DOUBLE_EQ(param.inp.mdp.dp_aparam[1], 1.2);
    EXPECT_FALSE(param.inp.mdp.dump_force);
    EXPECT_FALSE(param.inp.mdp.dump_vel);
    EXPECT_FALSE(param.inp.mdp.dump_virial);
    EXPECT_EQ(param.inp.sc_mag_switch, 0);
    EXPECT_TRUE(param.inp.decay_grad_switch);
    EXPECT_DOUBLE_EQ(param.inp.sc_thr, 1e-4);
    EXPECT_EQ(param.inp.nsc, 50);
    EXPECT_EQ(param.inp.nsc_min, 4);
    EXPECT_EQ(param.inp.sc_scf_nmin, 4);
    EXPECT_DOUBLE_EQ(param.inp.alpha_trial, 0.02);
    EXPECT_DOUBLE_EQ(param.inp.sccut, 4.0);
    EXPECT_EQ(param.inp.sc_scf_thr, 1e-3);
    EXPECT_EQ(param.inp.sc_drop_thr, 1e-3);
    EXPECT_EQ(param.inp.lr_nstates, 1);
    EXPECT_EQ(param.inp.nocc, param.inp.nbands);
    EXPECT_EQ(param.inp.nvirt, 1);
    EXPECT_EQ(param.inp.xc_kernel, "LDA");
    EXPECT_EQ(param.inp.lr_init_xc_kernel[0], "default");
    EXPECT_EQ(param.inp.lr_solver, "dav");
    EXPECT_DOUBLE_EQ(param.inp.lr_thr, 1e-2);
    EXPECT_FALSE(param.inp.lr_unrestricted);
    EXPECT_FALSE(param.inp.out_wfc_lr);
    EXPECT_EQ(param.inp.abs_wavelen_range.size(), 2);
    EXPECT_DOUBLE_EQ(param.inp.abs_wavelen_range[0], 0.0);
    EXPECT_DOUBLE_EQ(param.inp.abs_broadening, 0.01);
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
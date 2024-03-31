#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include <stdlib.h>

/************************************************
 *  unit test of Input::bcast
 ***********************************************/

/**
 * - Tested Functions:
 *   - bcast()
 *     - bcast input parameters to all processes
 */

#define private public
#include "module_io/input.h"

class InputParaTest : public ::testing::Test
{
  protected:
};

#ifdef __MPI
TEST_F(InputParaTest, Bcast)
{
    INPUT.Default();
    if (GlobalV::MY_RANK == 0)
    {
        INPUT.suffix = "BcastTest";
    }
    INPUT.Bcast();
    EXPECT_EQ(INPUT.suffix, "BcastTest");
    EXPECT_EQ(INPUT.stru_file, "");
    EXPECT_EQ(INPUT.kpoint_file, "");
    EXPECT_EQ(INPUT.pseudo_dir, "");
    EXPECT_EQ(INPUT.orbital_dir, "");
    EXPECT_EQ(INPUT.read_file_dir, "auto");
    EXPECT_EQ(INPUT.wannier_card, "none");
    EXPECT_EQ(INPUT.latname, "none");
    EXPECT_EQ(INPUT.calculation, "scf");
    EXPECT_EQ(INPUT.esolver_type, "ksdft");
    EXPECT_DOUBLE_EQ(INPUT.pseudo_rcut, 15.0);
    EXPECT_FALSE(INPUT.pseudo_mesh);
    EXPECT_EQ(INPUT.ntype, 0);
    EXPECT_EQ(INPUT.nbands, 0);
    EXPECT_EQ(INPUT.nbands_sto, 256);
    EXPECT_EQ(INPUT.nbands_istate, 5);
    EXPECT_EQ(INPUT.pw_seed, 1);
    EXPECT_EQ(INPUT.emin_sto, 0.0);
    EXPECT_EQ(INPUT.emax_sto, 0.0);
    EXPECT_EQ(INPUT.nche_sto, 100);
    EXPECT_EQ(INPUT.seed_sto, 0);
    EXPECT_EQ(INPUT.initsto_ecut, 0.0);
    EXPECT_EQ(INPUT.bndpar, 1);
    EXPECT_EQ(INPUT.kpar, 1);
    EXPECT_EQ(INPUT.initsto_freq, 0);
    EXPECT_EQ(INPUT.method_sto, 2);
    EXPECT_EQ(INPUT.npart_sto, 1);
    EXPECT_FALSE(INPUT.cal_cond);
    EXPECT_EQ(INPUT.dos_nche, 100);
    EXPECT_DOUBLE_EQ(INPUT.cond_che_thr, 1e-8);
    EXPECT_DOUBLE_EQ(INPUT.cond_dw, 0.1);
    EXPECT_DOUBLE_EQ(INPUT.cond_wcut, 10);
    EXPECT_EQ(INPUT.cond_dt, 0.02);
    EXPECT_EQ(INPUT.cond_dtbatch, 0);
    EXPECT_EQ(INPUT.cond_smear, 1);
    EXPECT_DOUBLE_EQ(INPUT.cond_fwhm, 0.4);
    EXPECT_TRUE(INPUT.cond_nonlocal);
    EXPECT_FALSE(INPUT.berry_phase);
    EXPECT_EQ(INPUT.gdir, 3);
    EXPECT_FALSE(INPUT.towannier90);
    EXPECT_EQ(INPUT.nnkpfile, "seedname.nnkp");
    EXPECT_EQ(INPUT.wannier_spin, "up");
    EXPECT_EQ(INPUT.wannier_method, 1);
    EXPECT_TRUE(INPUT.out_wannier_mmn);
    EXPECT_TRUE(INPUT.out_wannier_amn);
    EXPECT_FALSE(INPUT.out_wannier_unk);
    EXPECT_TRUE(INPUT.out_wannier_eig);
    EXPECT_TRUE(INPUT.out_wannier_wvfn_formatted);
    EXPECT_DOUBLE_EQ(INPUT.kspacing[0], 0.0);
    EXPECT_DOUBLE_EQ(INPUT.kspacing[1], 0.0);
    EXPECT_DOUBLE_EQ(INPUT.kspacing[2], 0.0);
    EXPECT_DOUBLE_EQ(INPUT.min_dist_coef, 0.2);
    EXPECT_EQ(INPUT.dft_functional, "default");
    EXPECT_DOUBLE_EQ(INPUT.xc_temperature, 0.0);
    EXPECT_EQ(INPUT.nspin, 1);
    EXPECT_DOUBLE_EQ(INPUT.nelec, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.nelec_delta, 0.0);
    EXPECT_EQ(INPUT.lmaxmax, 2);
    EXPECT_EQ(INPUT.basis_type, "pw");
    EXPECT_EQ(INPUT.ks_solver, "default");
    EXPECT_DOUBLE_EQ(INPUT.search_radius, -1.0);
    EXPECT_TRUE(INPUT.search_pbc);
    EXPECT_EQ(INPUT.symmetry, "default");
    EXPECT_FALSE(INPUT.init_vel);
    EXPECT_DOUBLE_EQ(INPUT.symmetry_prec, 1.0e-6);
    EXPECT_TRUE(INPUT.symmetry_autoclose);
    EXPECT_EQ(INPUT.cal_force, 0);
    EXPECT_DOUBLE_EQ(INPUT.force_thr, 1.0e-3);
    EXPECT_DOUBLE_EQ(INPUT.force_thr_ev2, 0);
    EXPECT_DOUBLE_EQ(INPUT.stress_thr, 0.5);
    EXPECT_DOUBLE_EQ(INPUT.press1, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.press2, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.press3, 0.0);
    EXPECT_FALSE(INPUT.cal_stress);
    EXPECT_EQ(INPUT.fixed_axes, "None");
    EXPECT_FALSE(INPUT.fixed_ibrav);
    EXPECT_FALSE(INPUT.fixed_atoms);
    EXPECT_EQ(INPUT.relax_method, "cg");
    EXPECT_DOUBLE_EQ(INPUT.relax_cg_thr, 0.5);
    EXPECT_EQ(INPUT.out_level, "ie");
    EXPECT_FALSE(INPUT.out_md_control);
    EXPECT_TRUE(INPUT.relax_new);
    EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w1, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w2, 0.5);
    EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmax, 0.8);
    EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmin, 1e-5);
    EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_init, 0.5);
    EXPECT_DOUBLE_EQ(INPUT.relax_scale_force, 0.5);
    EXPECT_EQ(INPUT.nbspline, -1);
    EXPECT_FALSE(INPUT.gamma_only);
    EXPECT_FALSE(INPUT.gamma_only_local);
    EXPECT_DOUBLE_EQ(INPUT.ecutwfc, 50.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_ecut, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_height, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_sigma, 0.1);
    EXPECT_EQ(INPUT.fft_mode, 0);
    EXPECT_EQ(INPUT.nx, 0);
    EXPECT_EQ(INPUT.ny, 0);
    EXPECT_EQ(INPUT.nz, 0);
    EXPECT_EQ(INPUT.bx, 0);
    EXPECT_EQ(INPUT.by, 0);
    EXPECT_EQ(INPUT.bz, 0);
    EXPECT_EQ(INPUT.ndx, 0);
    EXPECT_EQ(INPUT.ndy, 0);
    EXPECT_EQ(INPUT.ndz, 0);
    EXPECT_EQ(INPUT.diago_proc, 0);
    EXPECT_EQ(INPUT.pw_diag_nmax, 50);
    EXPECT_EQ(INPUT.diago_cg_prec, 1);
    EXPECT_EQ(INPUT.pw_diag_ndim, 4);
    EXPECT_DOUBLE_EQ(INPUT.pw_diag_thr, 1.0e-2);
    EXPECT_EQ(INPUT.nb2d, 0);
    EXPECT_EQ(INPUT.nurse, 0);
    EXPECT_EQ(INPUT.colour, 0);
    EXPECT_EQ(INPUT.t_in_h, 1);
    EXPECT_EQ(INPUT.vl_in_h, 1);
    EXPECT_EQ(INPUT.vnl_in_h, 1);
    EXPECT_EQ(INPUT.vh_in_h, 1);
    EXPECT_EQ(INPUT.vion_in_h, 1);
    EXPECT_EQ(INPUT.test_force, 0);
    EXPECT_EQ(INPUT.test_stress, 0);
    EXPECT_DOUBLE_EQ(INPUT.scf_thr, -1.0);
    EXPECT_EQ(INPUT.scf_thr_type, -1);
    EXPECT_EQ(INPUT.scf_nmax, 100);
    EXPECT_EQ(INPUT.relax_nmax, 0);
    EXPECT_EQ(INPUT.out_stru, 0);
    //      EXPECT_EQ(INPUT.occupations,"smearing");
    EXPECT_EQ(INPUT.smearing_method, "gauss");
    EXPECT_DOUBLE_EQ(INPUT.smearing_sigma, 0.015);
    EXPECT_EQ(INPUT.mixing_mode, "broyden");
    EXPECT_DOUBLE_EQ(INPUT.mixing_beta, -10.0);
    EXPECT_EQ(INPUT.mixing_ndim, 8);
    EXPECT_DOUBLE_EQ(INPUT.mixing_gg0, 1.00);
    EXPECT_EQ(INPUT.init_wfc, "atomic");
    EXPECT_EQ(INPUT.mem_saver, 0);
    EXPECT_EQ(INPUT.printe, 100);
    EXPECT_EQ(INPUT.init_chg, "atomic");
    EXPECT_EQ(INPUT.chg_extrap, "default");
    EXPECT_EQ(INPUT.out_freq_elec, 0);
    EXPECT_EQ(INPUT.out_freq_ion, 0);
    EXPECT_EQ(INPUT.out_chg, 0);
    EXPECT_EQ(INPUT.out_dm, 0);
    EXPECT_EQ(INPUT.out_dm1, 0);
    EXPECT_EQ(INPUT.deepks_out_labels, 0);
    EXPECT_EQ(INPUT.deepks_scf, 0);
    EXPECT_EQ(INPUT.deepks_bandgap, 0);
    EXPECT_EQ(INPUT.deepks_out_unittest, 0);
    EXPECT_EQ(INPUT.out_pot, 0);
    EXPECT_EQ(INPUT.out_wfc_pw, 0);
    EXPECT_EQ(INPUT.out_wfc_r, 0);
    EXPECT_EQ(INPUT.out_dos, 0);
    EXPECT_EQ(INPUT.out_band[0], 0);
    EXPECT_EQ(INPUT.out_band[1], 8);
    EXPECT_EQ(INPUT.out_proj_band, 0);
    EXPECT_EQ(INPUT.out_mat_hs[0], 0);
    EXPECT_EQ(INPUT.out_mat_hs[1], 8);
    EXPECT_EQ(INPUT.out_mat_hs2, 0);
    EXPECT_EQ(INPUT.out_mat_xc, 0);
    EXPECT_EQ(INPUT.out_interval, 1);
    EXPECT_EQ(INPUT.out_app_flag, 1);
    EXPECT_EQ(INPUT.out_mat_r, 0);
    EXPECT_FALSE(INPUT.out_wfc_lcao);
    EXPECT_FALSE(INPUT.out_alllog);
    EXPECT_DOUBLE_EQ(INPUT.dos_emin_ev, -15);
    EXPECT_DOUBLE_EQ(INPUT.dos_emax_ev, 15);
    EXPECT_DOUBLE_EQ(INPUT.dos_edelta_ev, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.dos_scale, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.dos_sigma, 0.07);
    EXPECT_FALSE(INPUT.out_element_info);
    EXPECT_DOUBLE_EQ(INPUT.lcao_ecut, 0);
    EXPECT_DOUBLE_EQ(INPUT.lcao_dk, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.lcao_dr, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.lcao_rmax, 30);
    EXPECT_TRUE(INPUT.bessel_nao_smooth);
    EXPECT_DOUBLE_EQ(INPUT.bessel_nao_sigma, 0.1);
    EXPECT_EQ(INPUT.bessel_nao_ecut, "default");
    EXPECT_DOUBLE_EQ(INPUT.bessel_nao_rcut, 6.0);
    /* I need to test INPUT.bessel_nao_rcuts here */
    EXPECT_EQ(INPUT.bessel_nao_rcuts.size(), 0);
    EXPECT_DOUBLE_EQ(INPUT.bessel_nao_tolerence, 1E-12);
    EXPECT_EQ(INPUT.bessel_descriptor_lmax, 2);
    EXPECT_TRUE(INPUT.bessel_descriptor_smooth);
    EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_sigma, 0.1);
    EXPECT_EQ(INPUT.bessel_descriptor_ecut, "default");
    EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_rcut, 6.0);
    EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_tolerence, 1E-12);

    EXPECT_FALSE(INPUT.efield_flag);
    EXPECT_FALSE(INPUT.dip_cor_flag);
    EXPECT_EQ(INPUT.efield_dir, 2);
    EXPECT_DOUBLE_EQ(INPUT.efield_pos_max, -1.0);
    EXPECT_DOUBLE_EQ(INPUT.efield_pos_dec, -1.0);
    EXPECT_DOUBLE_EQ(INPUT.efield_amp, 0.0);
    EXPECT_FALSE(INPUT.gate_flag);
    EXPECT_DOUBLE_EQ(INPUT.zgate, 0.5);
    EXPECT_FALSE(INPUT.relax);
    EXPECT_FALSE(INPUT.block);
    EXPECT_DOUBLE_EQ(INPUT.block_down, 0.45);
    EXPECT_DOUBLE_EQ(INPUT.block_up, 0.55);
    EXPECT_DOUBLE_EQ(INPUT.block_height, 0.1);
    EXPECT_EQ(INPUT.vdw_method, "none");
    EXPECT_EQ(INPUT.vdw_s6, "default");
    EXPECT_EQ(INPUT.vdw_s8, "default");
    EXPECT_EQ(INPUT.vdw_a1, "default");
    EXPECT_EQ(INPUT.vdw_a2, "default");
    EXPECT_DOUBLE_EQ(INPUT.vdw_d, 20);
    EXPECT_FALSE(INPUT.vdw_abc);
    EXPECT_EQ(INPUT.vdw_cutoff_radius, "default");
    EXPECT_EQ(INPUT.vdw_radius_unit, "Bohr");
    EXPECT_DOUBLE_EQ(INPUT.vdw_cn_thr, 40.0);
    EXPECT_EQ(INPUT.vdw_cn_thr_unit, "Bohr");
    EXPECT_EQ(INPUT.vdw_C6_file, "default");
    EXPECT_EQ(INPUT.vdw_C6_unit, "Jnm6/mol");
    EXPECT_EQ(INPUT.vdw_R0_file, "default");
    EXPECT_EQ(INPUT.vdw_R0_unit, "A");
    EXPECT_EQ(INPUT.vdw_cutoff_type, "radius");
    EXPECT_EQ(INPUT.vdw_cutoff_period[0], 3);
    EXPECT_EQ(INPUT.vdw_cutoff_period[1], 3);
    EXPECT_EQ(INPUT.vdw_cutoff_period[2], 3);
    EXPECT_EQ(INPUT.exx_hybrid_alpha, "default");
    EXPECT_EQ(INPUT.exx_real_number, "default");
    EXPECT_DOUBLE_EQ(INPUT.exx_hse_omega, 0.11);
    EXPECT_TRUE(INPUT.exx_separate_loop);
    EXPECT_EQ(INPUT.exx_hybrid_step, 100);
    EXPECT_DOUBLE_EQ(INPUT.exx_lambda, 0.3);
    EXPECT_DOUBLE_EQ(INPUT.exx_mixing_beta, 1.0);
    EXPECT_DOUBLE_EQ(INPUT.exx_pca_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(INPUT.exx_c_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(INPUT.exx_v_threshold, 1E-1);
    EXPECT_DOUBLE_EQ(INPUT.exx_dm_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(INPUT.exx_schwarz_threshold, 0);
    EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_threshold, 1E-7);
    EXPECT_DOUBLE_EQ(INPUT.exx_c_grad_threshold, 1E-4);
    EXPECT_DOUBLE_EQ(INPUT.exx_v_grad_threshold, 1E-1);
    EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_force_threshold, 1E-7);
    EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_stress_threshold, 1E-7);
    EXPECT_DOUBLE_EQ(INPUT.exx_ccp_threshold, 1E-8);
    EXPECT_EQ(INPUT.exx_ccp_rmesh_times, "default");
    EXPECT_EQ(INPUT.exx_distribute_type, "htime");
    EXPECT_EQ(INPUT.exx_opt_orb_lmax, 0);
    EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_ecut, 0.0);
    EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_tolerence, 0.0);
    EXPECT_FALSE(INPUT.noncolin);
    EXPECT_FALSE(INPUT.lspinorb);
    EXPECT_DOUBLE_EQ(INPUT.soc_lambda, 1.0);
    EXPECT_EQ(INPUT.input_error, 0);
    EXPECT_DOUBLE_EQ(INPUT.td_force_dt, 0.02);
    EXPECT_FALSE(INPUT.td_vext);
    EXPECT_EQ(INPUT.td_vext_dire, "1");
    EXPECT_EQ(INPUT.propagator, 0);
    EXPECT_EQ(INPUT.td_stype, 0);
    EXPECT_EQ(INPUT.td_ttype, "0");
    EXPECT_EQ(INPUT.td_tstart, 1);
    EXPECT_EQ(INPUT.td_tend, 1000);
    EXPECT_EQ(INPUT.td_lcut1, 0.05);
    EXPECT_EQ(INPUT.td_lcut2, 0.95);
    EXPECT_EQ(INPUT.td_gauss_amp, "0.25");
    EXPECT_EQ(INPUT.td_gauss_freq, "22.13");
    EXPECT_EQ(INPUT.td_gauss_phase, "0.0");
    EXPECT_EQ(INPUT.td_gauss_t0, "100.0");
    EXPECT_EQ(INPUT.td_gauss_sigma, "30.0");
    EXPECT_EQ(INPUT.td_trape_amp, "2.74");
    EXPECT_EQ(INPUT.td_trape_freq, "1.60");
    EXPECT_EQ(INPUT.td_trape_phase, "0.0");
    EXPECT_EQ(INPUT.td_trape_t1, "1875");
    EXPECT_EQ(INPUT.td_trape_t2, "5625");
    EXPECT_EQ(INPUT.td_trape_t3, "7500");
    EXPECT_EQ(INPUT.td_trigo_freq1, "1.164656");
    EXPECT_EQ(INPUT.td_trigo_freq2, "0.029116");
    EXPECT_EQ(INPUT.td_trigo_phase1, "0.0");
    EXPECT_EQ(INPUT.td_trigo_phase2, "0.0");
    EXPECT_EQ(INPUT.td_trigo_amp, "2.74");
    EXPECT_EQ(INPUT.td_heavi_t0, "100");
    EXPECT_EQ(INPUT.td_heavi_amp, "1.0");
    EXPECT_EQ(INPUT.out_dipole, 0);
    EXPECT_EQ(INPUT.out_efield, 0);
    EXPECT_EQ(INPUT.td_print_eij, -1.0);
    EXPECT_EQ(INPUT.td_edm, 0);
    EXPECT_DOUBLE_EQ(INPUT.cell_factor, 1.2);
    EXPECT_EQ(INPUT.out_mul, 0);
    EXPECT_FALSE(INPUT.restart_save);
    EXPECT_FALSE(INPUT.restart_load);
    EXPECT_FALSE(INPUT.test_skip_ewald);
    EXPECT_EQ(INPUT.dft_plus_u, 0);
    EXPECT_FALSE(INPUT.yukawa_potential);
    EXPECT_DOUBLE_EQ(INPUT.yukawa_lambda, -1.0);
    EXPECT_EQ(INPUT.omc, 0);
    EXPECT_FALSE(INPUT.dft_plus_dmft);
    EXPECT_FALSE(INPUT.rpa);
    //        EXPECT_EQ(INPUT.coulomb_type,"full");
    EXPECT_EQ(INPUT.imp_sol, 0);
    EXPECT_DOUBLE_EQ(INPUT.eb_k, 80.0);
    EXPECT_DOUBLE_EQ(INPUT.tau, 1.0798 * 1e-5);
    EXPECT_DOUBLE_EQ(INPUT.sigma_k, 0.6);
    EXPECT_DOUBLE_EQ(INPUT.nc_k, 0.00037);
    EXPECT_EQ(INPUT.of_kinetic, "wt");
    EXPECT_EQ(INPUT.of_method, "tn");
    EXPECT_EQ(INPUT.of_conv, "energy");
    EXPECT_DOUBLE_EQ(INPUT.of_tole, 1e-6);
    EXPECT_DOUBLE_EQ(INPUT.of_tolp, 1e-5);
    EXPECT_DOUBLE_EQ(INPUT.of_tf_weight, 1.);
    EXPECT_DOUBLE_EQ(INPUT.of_vw_weight, 1.);
    EXPECT_DOUBLE_EQ(INPUT.of_wt_alpha, 5. / 6.);
    EXPECT_DOUBLE_EQ(INPUT.of_wt_beta, 5. / 6.);
    EXPECT_DOUBLE_EQ(INPUT.of_wt_rho0, 0.);
    EXPECT_FALSE(INPUT.of_hold_rho0);
    EXPECT_DOUBLE_EQ(INPUT.of_lkt_a, 1.3);
    EXPECT_TRUE(INPUT.of_full_pw);
    EXPECT_EQ(INPUT.of_full_pw_dim, 0);
    EXPECT_FALSE(INPUT.of_read_kernel);
    EXPECT_EQ(INPUT.of_kernel_file, "WTkernel.txt");
    EXPECT_EQ(INPUT.device, "cpu");
    EXPECT_DOUBLE_EQ(INPUT.ecutrho, 0.0);
    EXPECT_EQ(INPUT.ncx, 0);
    EXPECT_EQ(INPUT.ncy, 0);
    EXPECT_EQ(INPUT.ncz, 0);
    EXPECT_NEAR(INPUT.mdp.lj_epsilon, 0.01032, 1e-7);
    EXPECT_NEAR(INPUT.mdp.lj_rcut, 8.5, 1e-7);
    EXPECT_NEAR(INPUT.mdp.lj_sigma, 3.405, 1e-7);
    EXPECT_EQ(INPUT.mdp.md_damp, 1);
    EXPECT_EQ(INPUT.mdp.md_dt, 1);
    EXPECT_EQ(INPUT.mdp.md_dumpfreq, 1);
    EXPECT_EQ(INPUT.mdp.md_nraise, 1);
    EXPECT_EQ(INPUT.mdp.md_nstep, 10);
    EXPECT_EQ(INPUT.mdp.md_pchain, 1);
    EXPECT_EQ(INPUT.mdp.md_pcouple, "none");
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfirst, -1);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfreq, 0);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_plast, -1);
    EXPECT_EQ(INPUT.mdp.md_pmode, "iso");
    EXPECT_EQ(INPUT.mdp.md_restart, 0);
    EXPECT_EQ(INPUT.mdp.md_restartfreq, 5);
    EXPECT_EQ(INPUT.mdp.md_seed, -1);
    EXPECT_EQ(INPUT.mdp.md_tchain, 1);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfirst, -1);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfreq, 0);
    EXPECT_EQ(INPUT.mdp.md_thermostat, "nhc");
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_tlast, -1);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_tolerance, 100);
    EXPECT_EQ(INPUT.mdp.md_type, "nvt");
    EXPECT_EQ(INPUT.mdp.msst_direction, 2);
    EXPECT_DOUBLE_EQ(INPUT.mdp.msst_qmass, -1);
    EXPECT_DOUBLE_EQ(INPUT.mdp.msst_tscale, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vel, 0);
    EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vis, 0);
    EXPECT_EQ(INPUT.mdp.pot_file, "graph.pb");
    EXPECT_TRUE(INPUT.mdp.dump_force);
    EXPECT_TRUE(INPUT.mdp.dump_vel);
    EXPECT_TRUE(INPUT.mdp.dump_virial);
    EXPECT_FALSE(INPUT.mixing_tau);
    EXPECT_FALSE(INPUT.mixing_dftu);
    EXPECT_EQ(INPUT.mixing_restart,0.0);
    EXPECT_EQ(INPUT.mixing_dmr,false);
    EXPECT_EQ(INPUT.out_bandgap, 0);
    EXPECT_EQ(INPUT.out_mat_t, 0);

    /* I need to test qo_switch, qo_strategy, qo_screening_coeff, qo_thr and qo_basis */
    EXPECT_EQ(INPUT.qo_switch, 0);
    EXPECT_EQ(INPUT.qo_strategy.size(), 0);
    EXPECT_EQ(INPUT.qo_screening_coeff.size(), 0);
    EXPECT_EQ(INPUT.qo_thr, 1e-6);
    EXPECT_EQ(INPUT.qo_basis, "szv");

    EXPECT_EQ(INPUT.pexsi_npole, 40);
    EXPECT_TRUE(INPUT.pexsi_inertia);
    EXPECT_EQ(INPUT.pexsi_nmax, 80);
    EXPECT_TRUE(INPUT.pexsi_comm);
    EXPECT_TRUE(INPUT.pexsi_storage);
    EXPECT_EQ(INPUT.pexsi_ordering, 0);
    EXPECT_EQ(INPUT.pexsi_row_ordering, 1);
    EXPECT_EQ(INPUT.pexsi_nproc, 1);
    EXPECT_TRUE(INPUT.pexsi_symm);
    EXPECT_FALSE(INPUT.pexsi_trans);
    EXPECT_EQ(INPUT.pexsi_method, 1);
    EXPECT_EQ(INPUT.pexsi_nproc_pole, 1);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_temp, 0.015);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_gap, 0);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_delta_e, 20);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu_lower, -10);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu_upper, 10);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu, 0);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu_thr, 0.05);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu_expand, 0.3);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_mu_guard, 0.2);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_elec_thr, 0.001);
    EXPECT_DOUBLE_EQ(INPUT.pexsi_zero_thr, 1e-10);
}

TEST_F(InputParaTest, Init)
{
    std::string input_file = "./support/INPUT";
    Input input_tmp;
    EXPECT_NO_THROW(input_tmp.Init(input_file));
    if (GlobalV::MY_RANK == 0)
    {
        int status = system("rm -r ./OUT.autotest/");
        EXPECT_EQ(status, 0);
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
#endif
#undef private

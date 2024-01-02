#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
/************************************************
 *  unit test of input.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Default()
 *     - set default value of input parameters
 *   - Read()
 *     - read input parameters from input files
 */

#define private public
#include "module_io/input.h"

class InputTest : public testing::Test
{
protected:
	std::string output;
};

TEST_F(InputTest, Default)
{
	INPUT.Default();
	EXPECT_EQ(INPUT.suffix,"ABACUS");
	EXPECT_EQ(INPUT.stru_file,"");
	EXPECT_EQ(INPUT.kpoint_file,"");
	EXPECT_EQ(INPUT.pseudo_dir,"");
	EXPECT_EQ(INPUT.orbital_dir,"");
	EXPECT_EQ(INPUT.read_file_dir,"auto");
	EXPECT_EQ(INPUT.wannier_card,"none");
	EXPECT_EQ(INPUT.latname,"none");
	EXPECT_EQ(INPUT.calculation,"scf");
	EXPECT_EQ(INPUT.esolver_type,"ksdft");
	EXPECT_DOUBLE_EQ(INPUT.pseudo_rcut,15.0);
	EXPECT_FALSE(INPUT.pseudo_mesh);
	EXPECT_EQ(INPUT.ntype,0);
	EXPECT_EQ(INPUT.nbands,0);
	EXPECT_EQ(INPUT.nbands_sto,256);
	EXPECT_EQ(INPUT.nbands_istate,5);
	EXPECT_EQ(INPUT.pw_seed,1);
	EXPECT_EQ(INPUT.emin_sto,0.0);
	EXPECT_EQ(INPUT.emax_sto,0.0);
	EXPECT_EQ(INPUT.nche_sto,100);
        EXPECT_EQ(INPUT.seed_sto,0);
		EXPECT_EQ(INPUT.initsto_ecut,0.0);
        EXPECT_EQ(INPUT.bndpar,1);
        EXPECT_EQ(INPUT.kpar,1);
        EXPECT_EQ(INPUT.initsto_freq,0);
        EXPECT_EQ(INPUT.method_sto,2);
        EXPECT_EQ(INPUT.npart_sto,1);
        EXPECT_FALSE(INPUT.cal_cond);
        EXPECT_EQ(INPUT.dos_nche,100);
        EXPECT_DOUBLE_EQ(INPUT.cond_che_thr,1e-8);
        EXPECT_DOUBLE_EQ(INPUT.cond_dw,0.1);
        EXPECT_DOUBLE_EQ(INPUT.cond_wcut,10);
        EXPECT_EQ(INPUT.cond_dt,0.02);
		EXPECT_EQ(INPUT.cond_dtbatch,0);
		EXPECT_EQ(INPUT.cond_smear,1);
        EXPECT_DOUBLE_EQ(INPUT.cond_fwhm,0.4);
        EXPECT_TRUE(INPUT.cond_nonlocal);
        EXPECT_FALSE(INPUT.berry_phase);
        EXPECT_EQ(INPUT.gdir,3);
        EXPECT_FALSE(INPUT.towannier90);
        EXPECT_EQ(INPUT.nnkpfile,"seedname.nnkp");
        EXPECT_EQ(INPUT.wannier_spin,"up");
        EXPECT_EQ(INPUT.wannier_method,1);
		EXPECT_TRUE(INPUT.out_wannier_amn);
		EXPECT_TRUE(INPUT.out_wannier_mmn);
		EXPECT_FALSE(INPUT.out_wannier_unk);
		EXPECT_TRUE(INPUT.out_wannier_eig);
        EXPECT_TRUE(INPUT.out_wannier_wvfn_formatted);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[0], 0.0);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[1],0.0);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[2],0.0);
        EXPECT_DOUBLE_EQ(INPUT.min_dist_coef,0.2);
        EXPECT_EQ(INPUT.dft_functional,"default");
        EXPECT_DOUBLE_EQ(INPUT.xc_temperature,0.0);
        EXPECT_EQ(INPUT.nspin,1);
        EXPECT_DOUBLE_EQ(INPUT.nelec,0.0);
        EXPECT_EQ(INPUT.lmaxmax,2);
        EXPECT_EQ(INPUT.basis_type,"pw");
        EXPECT_EQ(INPUT.ks_solver,"default");
        EXPECT_DOUBLE_EQ(INPUT.search_radius,-1.0);
        EXPECT_TRUE(INPUT.search_pbc);
        EXPECT_EQ(INPUT.symmetry,"default");
        EXPECT_FALSE(INPUT.init_vel);
        EXPECT_DOUBLE_EQ(INPUT.ref_cell_factor,1.0);
        EXPECT_DOUBLE_EQ(INPUT.symmetry_prec, 1.0e-6);
        EXPECT_TRUE(INPUT.symmetry_autoclose);
        EXPECT_EQ(INPUT.cal_force, 0);
        EXPECT_DOUBLE_EQ(INPUT.force_thr,1.0e-3);
        EXPECT_DOUBLE_EQ(INPUT.force_thr_ev2,0);
        EXPECT_DOUBLE_EQ(INPUT.stress_thr, 0.5);
        EXPECT_DOUBLE_EQ(INPUT.press1,0.0);
        EXPECT_DOUBLE_EQ(INPUT.press2,0.0);
        EXPECT_DOUBLE_EQ(INPUT.press3,0.0);
        EXPECT_FALSE(INPUT.cal_stress);
        EXPECT_EQ(INPUT.fixed_axes,"None");
        EXPECT_FALSE(INPUT.fixed_ibrav);
        EXPECT_FALSE(INPUT.fixed_atoms);
        EXPECT_EQ(INPUT.relax_method,"cg");
        EXPECT_DOUBLE_EQ(INPUT.relax_cg_thr,0.5);
        EXPECT_EQ(INPUT.out_level,"ie");
        EXPECT_FALSE(INPUT.out_md_control);
        EXPECT_TRUE(INPUT.relax_new);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w1,0.01);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w2,0.5);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmax,0.8);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmin,1e-5);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_init,0.5);
        EXPECT_DOUBLE_EQ(INPUT.relax_scale_force,0.5);
        EXPECT_EQ(INPUT.nbspline,-1);
        EXPECT_FALSE(INPUT.gamma_only);
        EXPECT_FALSE(INPUT.gamma_only_local);
        EXPECT_DOUBLE_EQ(INPUT.ecutwfc,50.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_ecut, 0.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_height, 0.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_sigma, 0.1);
		EXPECT_EQ(INPUT.fft_mode,0);
        EXPECT_EQ(INPUT.nx,0);
        EXPECT_EQ(INPUT.ny,0);
        EXPECT_EQ(INPUT.nz,0);
        EXPECT_EQ(INPUT.bx,0);
        EXPECT_EQ(INPUT.by,0);
        EXPECT_EQ(INPUT.bz,0);
        EXPECT_EQ(INPUT.ndx, 0);
        EXPECT_EQ(INPUT.ndy, 0);
        EXPECT_EQ(INPUT.ndz, 0);
        EXPECT_EQ(INPUT.diago_proc,0);
        EXPECT_EQ(INPUT.pw_diag_nmax,50);
        EXPECT_EQ(INPUT.diago_cg_prec,1);
        EXPECT_EQ(INPUT.pw_diag_ndim,4);
        EXPECT_DOUBLE_EQ(INPUT.pw_diag_thr,1.0e-2);
        EXPECT_EQ(INPUT.nb2d,0);
        EXPECT_EQ(INPUT.nurse,0);
        EXPECT_EQ(INPUT.colour,0);
        EXPECT_EQ(INPUT.t_in_h,1);
        EXPECT_EQ(INPUT.vl_in_h,1);
        EXPECT_EQ(INPUT.vnl_in_h,1);
        EXPECT_EQ(INPUT.vh_in_h,1);
        EXPECT_EQ(INPUT.vion_in_h,1);
        EXPECT_EQ(INPUT.test_force,0);
        EXPECT_EQ(INPUT.test_stress,0);
        EXPECT_DOUBLE_EQ(INPUT.scf_thr,-1.0);
        EXPECT_EQ(INPUT.scf_thr_type,-1);
        EXPECT_EQ(INPUT.scf_nmax,100);
        EXPECT_EQ(INPUT.relax_nmax,0);
        EXPECT_EQ(INPUT.out_stru,0);
        EXPECT_EQ(INPUT.occupations,"smearing");
        EXPECT_EQ(INPUT.smearing_method,"gauss");
        EXPECT_DOUBLE_EQ(INPUT.smearing_sigma,0.015);
        EXPECT_EQ(INPUT.mixing_mode,"broyden");
        EXPECT_DOUBLE_EQ(INPUT.mixing_beta,-10.0);
        EXPECT_EQ(INPUT.mixing_ndim,8);
        EXPECT_DOUBLE_EQ(INPUT.mixing_gg0,1.00);
		EXPECT_DOUBLE_EQ(INPUT.mixing_gg0_mag,0.00);
		EXPECT_DOUBLE_EQ(INPUT.mixing_beta_mag,-10.0);
		EXPECT_DOUBLE_EQ(INPUT.mixing_angle,-10.0);
        EXPECT_EQ(INPUT.init_wfc,"atomic");
        EXPECT_EQ(INPUT.mem_saver,0);
        EXPECT_EQ(INPUT.printe,100);
        EXPECT_EQ(INPUT.init_chg,"atomic");
        EXPECT_EQ(INPUT.chg_extrap, "default");
        EXPECT_EQ(INPUT.out_freq_elec,0);
        EXPECT_EQ(INPUT.out_freq_ion,0);
        EXPECT_EQ(INPUT.out_chg,0);
        EXPECT_EQ(INPUT.out_dm,0);
        EXPECT_EQ(INPUT.out_dm1,0);
        EXPECT_EQ(INPUT.deepks_out_labels,0);
        EXPECT_EQ(INPUT.deepks_scf,0);
        EXPECT_EQ(INPUT.deepks_bandgap,0);
        EXPECT_EQ(INPUT.deepks_out_unittest,0);
        EXPECT_EQ(INPUT.out_pot,0);
        EXPECT_EQ(INPUT.out_wfc_pw,0);
        EXPECT_EQ(INPUT.out_wfc_r,0);
        EXPECT_EQ(INPUT.out_dos,0);
        EXPECT_EQ(INPUT.out_band,0);
        EXPECT_EQ(INPUT.out_proj_band,0);
        EXPECT_EQ(INPUT.out_mat_hs,0);
        EXPECT_EQ(INPUT.out_mat_hs2,0);
        EXPECT_EQ(INPUT.out_interval,1);
        EXPECT_EQ(INPUT.out_app_flag,1);
        EXPECT_EQ(INPUT.out_mat_r,0);
        EXPECT_EQ(INPUT.out_wfc_lcao,0);
        EXPECT_FALSE(INPUT.out_alllog);
        EXPECT_DOUBLE_EQ(INPUT.dos_emin_ev,-15);
        EXPECT_DOUBLE_EQ(INPUT.dos_emax_ev,15);
        EXPECT_DOUBLE_EQ(INPUT.dos_edelta_ev,0.01);
        EXPECT_DOUBLE_EQ(INPUT.dos_scale,0.01);
        EXPECT_DOUBLE_EQ(INPUT.dos_sigma,0.07);
        EXPECT_FALSE(INPUT.out_element_info);
        EXPECT_DOUBLE_EQ(INPUT.lcao_ecut,0);
        EXPECT_DOUBLE_EQ(INPUT.lcao_dk,0.01);
        EXPECT_DOUBLE_EQ(INPUT.lcao_dr,0.01);
        EXPECT_DOUBLE_EQ(INPUT.lcao_rmax,30);
		EXPECT_TRUE(INPUT.bessel_nao_smooth);
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_sigma, 0.1);
		EXPECT_EQ(INPUT.bessel_nao_ecut, "default");
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_rcut, 6.0);
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_tolerence, 1E-12);
		EXPECT_EQ(INPUT.bessel_descriptor_lmax, 2);
		EXPECT_TRUE(INPUT.bessel_descriptor_smooth);
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_sigma, 0.1);
		EXPECT_EQ(INPUT.bessel_descriptor_ecut, "default");
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_rcut, 6.0);
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_tolerence, 1E-12);

        EXPECT_FALSE(INPUT.efield_flag);
        EXPECT_FALSE(INPUT.dip_cor_flag);
        EXPECT_EQ(INPUT.efield_dir,2);
        EXPECT_DOUBLE_EQ(INPUT.efield_pos_max,0.5);
        EXPECT_DOUBLE_EQ(INPUT.efield_pos_dec,0.1);
        EXPECT_DOUBLE_EQ(INPUT.efield_amp ,0.0);
        EXPECT_FALSE(INPUT.gate_flag);
        EXPECT_DOUBLE_EQ(INPUT.zgate,0.5);
        EXPECT_FALSE(INPUT.relax);
        EXPECT_FALSE(INPUT.block);
        EXPECT_DOUBLE_EQ(INPUT.block_down,0.45);
        EXPECT_DOUBLE_EQ(INPUT.block_up,0.55);
        EXPECT_DOUBLE_EQ(INPUT.block_height,0.1);
        EXPECT_EQ(INPUT.vdw_method,"none");
        EXPECT_EQ(INPUT.vdw_s6,"default");
        EXPECT_EQ(INPUT.vdw_s8,"default");
        EXPECT_EQ(INPUT.vdw_a1,"default");
        EXPECT_EQ(INPUT.vdw_a2,"default");
        EXPECT_DOUBLE_EQ(INPUT.vdw_d,20);
        EXPECT_FALSE(INPUT.vdw_abc);
        EXPECT_EQ(INPUT.vdw_cutoff_radius,"default");
        EXPECT_EQ(INPUT.vdw_radius_unit,"Bohr");
        EXPECT_DOUBLE_EQ(INPUT.vdw_cn_thr,40.0);
        EXPECT_EQ(INPUT.vdw_cn_thr_unit,"Bohr");
        EXPECT_EQ(INPUT.vdw_C6_file,"default");
        EXPECT_EQ(INPUT.vdw_C6_unit,"Jnm6/mol");
        EXPECT_EQ(INPUT.vdw_R0_file,"default");
        EXPECT_EQ(INPUT.vdw_R0_unit,"A");
        EXPECT_EQ(INPUT.vdw_cutoff_type,"radius");
        EXPECT_EQ(INPUT.vdw_cutoff_period[0],3);
        EXPECT_EQ(INPUT.vdw_cutoff_period[1],3);
        EXPECT_EQ(INPUT.vdw_cutoff_period[2],3);
        EXPECT_EQ(INPUT.exx_hybrid_alpha,"default");
        EXPECT_EQ(INPUT.exx_real_number,"default");
        EXPECT_DOUBLE_EQ(INPUT.exx_hse_omega,0.11);
        EXPECT_TRUE(INPUT.exx_separate_loop);
        EXPECT_EQ(INPUT.exx_hybrid_step,100);
        EXPECT_DOUBLE_EQ(INPUT.exx_lambda,0.3);
		EXPECT_DOUBLE_EQ(INPUT.exx_mixing_beta,1.0);
        EXPECT_DOUBLE_EQ(INPUT.exx_pca_threshold,1E-4);
        EXPECT_DOUBLE_EQ(INPUT.exx_c_threshold,1E-4);
        EXPECT_DOUBLE_EQ(INPUT.exx_v_threshold,1E-1);
        EXPECT_DOUBLE_EQ(INPUT.exx_dm_threshold,1E-4);
        EXPECT_DOUBLE_EQ(INPUT.exx_schwarz_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_threshold,1E-7);
        EXPECT_DOUBLE_EQ(INPUT.exx_c_grad_threshold,1E-4);
        EXPECT_DOUBLE_EQ(INPUT.exx_v_grad_threshold,1E-1);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_force_threshold,1E-7);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_stress_threshold,1E-7);
        EXPECT_DOUBLE_EQ(INPUT.exx_ccp_threshold,1E-8);
        EXPECT_EQ(INPUT.exx_ccp_rmesh_times,"default");
        EXPECT_EQ(INPUT.exx_distribute_type,"htime");
        EXPECT_EQ(INPUT.exx_opt_orb_lmax,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_ecut,0.0);
        EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_tolerence,0.0);
        EXPECT_FALSE(INPUT.noncolin);
        EXPECT_FALSE(INPUT.lspinorb);
        EXPECT_DOUBLE_EQ(INPUT.soc_lambda,1.0);
        EXPECT_EQ(INPUT.input_error,0);
        EXPECT_DOUBLE_EQ(INPUT.td_force_dt,0.02);
        EXPECT_FALSE(INPUT.td_vext);
        EXPECT_EQ(INPUT.td_vext_dire,"1");
		EXPECT_EQ(INPUT.propagator,0);
		EXPECT_EQ(INPUT.td_stype,0);
		EXPECT_EQ(INPUT.td_ttype,"0");
		EXPECT_EQ(INPUT.td_tstart,1);
		EXPECT_EQ(INPUT.td_tend,1000);
		EXPECT_EQ(INPUT.td_lcut1,0.05);
		EXPECT_EQ(INPUT.td_lcut2,0.95);
		EXPECT_EQ(INPUT.td_gauss_amp,"0.25");
		EXPECT_EQ(INPUT.td_gauss_freq,"22.13");
		EXPECT_EQ(INPUT.td_gauss_phase,"0.0");
		EXPECT_EQ(INPUT.td_gauss_t0,"100.0");
		EXPECT_EQ(INPUT.td_gauss_sigma,"30.0");
		EXPECT_EQ(INPUT.td_trape_amp,"2.74");
		EXPECT_EQ(INPUT.td_trape_freq,"1.60");
		EXPECT_EQ(INPUT.td_trape_phase,"0.0");
		EXPECT_EQ(INPUT.td_trape_t1,"1875");
		EXPECT_EQ(INPUT.td_trape_t2,"5625");
		EXPECT_EQ(INPUT.td_trape_t3,"7500");
		EXPECT_EQ(INPUT.td_trigo_freq1,"1.164656");
		EXPECT_EQ(INPUT.td_trigo_freq2,"0.029116");
		EXPECT_EQ(INPUT.td_trigo_phase1,"0.0");
		EXPECT_EQ(INPUT.td_trigo_phase2,"0.0");
		EXPECT_EQ(INPUT.td_trigo_amp,"2.74");
		EXPECT_EQ(INPUT.td_heavi_t0,"100");	
		EXPECT_EQ(INPUT.td_heavi_amp,"1.0");					
        EXPECT_EQ(INPUT.out_dipole,0);
		EXPECT_EQ(INPUT.out_efield,0);
		EXPECT_EQ(INPUT.td_print_eij,-1.0);
		EXPECT_EQ(INPUT.td_edm,0);
        EXPECT_DOUBLE_EQ(INPUT.cell_factor,1.2);
        EXPECT_EQ(INPUT.out_mul,0);
        EXPECT_FALSE(INPUT.restart_save);
        EXPECT_FALSE(INPUT.restart_load);
        EXPECT_FALSE(INPUT.test_skip_ewald);
        EXPECT_FALSE(INPUT.dft_plus_u);
        EXPECT_FALSE(INPUT.yukawa_potential);
        EXPECT_DOUBLE_EQ(INPUT.yukawa_lambda,-1.0);
        EXPECT_EQ(INPUT.omc,0);
        EXPECT_FALSE(INPUT.dft_plus_dmft);
        EXPECT_FALSE(INPUT.rpa);
        EXPECT_EQ(INPUT.coulomb_type,"full");
        EXPECT_EQ(INPUT.imp_sol,0);
        EXPECT_DOUBLE_EQ(INPUT.eb_k,80.0);
        EXPECT_DOUBLE_EQ(INPUT.tau,1.0798 * 1e-5);
        EXPECT_DOUBLE_EQ(INPUT.sigma_k,0.6);
        EXPECT_DOUBLE_EQ(INPUT.nc_k,0.00037);
        EXPECT_EQ(INPUT.of_kinetic,"wt");
        EXPECT_EQ(INPUT.of_method,"tn");
        EXPECT_EQ(INPUT.of_conv,"energy");
        EXPECT_DOUBLE_EQ(INPUT.of_tole,1e-6);
        EXPECT_DOUBLE_EQ(INPUT.of_tolp,1e-5);
        EXPECT_DOUBLE_EQ(INPUT.of_tf_weight,1.);
        EXPECT_DOUBLE_EQ(INPUT.of_vw_weight,1.);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_alpha,5./6.);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_beta,5./6.);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_rho0,0.);
        EXPECT_FALSE(INPUT.of_hold_rho0);
        EXPECT_DOUBLE_EQ(INPUT.of_lkt_a,1.3);
        EXPECT_TRUE(INPUT.of_full_pw);
        EXPECT_EQ(INPUT.of_full_pw_dim,0);
        EXPECT_FALSE(INPUT.of_read_kernel);
        EXPECT_EQ(INPUT.of_kernel_file,"WTkernel.txt");
        EXPECT_EQ(INPUT.device,"cpu");
        EXPECT_DOUBLE_EQ(INPUT.ecutrho,0.0);
        EXPECT_EQ(INPUT.ncx,0);
        EXPECT_EQ(INPUT.ncy,0);
        EXPECT_EQ(INPUT.ncz,0);
	EXPECT_NEAR(INPUT.mdp.lj_epsilon,0.01032,1e-7);
	EXPECT_NEAR(INPUT.mdp.lj_rcut,8.5,1e-7);
	EXPECT_NEAR(INPUT.mdp.lj_sigma,3.405,1e-7);
	EXPECT_EQ(INPUT.mdp.md_damp,1);
	EXPECT_EQ(INPUT.mdp.md_dt,1);
	EXPECT_EQ(INPUT.mdp.md_dumpfreq,1);
	EXPECT_EQ(INPUT.mdp.md_nraise,1);
	EXPECT_EQ(INPUT.cal_syns,0);
	EXPECT_EQ(INPUT.dmax,0.01);
	EXPECT_EQ(INPUT.mdp.md_nstep,10);
	EXPECT_EQ(INPUT.mdp.md_pchain,1);
	EXPECT_EQ(INPUT.mdp.md_pcouple,"none");
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfirst,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfreq,0);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_plast,-1);
	EXPECT_EQ(INPUT.mdp.md_pmode,"iso");
	EXPECT_EQ(INPUT.mdp.md_restart,0);
	EXPECT_EQ(INPUT.mdp.md_restartfreq,5);
	EXPECT_EQ(INPUT.mdp.md_seed,-1);
    EXPECT_EQ(INPUT.mdp.md_prec_level,0);
	EXPECT_EQ(INPUT.mdp.md_tchain,1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfirst,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfreq,0);
	EXPECT_EQ(INPUT.mdp.md_thermostat,"nhc");
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tlast,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tolerance,100);
	EXPECT_EQ(INPUT.mdp.md_type,"nvt");
	EXPECT_EQ(INPUT.mdp.msst_direction,2);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_qmass,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_tscale,0.01);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vel,0);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vis,0);
	EXPECT_EQ(INPUT.mdp.pot_file,"graph.pb");
    EXPECT_TRUE(INPUT.mdp.dump_force);
    EXPECT_TRUE(INPUT.mdp.dump_vel);
    EXPECT_TRUE(INPUT.mdp.dump_virial);
    EXPECT_EQ(INPUT.sc_mag_switch,0);
    EXPECT_FALSE(INPUT.decay_grad_switch);
    EXPECT_DOUBLE_EQ(INPUT.sc_thr, 1e-6);
    EXPECT_EQ(INPUT.nsc, 100);
    EXPECT_EQ(INPUT.nsc_min, 2);
	EXPECT_EQ(INPUT.sc_scf_nmin, 2);
    EXPECT_DOUBLE_EQ(INPUT.alpha_trial, 0.01);
    EXPECT_DOUBLE_EQ(INPUT.sccut, 3.0);
    EXPECT_EQ(INPUT.sc_file, "none");
}

TEST_F(InputTest, Read)
{
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	EXPECT_EQ(INPUT.suffix,"autotest");
	EXPECT_EQ(INPUT.stru_file,"./support/STRU");
	EXPECT_EQ(INPUT.kpoint_file,"KPT");
	EXPECT_EQ(INPUT.pseudo_dir,"../../PP_ORB/");
	EXPECT_EQ(INPUT.orbital_dir,"../../PP_ORB/");
	EXPECT_EQ(INPUT.read_file_dir,"auto");
	EXPECT_EQ(INPUT.wannier_card,"none");
	EXPECT_EQ(INPUT.latname,"none");
	EXPECT_EQ(INPUT.calculation,"scf");
	EXPECT_EQ(INPUT.esolver_type,"ksdft");
	EXPECT_DOUBLE_EQ(INPUT.pseudo_rcut,15.0);
	EXPECT_FALSE(INPUT.pseudo_mesh);
	EXPECT_EQ(INPUT.ntype,1);
	EXPECT_EQ(INPUT.nbands,8);
	EXPECT_EQ(INPUT.nbands_sto,256);
	EXPECT_EQ(INPUT.nbands_istate,5);
	EXPECT_EQ(INPUT.pw_seed,1);
	EXPECT_EQ(INPUT.emin_sto,0.0);
	EXPECT_EQ(INPUT.emax_sto,0.0);
	EXPECT_EQ(INPUT.nche_sto,100);
        EXPECT_EQ(INPUT.seed_sto,0);
		EXPECT_EQ(INPUT.initsto_ecut,0.0);
        EXPECT_EQ(INPUT.bndpar,1);
        EXPECT_EQ(INPUT.kpar,1);
        EXPECT_EQ(INPUT.initsto_freq,0);
        EXPECT_EQ(INPUT.method_sto,3);
        EXPECT_EQ(INPUT.npart_sto,1);
        EXPECT_FALSE(INPUT.cal_cond);
        EXPECT_EQ(INPUT.dos_nche,100);
        EXPECT_DOUBLE_EQ(INPUT.cond_che_thr,1e-8);
        EXPECT_DOUBLE_EQ(INPUT.cond_dw,0.1);
        EXPECT_DOUBLE_EQ(INPUT.cond_wcut,10);
        EXPECT_EQ(INPUT.cond_dt,0.07);
		EXPECT_EQ(INPUT.cond_dtbatch,2);
        EXPECT_DOUBLE_EQ(INPUT.cond_fwhm,0.3);
        EXPECT_TRUE(INPUT.cond_nonlocal);
        EXPECT_FALSE(INPUT.berry_phase);
        EXPECT_EQ(INPUT.gdir,3);
        EXPECT_FALSE(INPUT.towannier90);
        EXPECT_EQ(INPUT.nnkpfile,"seedname.nnkp");
        EXPECT_EQ(INPUT.wannier_spin,"up");
        EXPECT_EQ(INPUT.wannier_method,1);
		EXPECT_TRUE(INPUT.out_wannier_amn);
		EXPECT_TRUE(INPUT.out_wannier_mmn);
		EXPECT_TRUE(INPUT.out_wannier_unk);
		EXPECT_TRUE(INPUT.out_wannier_eig);
        EXPECT_TRUE(INPUT.out_wannier_wvfn_formatted);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[0], 0.0);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[1],0.0);
        EXPECT_DOUBLE_EQ(INPUT.kspacing[2],0.0);
        EXPECT_DOUBLE_EQ(INPUT.min_dist_coef,0.2);
        EXPECT_EQ(INPUT.dft_functional,"hse");
        EXPECT_DOUBLE_EQ(INPUT.xc_temperature,0.0);
        EXPECT_EQ(INPUT.nspin,1);
        EXPECT_DOUBLE_EQ(INPUT.nelec,0.0);
        EXPECT_EQ(INPUT.lmaxmax,2);
        EXPECT_EQ(INPUT.basis_type,"lcao");
        EXPECT_EQ(INPUT.ks_solver,"genelpa");
        EXPECT_DOUBLE_EQ(INPUT.search_radius,-1.0);
        EXPECT_TRUE(INPUT.search_pbc);
        EXPECT_EQ(INPUT.symmetry,"1");
        EXPECT_FALSE(INPUT.init_vel);
        EXPECT_DOUBLE_EQ(INPUT.symmetry_prec, 1.0e-6);
        EXPECT_TRUE(INPUT.symmetry_autoclose);
        EXPECT_EQ(INPUT.cal_force, 0);
        EXPECT_NEAR(INPUT.force_thr,1.0e-3,1.0e-7);
        EXPECT_DOUBLE_EQ(INPUT.force_thr_ev2,0);
        EXPECT_DOUBLE_EQ(INPUT.stress_thr,1.0e-2);
        EXPECT_DOUBLE_EQ(INPUT.press1,0.0);
        EXPECT_DOUBLE_EQ(INPUT.press2,0.0);
        EXPECT_DOUBLE_EQ(INPUT.press3,0.0);
        EXPECT_FALSE(INPUT.cal_stress);
        EXPECT_EQ(INPUT.fixed_axes,"None");
        EXPECT_FALSE(INPUT.fixed_ibrav);
        EXPECT_FALSE(INPUT.fixed_atoms);
        EXPECT_EQ(INPUT.relax_method,"cg");
        EXPECT_DOUBLE_EQ(INPUT.relax_cg_thr,0.5);
        EXPECT_EQ(INPUT.out_level,"ie");
        EXPECT_TRUE(INPUT.out_md_control);
        EXPECT_TRUE(INPUT.relax_new);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w1,0.01);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_w2,0.5);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmax,0.8);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_rmin,1e-5);
        EXPECT_DOUBLE_EQ(INPUT.relax_bfgs_init,0.5);
        EXPECT_DOUBLE_EQ(INPUT.relax_scale_force,0.5);
        EXPECT_EQ(INPUT.nbspline,-1);
        EXPECT_TRUE(INPUT.gamma_only);
        EXPECT_TRUE(INPUT.gamma_only_local);
        EXPECT_DOUBLE_EQ(INPUT.ecutwfc,20.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_ecut, 20.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_height, 20.0);
        EXPECT_DOUBLE_EQ(INPUT.erf_sigma, 4.0);
        EXPECT_DOUBLE_EQ(INPUT.ecutrho, 0.0);
		EXPECT_EQ(INPUT.fft_mode,0);
        EXPECT_EQ(INPUT.ncx,0);
        EXPECT_EQ(INPUT.ncy,0);
        EXPECT_EQ(INPUT.ncz,0);
        EXPECT_EQ(INPUT.nx,0);
        EXPECT_EQ(INPUT.ny,0);
        EXPECT_EQ(INPUT.nz,0);
        EXPECT_EQ(INPUT.bx,2);
        EXPECT_EQ(INPUT.by,2);
        EXPECT_EQ(INPUT.bz,2);
        EXPECT_EQ(INPUT.ndx, 0);
        EXPECT_EQ(INPUT.ndy, 0);
        EXPECT_EQ(INPUT.ndz, 0);
        EXPECT_EQ(INPUT.diago_proc,4);
        EXPECT_EQ(INPUT.pw_diag_nmax,50);
        EXPECT_EQ(INPUT.diago_cg_prec,1);
        EXPECT_EQ(INPUT.pw_diag_ndim,4);
        EXPECT_DOUBLE_EQ(INPUT.pw_diag_thr,1.0e-2);
        EXPECT_EQ(INPUT.nb2d,0);
        EXPECT_EQ(INPUT.nurse,0);
        EXPECT_EQ(INPUT.colour,0);
        EXPECT_EQ(INPUT.t_in_h,1);
        EXPECT_EQ(INPUT.vl_in_h,1);
        EXPECT_EQ(INPUT.vnl_in_h,1);
        EXPECT_EQ(INPUT.vh_in_h,1);
        EXPECT_EQ(INPUT.vion_in_h,1);
        EXPECT_EQ(INPUT.test_force,0);
        EXPECT_EQ(INPUT.test_stress,0);
        EXPECT_NEAR(INPUT.scf_thr,1.0e-8,1.0e-15);
        EXPECT_EQ(INPUT.scf_nmax,50);
        EXPECT_EQ(INPUT.relax_nmax,1);
        EXPECT_EQ(INPUT.out_stru,0);
        EXPECT_EQ(INPUT.occupations,"smearing");
        EXPECT_EQ(INPUT.smearing_method,"gauss");
        EXPECT_DOUBLE_EQ(INPUT.smearing_sigma,0.002);
        EXPECT_EQ(INPUT.mixing_mode,"broyden");
        EXPECT_DOUBLE_EQ(INPUT.mixing_beta,0.7);
        EXPECT_EQ(INPUT.mixing_ndim,8);
        EXPECT_DOUBLE_EQ(INPUT.mixing_gg0,0.00);
        EXPECT_EQ(INPUT.init_wfc,"atomic");
        EXPECT_EQ(INPUT.mem_saver,0);
        EXPECT_EQ(INPUT.printe,100);
        EXPECT_EQ(INPUT.init_chg,"atomic");
        EXPECT_EQ(INPUT.chg_extrap,"atomic");
        EXPECT_EQ(INPUT.out_freq_elec,0);
        EXPECT_EQ(INPUT.out_freq_ion,0);
        EXPECT_EQ(INPUT.out_chg,0);
        EXPECT_EQ(INPUT.out_dm,0);
        EXPECT_EQ(INPUT.out_dm1,0);
        EXPECT_EQ(INPUT.deepks_out_labels,0);
        EXPECT_EQ(INPUT.deepks_scf,0);
        EXPECT_EQ(INPUT.deepks_bandgap,0);
        EXPECT_EQ(INPUT.deepks_out_unittest,0);
        EXPECT_EQ(INPUT.out_pot,2);
        EXPECT_EQ(INPUT.out_wfc_pw,0);
        EXPECT_EQ(INPUT.out_wfc_r,0);
        EXPECT_EQ(INPUT.out_dos,0);
        EXPECT_EQ(INPUT.out_band,0);
        EXPECT_EQ(INPUT.out_proj_band,0);
        EXPECT_EQ(INPUT.out_mat_hs,0);
        EXPECT_EQ(INPUT.out_mat_hs2,0);
        EXPECT_EQ(INPUT.out_interval,1);
        EXPECT_EQ(INPUT.out_app_flag,0);
        EXPECT_EQ(INPUT.out_mat_r,0);
        EXPECT_FALSE(INPUT.out_wfc_lcao);
        EXPECT_FALSE(INPUT.out_alllog);
        EXPECT_DOUBLE_EQ(INPUT.dos_emin_ev,-15);
        EXPECT_DOUBLE_EQ(INPUT.dos_emax_ev,15);
        EXPECT_DOUBLE_EQ(INPUT.dos_edelta_ev,0.01);
        EXPECT_DOUBLE_EQ(INPUT.dos_scale,0.01);
        EXPECT_DOUBLE_EQ(INPUT.dos_sigma,0.07);
        EXPECT_FALSE(INPUT.out_element_info);
        EXPECT_DOUBLE_EQ(INPUT.lcao_ecut,20);
        EXPECT_DOUBLE_EQ(INPUT.lcao_dk,0.01);
        EXPECT_DOUBLE_EQ(INPUT.lcao_dr,0.01);
        EXPECT_DOUBLE_EQ(INPUT.lcao_rmax,30);
		EXPECT_TRUE(INPUT.bessel_nao_smooth);
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_sigma, 0.1);
		EXPECT_EQ(INPUT.bessel_nao_ecut, "default");
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_rcut, 6.0);
		EXPECT_DOUBLE_EQ(INPUT.bessel_nao_tolerence, 1E-12);
		EXPECT_EQ(INPUT.bessel_descriptor_lmax, 2);
		EXPECT_TRUE(INPUT.bessel_descriptor_smooth);
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_sigma, 0.1);
		EXPECT_EQ(INPUT.bessel_descriptor_ecut, "default");
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_rcut, 6.0);
		EXPECT_DOUBLE_EQ(INPUT.bessel_descriptor_tolerence, 1E-12);
        EXPECT_FALSE(INPUT.efield_flag);
        EXPECT_FALSE(INPUT.dip_cor_flag);
        EXPECT_EQ(INPUT.efield_dir,2);
        EXPECT_DOUBLE_EQ(INPUT.efield_pos_max,0.5);
        EXPECT_DOUBLE_EQ(INPUT.efield_pos_dec,0.1);
        EXPECT_DOUBLE_EQ(INPUT.efield_amp ,0.0);
        EXPECT_FALSE(INPUT.gate_flag);
        EXPECT_DOUBLE_EQ(INPUT.zgate,0.5);
        EXPECT_FALSE(INPUT.relax);
        EXPECT_FALSE(INPUT.block);
        EXPECT_DOUBLE_EQ(INPUT.block_down,0.45);
        EXPECT_DOUBLE_EQ(INPUT.block_up,0.55);
        EXPECT_DOUBLE_EQ(INPUT.block_height,0.1);
        EXPECT_EQ(INPUT.vdw_method,"d2");
        EXPECT_EQ(INPUT.vdw_s6,"default");
        EXPECT_EQ(INPUT.vdw_s8,"default");
        EXPECT_EQ(INPUT.vdw_a1,"default");
        EXPECT_EQ(INPUT.vdw_a2,"default");
        EXPECT_DOUBLE_EQ(INPUT.vdw_d,20);
        EXPECT_FALSE(INPUT.vdw_abc);
        EXPECT_EQ(INPUT.vdw_cutoff_radius,"default");
        EXPECT_EQ(INPUT.vdw_radius_unit,"Bohr");
        EXPECT_DOUBLE_EQ(INPUT.vdw_cn_thr,40.0);
        EXPECT_EQ(INPUT.vdw_cn_thr_unit,"Bohr");
        EXPECT_EQ(INPUT.vdw_C6_file,"default");
        EXPECT_EQ(INPUT.vdw_C6_unit,"Jnm6/mol");
        EXPECT_EQ(INPUT.vdw_R0_file,"default");
        EXPECT_EQ(INPUT.vdw_R0_unit,"A");
        EXPECT_EQ(INPUT.vdw_cutoff_type,"radius");
        EXPECT_EQ(INPUT.vdw_cutoff_period[0],3);
        EXPECT_EQ(INPUT.vdw_cutoff_period[1],3);
        EXPECT_EQ(INPUT.vdw_cutoff_period[2],3);
        EXPECT_EQ(INPUT.exx_hybrid_alpha,"default");
        EXPECT_EQ(INPUT.exx_real_number,"default");
        EXPECT_DOUBLE_EQ(INPUT.exx_hse_omega,0.11);
        EXPECT_TRUE(INPUT.exx_separate_loop);
        EXPECT_EQ(INPUT.exx_hybrid_step,100);
        EXPECT_DOUBLE_EQ(INPUT.exx_lambda,0.3);
		EXPECT_DOUBLE_EQ(INPUT.exx_mixing_beta,1.0);
        EXPECT_DOUBLE_EQ(INPUT.exx_pca_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_c_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_v_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_dm_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_schwarz_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_c_grad_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_v_grad_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_force_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_cauchy_stress_threshold,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_ccp_threshold,1E-8);
        EXPECT_EQ(INPUT.exx_ccp_rmesh_times,"default");
        EXPECT_EQ(INPUT.exx_distribute_type,"htime");
        EXPECT_EQ(INPUT.exx_opt_orb_lmax,0);
        EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_ecut,0.0);
        EXPECT_DOUBLE_EQ(INPUT.exx_opt_orb_tolerence,0.0);
        EXPECT_FALSE(INPUT.noncolin);
        EXPECT_FALSE(INPUT.lspinorb);
        EXPECT_DOUBLE_EQ(INPUT.soc_lambda,1.0);
        EXPECT_EQ(INPUT.input_error,0);
        EXPECT_DOUBLE_EQ(INPUT.td_force_dt,0.02);
        EXPECT_EQ(INPUT.td_vext,0);
        // EXPECT_EQ(INPUT.td_vext_dire,"1");
		EXPECT_EQ(INPUT.propagator,0);
		EXPECT_EQ(INPUT.td_stype,0);
		// EXPECT_EQ(INPUT.td_ttype,"0");
		EXPECT_EQ(INPUT.td_tstart,1);
		EXPECT_EQ(INPUT.td_tend,1000);
		EXPECT_EQ(INPUT.td_lcut1,0.05);
		EXPECT_EQ(INPUT.td_lcut2,0.95);
		// EXPECT_EQ(INPUT.td_gauss_amp,"0.25");
		// EXPECT_EQ(INPUT.td_gauss_freq,"22.13");
		// EXPECT_EQ(INPUT.td_gauss_phase,"0.0");
		// EXPECT_EQ(INPUT.td_gauss_t0,"100.0");
		// EXPECT_EQ(INPUT.td_gauss_sigma,"30.0");
		// EXPECT_EQ(INPUT.td_trape_amp,"2.74");
		// EXPECT_EQ(INPUT.td_trape_freq,"1.60");
		// EXPECT_EQ(INPUT.td_trape_phase,"0.0");
		// EXPECT_EQ(INPUT.td_trape_t1,"1875");
		// EXPECT_EQ(INPUT.td_trape_t2,"5625");
		// EXPECT_EQ(INPUT.td_trape_t3,"7500");
		// EXPECT_EQ(INPUT.td_trigo_freq1,"1.164656");
		// EXPECT_EQ(INPUT.td_trigo_freq2,"0.029116");
		// EXPECT_EQ(INPUT.td_trigo_phase1,"0.0");
		// EXPECT_EQ(INPUT.td_trigo_phase2,"0.0");
		// EXPECT_EQ(INPUT.td_trigo_amp,"2.74");
		// EXPECT_EQ(INPUT.td_heavi_t0,"100");	
		// EXPECT_EQ(INPUT.td_heavi_amp,"1.0");		
        EXPECT_EQ(INPUT.out_dipole,0);
		EXPECT_EQ(INPUT.out_efield,0);
		EXPECT_EQ(INPUT.td_print_eij,-1.0);
		EXPECT_EQ(INPUT.td_edm,0);
        EXPECT_DOUBLE_EQ(INPUT.cell_factor,1.2);
        EXPECT_EQ(INPUT.out_mul,0);
        EXPECT_FALSE(INPUT.restart_save);
        EXPECT_FALSE(INPUT.restart_load);
        EXPECT_FALSE(INPUT.test_skip_ewald);
        EXPECT_FALSE(INPUT.dft_plus_u);
        EXPECT_FALSE(INPUT.yukawa_potential);
        EXPECT_DOUBLE_EQ(INPUT.yukawa_lambda,-1.0);
        EXPECT_EQ(INPUT.omc,0);
        EXPECT_FALSE(INPUT.dft_plus_dmft);
        EXPECT_FALSE(INPUT.rpa);
        EXPECT_EQ(INPUT.coulomb_type,"full");
        EXPECT_EQ(INPUT.imp_sol,0);
        EXPECT_DOUBLE_EQ(INPUT.eb_k,80.0);
        EXPECT_DOUBLE_EQ(INPUT.tau,1.0798 * 1e-5);
        EXPECT_DOUBLE_EQ(INPUT.sigma_k,0.6);
        EXPECT_DOUBLE_EQ(INPUT.nc_k,0.00037);
        EXPECT_EQ(INPUT.of_kinetic,"vw");
        EXPECT_EQ(INPUT.of_method,"tn");
        EXPECT_EQ(INPUT.of_conv,"energy");
        EXPECT_DOUBLE_EQ(INPUT.of_tole,1e-6);
        EXPECT_DOUBLE_EQ(INPUT.of_tolp,1e-5);
        EXPECT_DOUBLE_EQ(INPUT.of_tf_weight,1.);
        EXPECT_DOUBLE_EQ(INPUT.of_vw_weight,1.);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_alpha,0.833333);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_beta,0.833333);
        EXPECT_DOUBLE_EQ(INPUT.of_wt_rho0,1.);
        EXPECT_FALSE(INPUT.of_hold_rho0);
        EXPECT_DOUBLE_EQ(INPUT.of_lkt_a, 1.3);
        EXPECT_FALSE(INPUT.of_full_pw);
        EXPECT_EQ(INPUT.of_full_pw_dim,0);
        EXPECT_FALSE(INPUT.of_read_kernel);
        EXPECT_EQ(INPUT.of_kernel_file,"WTkernel.txt");
        EXPECT_EQ(INPUT.device, "cpu");
        EXPECT_EQ(INPUT.ncx,0);
        EXPECT_EQ(INPUT.ncy,0);
        EXPECT_EQ(INPUT.ncz,0);
        //EXPECT_NEAR(INPUT.force_thr_ev,0.0257112,1e-8);
        EXPECT_DOUBLE_EQ(INPUT.hubbard_u[0],0);
	EXPECT_EQ(INPUT.orbital_corr[0],-1);
	EXPECT_NEAR(INPUT.mdp.lj_epsilon,0.01032,1e-7);
	EXPECT_NEAR(INPUT.mdp.lj_rcut,8.5,1e-7);
	EXPECT_NEAR(INPUT.mdp.lj_sigma,3.405,1e-7);
	EXPECT_EQ(INPUT.mdp.md_damp,1);
	EXPECT_EQ(INPUT.mdp.md_dt,1);
	EXPECT_EQ(INPUT.mdp.md_dumpfreq,1);
	EXPECT_EQ(INPUT.mdp.md_nraise,1);
	EXPECT_EQ(INPUT.cal_syns,0);
	EXPECT_EQ(INPUT.dmax,0.01);
	EXPECT_EQ(INPUT.mdp.md_nstep,10);
	EXPECT_EQ(INPUT.mdp.md_pchain,1);
	EXPECT_EQ(INPUT.mdp.md_pcouple,"none");
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfirst,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfreq,0);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_plast,-1);
	EXPECT_EQ(INPUT.mdp.md_pmode,"iso");
	EXPECT_EQ(INPUT.mdp.md_restart,0);
	EXPECT_EQ(INPUT.mdp.md_restartfreq,5);
	EXPECT_EQ(INPUT.mdp.md_seed,-1);
    EXPECT_EQ(INPUT.mdp.md_prec_level, 2);
    EXPECT_DOUBLE_EQ(INPUT.ref_cell_factor,1.2);
	EXPECT_EQ(INPUT.mdp.md_tchain,1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfirst,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfreq,0);
	EXPECT_EQ(INPUT.mdp.md_thermostat,"nhc");
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tlast,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_tolerance,100);
	EXPECT_EQ(INPUT.mdp.md_type,"nvt");
	EXPECT_EQ(INPUT.mdp.msst_direction,2);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_qmass,-1);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_tscale,0.01);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vel,0);
	EXPECT_DOUBLE_EQ(INPUT.mdp.msst_vis,0);
	EXPECT_EQ(INPUT.mdp.pot_file,"graph.pb");
    EXPECT_FALSE(INPUT.mdp.dump_force);
    EXPECT_FALSE(INPUT.mdp.dump_vel);
    EXPECT_FALSE(INPUT.mdp.dump_virial);
    EXPECT_EQ(INPUT.sc_mag_switch, 0);
    EXPECT_TRUE(INPUT.decay_grad_switch);
    EXPECT_DOUBLE_EQ(INPUT.sc_thr, 1e-4);
    EXPECT_EQ(INPUT.nsc, 50);
	EXPECT_EQ(INPUT.nsc_min, 4);
	EXPECT_EQ(INPUT.sc_scf_nmin, 4);
    EXPECT_DOUBLE_EQ(INPUT.alpha_trial, 0.02);
	EXPECT_DOUBLE_EQ(INPUT.sccut, 4.0);
    EXPECT_EQ(INPUT.sc_file, "sc.json");
}

TEST_F(InputTest, Default_2)
{
	//==================================================
	// prepare default parameters for the 1st calling
	EXPECT_EQ(INPUT.vdw_method,"d2");
	EXPECT_EQ(INPUT.vdw_s6,"default");
        EXPECT_EQ(INPUT.vdw_s8,"default");
        EXPECT_EQ(INPUT.vdw_a1,"default");
        EXPECT_EQ(INPUT.vdw_a2,"default");
        EXPECT_EQ(INPUT.vdw_cutoff_radius,"default");
	EXPECT_NE(INPUT.esolver_type,"sdft");
	EXPECT_NE(INPUT.method_sto,1);
	EXPECT_NE(INPUT.method_sto,2);
	EXPECT_NE(INPUT.of_wt_rho0,0.0);
	EXPECT_EQ(INPUT.exx_hybrid_alpha,"default");
	EXPECT_EQ(INPUT.dft_functional,"hse");
	EXPECT_EQ(INPUT.exx_real_number,"default");
	EXPECT_TRUE(INPUT.gamma_only);
        EXPECT_EQ(INPUT.exx_ccp_rmesh_times,"default");
	EXPECT_EQ(INPUT.diago_proc,4);
	EXPECT_EQ(GlobalV::NPROC,1);
	EXPECT_EQ(INPUT.calculation,"scf");
	EXPECT_EQ(INPUT.basis_type,"lcao");
	INPUT.ks_solver = "default";
	INPUT.lcao_ecut = 0;
	INPUT.scf_thr = -1.0;
	INPUT.scf_thr_type = -1;
    EXPECT_DOUBLE_EQ(INPUT.ecutwfc, 20.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_ecut, 20.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_height, 20.0);
    EXPECT_DOUBLE_EQ(INPUT.erf_sigma, 4.0);
    INPUT.nbndsto_str = "all";
    INPUT.nx = INPUT.ny = INPUT.nz = 4;
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 0;
    // the 1st calling
    INPUT.Default_2();
    // ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.ndx, 4);
    EXPECT_EQ(INPUT.ndy, 4);
    EXPECT_EQ(INPUT.ndz, 4);
    EXPECT_FALSE(GlobalV::double_grid);
    EXPECT_DOUBLE_EQ(INPUT.ecutrho, 80.0);
    EXPECT_EQ(INPUT.vdw_s6, "0.75");
    EXPECT_EQ(INPUT.vdw_cutoff_radius, "56.6918");
    EXPECT_EQ(INPUT.bndpar,1);
	EXPECT_EQ(INPUT.method_sto,2);
	EXPECT_TRUE(INPUT.of_hold_rho0);
        EXPECT_EQ(INPUT.of_full_pw_dim,0);
	EXPECT_EQ(INPUT.exx_hybrid_alpha,"0.25");
	EXPECT_EQ(INPUT.exx_real_number,"1");
        EXPECT_EQ(INPUT.exx_ccp_rmesh_times,"1.5");
	EXPECT_EQ(INPUT.diago_proc,1);
	EXPECT_EQ(INPUT.mem_saver,0);
	EXPECT_EQ(INPUT.relax_nmax,1);
	EXPECT_DOUBLE_EQ(INPUT.scf_thr,1.0e-7);
	EXPECT_EQ(INPUT.scf_thr_type,2);
#ifdef __ELPA
	EXPECT_EQ(INPUT.ks_solver,"genelpa");
#else
	EXPECT_EQ(INPUT.ks_solver,"scalapack_gvx");
#endif
        EXPECT_DOUBLE_EQ(INPUT.lcao_ecut,20.0);
	EXPECT_EQ(INPUT.nbands_sto, 0);
	//==================================================
	// prepare default parameters for the 2nd calling
	INPUT.vdw_method = "d3_0";
	INPUT.vdw_s6 = "default";
	INPUT.vdw_s8 = "default";
	INPUT.vdw_a1 = "default";
	INPUT.vdw_a2 = "default";
	INPUT.vdw_cutoff_radius = "default";
	INPUT.exx_hybrid_alpha = "default";
	INPUT.dft_functional = "hf";
	INPUT.exx_real_number = "default";
	INPUT.gamma_only = 0;
        INPUT.exx_ccp_rmesh_times = "default";
	INPUT.diago_proc = 0;
	INPUT.calculation = "relax";
    INPUT.chg_extrap = "default";
    INPUT.relax_nmax = 0;
    INPUT.basis_type = "pw";
    INPUT.ks_solver = "default";
    INPUT.gamma_only_local = 1;
	INPUT.scf_thr = -1.0;
	INPUT.scf_thr_type = -1;
    INPUT.nbndsto_str = "0";
    INPUT.esolver_type = "sdft";
    INPUT.nx = INPUT.ny = INPUT.nz = 0;
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 4;
    // the 2nd calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.nx, 4);
    EXPECT_EQ(INPUT.ny, 4);
    EXPECT_EQ(INPUT.nz, 4);
    EXPECT_FALSE(GlobalV::double_grid);
    EXPECT_EQ(INPUT.chg_extrap, "first-order");
    EXPECT_EQ(INPUT.vdw_s6, "1.0");
    EXPECT_EQ(INPUT.vdw_s8, "0.722");
    EXPECT_EQ(INPUT.vdw_a1, "1.217");
    EXPECT_EQ(INPUT.vdw_a2,"1.0");
	EXPECT_EQ(INPUT.vdw_cutoff_radius,"95");
	EXPECT_EQ(INPUT.exx_hybrid_alpha,"1");
	EXPECT_EQ(INPUT.exx_real_number,"0");
        EXPECT_EQ(INPUT.exx_ccp_rmesh_times,"5");
	EXPECT_EQ(INPUT.diago_proc,1);
	EXPECT_EQ(INPUT.mem_saver,0);
	EXPECT_EQ(INPUT.cal_force,1);
	EXPECT_EQ(INPUT.relax_nmax,50);
	EXPECT_EQ(INPUT.ks_solver,"cg");
	EXPECT_EQ(INPUT.gamma_only_local,0);
	EXPECT_EQ(INPUT.bx,1);
	EXPECT_EQ(INPUT.by,1);
	EXPECT_EQ(INPUT.bz,1);
	EXPECT_DOUBLE_EQ(INPUT.scf_thr,1.0e-9);
	EXPECT_EQ(INPUT.scf_thr_type,1);
	EXPECT_EQ(INPUT.esolver_type, "ksdft");
	//==================================================
	// prepare default parameters for the 3rd calling
	INPUT.vdw_method = "d3_bj";
	INPUT.vdw_s6 = "default";
	INPUT.vdw_s8 = "default";
	INPUT.vdw_a1 = "default";
	INPUT.vdw_a2 = "default";
	INPUT.vdw_cutoff_radius = "default";
	INPUT.calculation = "get_S";
    INPUT.chg_extrap = "default";
    INPUT.basis_type = "pw";
    INPUT.pw_diag_thr = 1.0e-2;
    INPUT.cal_force = 1;
	INPUT.init_chg = "atomic";
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "cg";
	GlobalV::NPROC = 8;
	INPUT.diago_proc = 1;
    INPUT.nx = INPUT.ny = INPUT.nz = 4;
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 6;
    // the 3rd calling
    INPUT.Default_2();
    // ^^^^^^^^^^^^^^
    EXPECT_TRUE(GlobalV::double_grid);
    EXPECT_EQ(INPUT.chg_extrap, "atomic");
    EXPECT_EQ(INPUT.vdw_s6, "1.0");
    EXPECT_EQ(INPUT.vdw_s8, "0.7875");
    EXPECT_EQ(INPUT.vdw_a1,"0.4289");
	EXPECT_EQ(INPUT.vdw_a2,"4.4407");
	EXPECT_EQ(INPUT.vdw_cutoff_radius,"95");
	EXPECT_EQ(GlobalV::CALCULATION,"nscf");
	EXPECT_EQ(INPUT.relax_nmax,1);
	EXPECT_EQ(INPUT.out_stru,0);
	EXPECT_DOUBLE_EQ(INPUT.pw_diag_thr,1.0e-5);
	EXPECT_FALSE(INPUT.cal_force);
	EXPECT_EQ(INPUT.init_chg,"file");
	EXPECT_EQ(INPUT.diago_proc,8);
	//==================================================
	// prepare default parameters for the 4th calling
	INPUT.calculation = "get_pchg";
    INPUT.chg_extrap = "default";
    INPUT.symmetry = "default";
    INPUT.ecutwfc = 10;
    INPUT.ecutrho = 100;
    INPUT.nx = INPUT.ny = INPUT.nz = 0;
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 0;
    GlobalV::double_grid = false;
    // the 4th calling
    INPUT.Default_2();
    // ^^^^^^^^^^^^^^
    EXPECT_TRUE(GlobalV::double_grid);
    EXPECT_EQ(GlobalV::CALCULATION, "get_pchg");
    EXPECT_EQ(INPUT.relax_nmax, 1);
    EXPECT_EQ(INPUT.out_stru, 0);
    EXPECT_EQ(INPUT.symmetry, "0");
	EXPECT_EQ(INPUT.out_band,0);
	EXPECT_EQ(INPUT.out_proj_band,0);
	EXPECT_EQ(INPUT.cal_force,0);
	EXPECT_EQ(INPUT.init_wfc,"file");
	EXPECT_EQ(INPUT.init_chg,"atomic");
	EXPECT_EQ(INPUT.chg_extrap,"atomic");
	EXPECT_EQ(INPUT.out_chg,1);
	EXPECT_EQ(INPUT.out_dm,0);
	EXPECT_EQ(INPUT.out_dm1,0);
	EXPECT_EQ(INPUT.out_pot,0);
	//==================================================
	// prepare default parameters for the 5th calling
	INPUT.calculation = "get_wf";
    INPUT.symmetry = "default";
    INPUT.chg_extrap = "default";
    // the 5th calling
    INPUT.Default_2();
    // ^^^^^^^^^^^^^^
	EXPECT_EQ(GlobalV::CALCULATION,"get_wf");
    EXPECT_EQ(INPUT.relax_nmax, 1);
    EXPECT_EQ(INPUT.symmetry, "0");
    EXPECT_EQ(INPUT.out_stru, 0);
	EXPECT_EQ(INPUT.out_band,0);
	EXPECT_EQ(INPUT.out_proj_band,0);
	EXPECT_EQ(INPUT.cal_force,0);
	EXPECT_EQ(INPUT.init_wfc,"file");
	EXPECT_EQ(INPUT.init_chg,"atomic");
	EXPECT_EQ(INPUT.chg_extrap,"atomic");
	EXPECT_EQ(INPUT.out_chg,1);
	EXPECT_EQ(INPUT.out_dm,0);
	EXPECT_EQ(INPUT.out_dm1,0);
	EXPECT_EQ(INPUT.out_pot,0);
	//==================================================
	// prepare default parameters for the 6th calling
	INPUT.calculation = "md";
    INPUT.chg_extrap = "default";
    INPUT.mdp.md_nstep = 0;
	INPUT.out_md_control = 0;
	INPUT.mdp.md_tlast = -1.0;
	INPUT.mdp.md_plast = -1.0;
	INPUT.mdp.md_tfreq = 0;
	INPUT.mdp.md_pfreq = 0;
	INPUT.mdp.md_restart = 1;
	INPUT.mdp.md_type = "npt";
	INPUT.mdp.md_pmode = "iso";
	// the 6th calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
	EXPECT_EQ(GlobalV::CALCULATION,"md");
    EXPECT_EQ(INPUT.chg_extrap, "second-order");
    EXPECT_EQ(INPUT.symmetry, "0");
    EXPECT_EQ(INPUT.cal_force, 1);
    EXPECT_EQ(INPUT.mdp.md_nstep,50);
    EXPECT_EQ(INPUT.out_level, "m");
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_plast, INPUT.mdp.md_pfirst);
    EXPECT_DOUBLE_EQ(INPUT.mdp.md_tfreq,1.0/40/INPUT.mdp.md_dt);
	EXPECT_DOUBLE_EQ(INPUT.mdp.md_pfreq,1.0/400/INPUT.mdp.md_dt);
	EXPECT_EQ(INPUT.init_vel,1);
	EXPECT_EQ(INPUT.cal_stress,1);
	//==================================================
	// prepare default parameters for the 7th calling
	INPUT.calculation = "cell-relax";
    INPUT.chg_extrap = "default";
    INPUT.relax_nmax = 0;
    // the 7th calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.chg_extrap, "first-order");
    EXPECT_EQ(INPUT.cal_force, 1);
    EXPECT_EQ(INPUT.cal_stress,1);
	EXPECT_EQ(INPUT.relax_nmax,50);
	//==================================================
	// prepare default parameters for the 8th calling
	INPUT.calculation = "test_memory";
    INPUT.chg_extrap = "default";
    // the 8th calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.chg_extrap, "atomic");
    EXPECT_EQ(INPUT.relax_nmax,1);
	//==================================================
	// prepare default parameters for the 9th calling
	INPUT.calculation = "test_neighbour";
    INPUT.chg_extrap = "default";
    // the 9th calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.chg_extrap, "atomic");
    EXPECT_EQ(INPUT.relax_nmax,1);
	//==================================================
	// prepare default parameters for the 10th calling
	INPUT.calculation = "gen_bessel";
    INPUT.chg_extrap = "default";
    // the 10th calling
	INPUT.Default_2();
	// ^^^^^^^^^^^^^^
    EXPECT_EQ(INPUT.chg_extrap, "atomic");
    EXPECT_EQ(INPUT.relax_nmax,1);
	//==================================================
	remove("INPUT");
	remove("STRU");
}

TEST_F(InputTest, Check)
{
    INPUT.ecutwfc = 20.0;
    INPUT.ecutrho = 10;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(INPUT.Check(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("ecutrho/ecutwfc must >= 4"));
    INPUT.ecutrho = 100.0;
    INPUT.nx = INPUT.ny = INPUT.nz = 10;
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 8;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(INPUT.Check(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("smooth grids is denser than dense grids"));
    INPUT.ndx = INPUT.ndy = INPUT.ndz = 11;
    //
    INPUT.nbands = -1;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("NBANDS must >= 0"));
	INPUT.nbands = 2;
	//
	INPUT.nb2d = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nb2d must > 0"));
	INPUT.nb2d = 1;
	//
	INPUT.ntype = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("ntype must > 0"));
	INPUT.ntype = 1;
	//
	INPUT.basis_type = "lcao";
	INPUT.diago_proc = 2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("please don't set diago_proc with lcao base"));
	INPUT.diago_proc = 1;
	//
	INPUT.kspacing[0] = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("kspacing must > 0"));
	INPUT.kspacing[0] = INPUT.kspacing[1] = INPUT.kspacing[2] = 0.8;
	//
	INPUT.nelec = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nelec < 0 is not allowed !"));
	INPUT.nelec = 100;
	//
	INPUT.efield_flag = 0;
	INPUT.dip_cor_flag = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("dipole correction is not active if efield_flag=false !"));
	INPUT.dip_cor_flag = 0;
	//
	INPUT.efield_flag = 1;
	INPUT.gate_flag = 1;
	INPUT.dip_cor_flag = 0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("gate field cannot be used with efield if dip_cor_flag=false !"));
	INPUT.gate_flag = 0;
	//
	INPUT.calculation = "nscf";
	INPUT.out_dos = 3;
	INPUT.symmetry = "1";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("symmetry can't be used for out_dos==3(Fermi Surface Plotting) by now."));
	INPUT.symmetry = "0";
	INPUT.out_dos = 0;
	//
	INPUT.calculation = "get_pchg";
	INPUT.basis_type = "pw";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculate = get_pchg is only availble for LCAO"));
	INPUT.basis_type = "lcao";
	//
	INPUT.calculation = "get_wf";
	INPUT.basis_type = "pw";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculate = get_wf is only availble for LCAO"));
	INPUT.basis_type = "lcao";
	//
	INPUT.calculation = "md";
	INPUT.mdp.md_dt = -1.0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("time interval of MD calculation should be set!"));
	INPUT.mdp.md_dt = 1.0;
    //
    INPUT.mdp.md_type = "npt";
	INPUT.mdp.md_pmode = "iso";
	INPUT.mdp.md_pfirst = -1.0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("pressure of MD calculation should be set!"));
	INPUT.mdp.md_pfirst = 1.0;
	//
	INPUT.mdp.md_type = "msst";
	INPUT.mdp.msst_qmass = -1.0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("msst_qmass must be greater than 0!"));
	INPUT.mdp.msst_qmass = 1.0;
	//
	INPUT.esolver_type = "dp";
	INPUT.mdp.pot_file = "graph.pb";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Can not find DP model !"));
	INPUT.esolver_type = "ksdft";
	//
	INPUT.calculation = "gen_bessel";
	INPUT.basis_type = "lcao";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("to generate descriptors, please use pw basis"));
	INPUT.basis_type = "pw";
	//
	INPUT.calculation = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("check 'calculation' !"));
	INPUT.calculation = "scf";
	//
	INPUT.init_chg = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("wrong 'init_chg',not 'atomic', 'file',please check"));
	INPUT.init_chg = "atomic";
	//
	INPUT.gamma_only_local = 0;
	INPUT.out_dm = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("out_dm with k-point algorithm is not implemented yet."));
	INPUT.out_dm = 0;
	//
	INPUT.gamma_only_local = 1;
	INPUT.out_dm1 = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("out_dm1 is only for multi-k"));
	INPUT.gamma_only_local = 0;
	INPUT.out_dm1 = 0;
	//
	INPUT.basis_type = "pw";
	INPUT.chg_extrap = "dm";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("wrong 'chg_extrap=dm' is only available for local orbitals."));
	INPUT.chg_extrap = "atomic";
	//
	INPUT.nbands = 100001;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nbnd >100000, out of range"));
	INPUT.nbands = 100;
	//
	INPUT.nelec = 2*INPUT.nbands + 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nelec > 2*nbnd , bands not enough!"));
	INPUT.nelec = INPUT.nbands;
	//
	INPUT.nspin = 3;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nspin does not equal to 1, 2, or 4!"));
	INPUT.nspin = 1;
	//
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "genelpa";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("genelpa can not be used with plane wave basis."));
	//
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "scalapack_gvx";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("scalapack_gvx can not be used with plane wave basis."));
	//
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "lapack";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("lapack can not be used with plane wave basis."));
	//
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("please check the ks_solver parameter!"));
	INPUT.ks_solver = "cg";
	//
	INPUT.basis_type = "pw";
	INPUT.gamma_only = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("gamma_only not implemented for plane wave now."));
	INPUT.gamma_only = 0;
	//
	INPUT.basis_type = "pw";
	INPUT.out_proj_band = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("out_proj_band not implemented for plane wave now."));
	INPUT.out_proj_band = 0;
	//
	INPUT.basis_type = "pw";
	INPUT.out_dos = 3;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Fermi Surface Plotting not implemented for plane wave now."));
	INPUT.out_dos = 0;
	//
	INPUT.basis_type = "pw";
	INPUT.sc_mag_switch = 3;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Non-colliner Spin-constrained DFT not implemented for plane wave now."));
	INPUT.sc_mag_switch = 0;
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "cg";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("not ready for cg method in lcao ."));
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "genelpa";
#ifndef __MPI
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("genelpa can not be used for series version."));
#endif
#ifndef __ELPA
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Can not use genelpa if abacus is not compiled with ELPA. Please change ks_solver to scalapack_gvx."));
#endif
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "scalapack_gvx";
#ifndef __MPI
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("scalapack_gvx can not be used for series version."));
#endif
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "lapack";
#ifdef __MPI
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("ks_solver=lapack is not an option for parallel version of ABACUS (try genelpa)"));
#endif
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "cusolver";
#ifndef __MPI
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("Cusolver can not be used for series version."));
#endif
	//
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("please check the ks_solver parameter!"));
	INPUT.ks_solver = "genelpa";
	//
	INPUT.basis_type = "lcao";
	INPUT.kpar = 2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("kpar > 1 has not been supported for lcao calculation."));
	INPUT.kpar = 1;
	//
	INPUT.basis_type = "lcao";
	INPUT.out_wfc_lcao = 3;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("out_wfc_lcao must be 0, 1, or 2"));
	INPUT.out_wfc_lcao = 0;
	//
	INPUT.basis_type = "lcao_in_pw";
	INPUT.ks_solver = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("LCAO in plane wave can only done with lapack."));
	INPUT.ks_solver = "default";
	//
	INPUT.basis_type = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("please check the basis_type parameter!"));
	INPUT.basis_type = "pw";
	//
	INPUT.relax_method = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("relax_method can only be sd, cg, bfgs or cg_bfgs."));
	INPUT.relax_method = "cg";
	//
	INPUT.bx = 11; INPUT.by = 1; INPUT.bz = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("bx, or by, or bz is larger than 10!"));
	INPUT.bx = 1;
	//
	INPUT.vdw_method = "d2";
	INPUT.vdw_C6_unit = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_C6_unit must be Jnm6/mol or eVA6"));
	INPUT.vdw_C6_unit = "eVA6";
	//
	INPUT.vdw_R0_unit = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_R0_unit must be A or Bohr"));
	INPUT.vdw_R0_unit = "A";
	//
	INPUT.vdw_cutoff_type = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_type must be radius or period"));
	INPUT.vdw_cutoff_type = "radius";
	//
	INPUT.vdw_cutoff_period.x = 0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_period <= 0 is not allowd"));
	INPUT.vdw_cutoff_period.x = 3;
	//
	INPUT.vdw_cutoff_radius = "0";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_cutoff_radius <= 0 is not allowd"));
	INPUT.vdw_cutoff_radius = "1.0";
	//
	INPUT.vdw_radius_unit = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_radius_unit must be A or Bohr"));
	INPUT.vdw_radius_unit = "A";
	//
	INPUT.vdw_cn_thr = 0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_cn_thr <= 0 is not allowd"));
	INPUT.vdw_cn_thr = 1.0;
	//
	INPUT.vdw_cn_thr_unit = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("vdw_cn_thr_unit must be A or Bohr"));
	INPUT.vdw_cn_thr_unit = "A";
	//
	INPUT.dft_functional = "scan0";
	INPUT.exx_hybrid_alpha = "1.25";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("must 0 <= exx_hybrid_alpha <= 1"));
	INPUT.exx_hybrid_alpha = "0.25";
	//
	INPUT.exx_hybrid_step = 0;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("must exx_hybrid_step > 0"));
	INPUT.exx_hybrid_step = 1;
	//
	INPUT.exx_ccp_rmesh_times = "-1";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("must exx_ccp_rmesh_times >= 1"));
	INPUT.exx_ccp_rmesh_times = "1.5";
	//
	INPUT.exx_distribute_type = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("exx_distribute_type must be htime or kmeans2 or kmeans1"));
	INPUT.exx_distribute_type = "htime";
	//
	INPUT.dft_functional = "opt_orb";
	INPUT.exx_opt_orb_lmax = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_lmax must >=0"));
	INPUT.exx_opt_orb_lmax = 0;
	//
	INPUT.exx_opt_orb_ecut = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_ecut must >=0"));
	INPUT.exx_opt_orb_ecut = 0;
	//
	INPUT.exx_opt_orb_tolerence = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("exx_opt_orb_tolerence must >=0"));
	INPUT.exx_opt_orb_tolerence = 0;
	//
	INPUT.berry_phase = 1;
	INPUT.basis_type = "lcao_in_pw";
	INPUT.ks_solver = "lapack";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set basis_type = pw or lcao"));
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "cg";
	//
	INPUT.calculation = "scf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set calculation = nscf"));
	INPUT.calculation = "nscf";
	//
	INPUT.gdir = 4;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculate berry phase, please set gdir = 1 or 2 or 3"));
	INPUT.gdir = 3;
	INPUT.berry_phase = 0;
	//
	INPUT.towannier90 = 1;
	// due to the repair of lcao_in_pw, original warning has been deprecated, 2023/12/23, ykhuang
	// INPUT.basis_type = "lcao_in_pw";
	// INPUT.ks_solver = "lapack";
	// testing::internal::CaptureStdout();
	// EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	// output = testing::internal::GetCapturedStdout();
	// EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set basis_type = pw or lcao"));
	INPUT.basis_type = "pw";
	INPUT.ks_solver = "cg";
	//
	INPUT.calculation = "scf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set calculation = nscf"));
	INPUT.calculation = "nscf";
	//
	INPUT.nspin = 2;
	INPUT.wannier_spin = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("to use towannier90, please set wannier_spin = up or down"));
	INPUT.wannier_spin = "up";
	INPUT.towannier90 = 0;
	//
	INPUT.read_file_dir = "arbitrary";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("please set right files directory for reading in."));
	INPUT.read_file_dir = "auto";
	// Start to check deltaspin parameters
	INPUT.sc_mag_switch = 1;
	INPUT.sc_file = "none";
	INPUT.basis_type = "lcao";
	INPUT.ks_solver = "genelpa";
	// warning 1 of Deltaspin
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("sc_file (json format) must be set when sc_mag_switch > 0"));
	// warning 2 of Deltaspin
	INPUT.sc_file = "sc.json";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("sc_file does not exist"));
	INPUT.sc_file = "./support/sc.json";
	// warning 3 of Deltaspin
	INPUT.nspin = 1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nspin must be 2 or 4 when sc_mag_switch > 0"));
	INPUT.nspin = 4;
	// warning 4 of Deltaspin
	INPUT.calculation = "nscf";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("calculation must be scf when sc_mag_switch > 0"));
	INPUT.calculation = "scf";
	// warning 5 of Deltaspin
	INPUT.sc_thr = -1;
		testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("sc_thr must > 0"));
	INPUT.sc_thr = 1e-6;
	// warning 6 of Deltaspin
	INPUT.nsc = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nsc must > 0"));
	INPUT.nsc = 100;
	// warning 7 of Deltaspin
	INPUT.nsc_min = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("nsc_min must > 0"));
	INPUT.nsc_min = 2;
	// warning 8 of Deltapsin
    INPUT.alpha_trial = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("alpha_trial must > 0"));
	INPUT.alpha_trial = 0.01;
	// warning 9 of Deltapsin
    INPUT.sccut = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("sccut must > 0"));
	INPUT.sccut = 3.0;
	// warning 10 of Deltaspin
	INPUT.sc_scf_nmin = -1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("sc_scf_nmin must >= 2"));
	INPUT.sc_scf_nmin = 2;
    // restore to default values
    INPUT.nspin = 1;
	INPUT.sc_file = "none";
	INPUT.sc_mag_switch = 0;
	INPUT.ks_solver = "default";
	INPUT.basis_type = "pw";
	// End of checking Deltaspin parameters

	/*
	testing::internal::CaptureStdout();
	EXPECT_EXIT(INPUT.Check(),::testing::ExitedWithCode(0), "");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr(""));
	*/
}


#undef private


class ReadKSpacingTest : public ::testing::Test {
protected:
    void SetUp() override 
	{
        // create a temporary file for testing
        char tmpname[] = "tmpfile.tmp";
        int fd = mkstemp(tmpname);
        tmpfile = tmpname;
        std::ofstream ofs(tmpfile);
        ofs << "1.0"; // valid input
        ofs.close();
    }

    void TearDown() override {
        close(fd);
        unlink(tmpfile.c_str());
    }

    std::string tmpfile;
    int fd;
};

TEST_F(ReadKSpacingTest, ValidInputOneValue) {
    std::ifstream ifs(tmpfile);
    EXPECT_NO_THROW(INPUT.read_kspacing(ifs));
    EXPECT_EQ(INPUT.kspacing[0], 1.0);
    EXPECT_EQ(INPUT.kspacing[1], 1.0);
    EXPECT_EQ(INPUT.kspacing[2], 1.0);
}

TEST_F(ReadKSpacingTest, ValidInputThreeValue) {
	std::ofstream ofs(tmpfile);
    ofs << "1.0 2.0 3.0"; // invalid input
    ofs.close();

    std::ifstream ifs(tmpfile);
    EXPECT_NO_THROW(INPUT.read_kspacing(ifs));
    EXPECT_EQ(INPUT.kspacing[0], 1.0);
    EXPECT_EQ(INPUT.kspacing[1], 2.0);
    EXPECT_EQ(INPUT.kspacing[2], 3.0);
}

TEST_F(ReadKSpacingTest, InvalidInput) {
    std::ofstream ofs(tmpfile);
    ofs << "1.0 2.0"; // invalid input
    ofs.close();

    std::ifstream ifs(tmpfile);
	testing::internal::CaptureStdout();
	INPUT.read_kspacing(ifs);
	std::string output;
	output = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(ifs.fail());
}

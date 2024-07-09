#include <sstream>
#include <string>

#include "input.h"
#include "input_conv.h"
#include "module_parameter/parameter.h"

void Input_Conv::tmp_convert()
{
    INPUT.stru_file = PARAM.inp.stru_file;
    INPUT.pseudo_dir = PARAM.inp.pseudo_dir;
    INPUT.orbital_dir = PARAM.inp.orbital_dir;
    INPUT.read_file_dir = PARAM.inp.read_file_dir;
    INPUT.kpoint_file = PARAM.inp.kpoint_file;
    INPUT.wannier_card = PARAM.inp.wannier_card;
    INPUT.latname = PARAM.inp.latname;
    INPUT.calculation = PARAM.inp.calculation;
    INPUT.esolver_type = PARAM.inp.esolver_type;
    INPUT.pseudo_rcut = PARAM.inp.pseudo_rcut;
    INPUT.pseudo_mesh = PARAM.inp.pseudo_mesh;
    INPUT.ntype = PARAM.inp.ntype;
    INPUT.nbands = PARAM.inp.nbands;
    INPUT.nbands_istate = PARAM.inp.nbands_istate;
    INPUT.pw_seed = PARAM.inp.pw_seed;
    INPUT.init_vel = PARAM.inp.init_vel;
    INPUT.ref_cell_factor = PARAM.inp.ref_cell_factor;

    INPUT.kpar = PARAM.inp.kpar;
    INPUT.berry_phase = PARAM.inp.berry_phase;
    INPUT.gdir = PARAM.inp.gdir;
    INPUT.kspacing[0] = PARAM.inp.kspacing[0];
    INPUT.kspacing[1] = PARAM.inp.kspacing[1];
    INPUT.kspacing[2] = PARAM.inp.kspacing[2];
    INPUT.min_dist_coef = PARAM.inp.min_dist_coef;
    INPUT.towannier90 = PARAM.inp.towannier90;
    INPUT.nnkpfile = PARAM.inp.nnkpfile;
    INPUT.wannier_spin = PARAM.inp.wannier_spin;
    INPUT.wannier_method = PARAM.inp.wannier_method;
    INPUT.out_wannier_mmn = PARAM.inp.out_wannier_mmn;
    INPUT.out_wannier_amn = PARAM.inp.out_wannier_amn;
    INPUT.out_wannier_unk = PARAM.inp.out_wannier_unk;
    INPUT.out_wannier_eig = PARAM.inp.out_wannier_eig;
    INPUT.out_wannier_wvfn_formatted = PARAM.inp.out_wannier_wvfn_formatted;

    INPUT.nche_sto = PARAM.inp.nche_sto;
    INPUT.nbands_sto = PARAM.inp.nbands_sto;
    INPUT.seed_sto = PARAM.inp.seed_sto;
    INPUT.initsto_ecut = PARAM.inp.initsto_ecut;
    INPUT.emax_sto = PARAM.inp.emax_sto;
    INPUT.emin_sto = PARAM.inp.emin_sto;
    INPUT.bndpar = PARAM.inp.bndpar;
    INPUT.initsto_freq = PARAM.inp.initsto_freq;
    INPUT.method_sto = PARAM.inp.method_sto;
    INPUT.npart_sto = PARAM.inp.npart_sto;
    INPUT.cal_cond = PARAM.inp.cal_cond;
    INPUT.cond_che_thr = PARAM.inp.cond_che_thr;
    INPUT.cond_smear = PARAM.inp.cond_smear;
    INPUT.cond_dw = PARAM.inp.cond_dw;
    INPUT.cond_wcut = PARAM.inp.cond_wcut;
    INPUT.cond_dt = PARAM.inp.cond_dt;
    INPUT.cond_dtbatch = PARAM.inp.cond_dtbatch;
    INPUT.cond_fwhm = PARAM.inp.cond_fwhm;
    INPUT.cond_nonlocal = PARAM.inp.cond_nonlocal;

    INPUT.dft_functional = PARAM.inp.dft_functional;
    INPUT.xc_temperature = PARAM.inp.xc_temperature;
    INPUT.nspin = PARAM.inp.nspin;
    INPUT.lmaxmax = PARAM.inp.lmaxmax;

    INPUT.basis_type = PARAM.inp.basis_type;
    INPUT.ks_solver = PARAM.inp.ks_solver;
    INPUT.cal_force = PARAM.inp.cal_force;
    INPUT.force_thr = PARAM.inp.force_thr;
    INPUT.stress_thr = PARAM.inp.stress_thr;
    INPUT.press1 = PARAM.inp.press1;
    INPUT.press2 = PARAM.inp.press2;
    INPUT.press3 = PARAM.inp.press3;
    INPUT.cal_stress = PARAM.inp.cal_stress;
    INPUT.nstream = PARAM.inp.nstream;
    INPUT.fixed_axes = PARAM.inp.fixed_axes;
    INPUT.fixed_ibrav = PARAM.inp.fixed_ibrav;
    INPUT.fixed_atoms = PARAM.inp.fixed_atoms;
    INPUT.relax_method = PARAM.inp.relax_method;
    INPUT.relax_new = PARAM.inp.relax_new;
    INPUT.relax_cg_thr = PARAM.inp.relax_cg_thr;
    INPUT.relax_scale_force = PARAM.inp.relax_scale_force;
    INPUT.gamma_only = PARAM.inp.gamma_only;
    INPUT.gamma_only_local = PARAM.inp.sup.gamma_only_local;
    INPUT.fft_mode = PARAM.inp.fft_mode;
    INPUT.ecutwfc = PARAM.inp.ecutwfc;
    INPUT.ecutrho = PARAM.inp.ecutrho;
    INPUT.erf_ecut = PARAM.inp.erf_ecut;
    INPUT.erf_height = PARAM.inp.erf_height;
    INPUT.erf_sigma = PARAM.inp.erf_sigma;
    INPUT.ncx = PARAM.inp.sup.ncx;
    INPUT.ncy = PARAM.inp.sup.ncy;
    INPUT.ncz = PARAM.inp.sup.ncz;
    INPUT.nx = PARAM.inp.nx;
    INPUT.ny = PARAM.inp.ny;
    INPUT.nz = PARAM.inp.nz;
    INPUT.bx = PARAM.inp.bx;
    INPUT.by = PARAM.inp.by;
    INPUT.bz = PARAM.inp.bz;
    INPUT.ndx = PARAM.inp.ndx;
    INPUT.ndy = PARAM.inp.ndy;
    INPUT.ndz = PARAM.inp.ndz;
    INPUT.diago_proc = PARAM.inp.diago_proc;
    INPUT.pw_diag_nmax = PARAM.inp.pw_diag_nmax;
    INPUT.diago_cg_prec = PARAM.inp.diago_cg_prec;
    INPUT.pw_diag_ndim = PARAM.inp.pw_diag_ndim;
    INPUT.pw_diag_thr = PARAM.inp.pw_diag_thr;
    INPUT.nb2d = PARAM.inp.nb2d;
    INPUT.nurse = PARAM.inp.nurse;
    INPUT.nbspline = PARAM.inp.nbspline;
    INPUT.colour = PARAM.inp.colour;
    INPUT.t_in_h = PARAM.inp.t_in_h;
    INPUT.vl_in_h = PARAM.inp.vl_in_h;
    INPUT.vnl_in_h = PARAM.inp.vnl_in_h;
    INPUT.vh_in_h = PARAM.inp.vh_in_h;
    INPUT.vion_in_h = PARAM.inp.vion_in_h;
    INPUT.test_force = PARAM.inp.test_force;
    INPUT.test_stress = PARAM.inp.test_stress;
    INPUT.scf_thr = PARAM.inp.scf_thr;
    INPUT.scf_thr_type = PARAM.inp.scf_thr_type;
    INPUT.scf_nmax = PARAM.inp.scf_nmax;
    INPUT.relax_nmax = PARAM.inp.relax_nmax;
    INPUT.out_level = PARAM.inp.out_level;
    INPUT.out_md_control = PARAM.inp.sup.out_md_control;
    INPUT.smearing_method = PARAM.inp.smearing_method;
    INPUT.smearing_sigma = PARAM.inp.smearing_sigma;

    INPUT.init_wfc = PARAM.inp.init_wfc;
    INPUT.init_chg = PARAM.inp.init_chg;
    INPUT.psi_initializer = PARAM.inp.psi_initializer;
    INPUT.chg_extrap = PARAM.inp.chg_extrap;
    INPUT.mem_saver = PARAM.inp.mem_saver;
    INPUT.printe = PARAM.inp.printe;
    INPUT.out_freq_elec = PARAM.inp.out_freq_elec;
    INPUT.out_freq_ion = PARAM.inp.out_freq_ion;
    INPUT.out_chg = PARAM.inp.out_chg;
    INPUT.out_pot = PARAM.inp.out_pot;
    INPUT.out_wfc_pw = PARAM.inp.out_wfc_pw;
    INPUT.out_wfc_r = PARAM.inp.out_wfc_r;
    INPUT.out_dos = PARAM.inp.out_dos;
    INPUT.out_band = PARAM.inp.out_band;
    INPUT.out_proj_band = PARAM.inp.out_proj_band;
    INPUT.cal_syns = PARAM.inp.cal_syns;
    INPUT.dmax = PARAM.inp.dmax;
    INPUT.out_interval = PARAM.inp.out_interval;
    INPUT.out_app_flag = PARAM.inp.out_app_flag;
    INPUT.out_ndigits = PARAM.inp.out_ndigits;
    INPUT.out_mat_r = PARAM.inp.out_mat_r;
    INPUT.out_alllog = PARAM.inp.out_alllog;
    INPUT.out_element_info = PARAM.inp.out_element_info;
    INPUT.out_bandgap = PARAM.inp.out_bandgap;
    INPUT.dos_emin_ev = PARAM.inp.dos_emin_ev;
    INPUT.dos_emax_ev = PARAM.inp.dos_emax_ev;
    INPUT.dos_edelta_ev = PARAM.inp.dos_edelta_ev;
    INPUT.dos_scale = PARAM.inp.dos_scale;
    INPUT.dos_nche = PARAM.inp.dos_nche;
    INPUT.dos_setemin = PARAM.inp.sup.dos_setemin;
    INPUT.dos_setemax = PARAM.inp.sup.dos_setemax;
    INPUT.dos_sigma = PARAM.inp.dos_sigma;
    INPUT.lcao_ecut = PARAM.inp.lcao_ecut;
    INPUT.lcao_dk = PARAM.inp.lcao_dk;
    INPUT.lcao_dr = PARAM.inp.lcao_dr;
    INPUT.lcao_rmax = PARAM.inp.lcao_rmax;
    INPUT.search_radius = PARAM.inp.search_radius;
    INPUT.search_pbc = PARAM.inp.search_pbc;
    INPUT.onsite_radius = PARAM.inp.onsite_radius;

    INPUT.mdp = PARAM.mdp;

    INPUT.vdw_method = PARAM.inp.vdw_method;
    INPUT.vdw_s6 = PARAM.inp.vdw_s6;
    INPUT.vdw_s8 = PARAM.inp.vdw_s8;
    INPUT.vdw_a1 = PARAM.inp.vdw_a1;
    INPUT.vdw_a2 = PARAM.inp.vdw_a2;
    INPUT.vdw_d = PARAM.inp.vdw_d;
    INPUT.vdw_abc = PARAM.inp.vdw_abc;
    INPUT.vdw_cutoff_radius = PARAM.inp.vdw_cutoff_radius;
    INPUT.vdw_radius_unit = PARAM.inp.vdw_radius_unit;
    INPUT.vdw_cn_thr = PARAM.inp.vdw_cn_thr;
    INPUT.vdw_cn_thr_unit = PARAM.inp.vdw_cn_thr_unit;
    INPUT.vdw_C6_file = PARAM.inp.vdw_C6_file;
    INPUT.vdw_C6_unit = PARAM.inp.vdw_C6_unit;
    INPUT.vdw_R0_file = PARAM.inp.vdw_R0_file;
    INPUT.vdw_R0_unit = PARAM.inp.vdw_R0_unit;
    INPUT.vdw_cutoff_type = PARAM.inp.vdw_cutoff_type;
    INPUT.vdw_cutoff_period = PARAM.inp.vdw_cutoff_period;

    INPUT.ocp = PARAM.inp.ocp;
    INPUT.ocp_set = PARAM.inp.ocp_set;
    INPUT.out_mul = PARAM.inp.out_mul;
    INPUT.noncolin = PARAM.inp.noncolin;
    INPUT.lspinorb = PARAM.inp.lspinorb;
    INPUT.soc_lambda = PARAM.inp.soc_lambda;

    INPUT.restart_save = PARAM.inp.restart_save;
    INPUT.restart_load = PARAM.inp.restart_load;
    INPUT.cell_factor = PARAM.inp.cell_factor;
    INPUT.dft_plus_u = PARAM.inp.dft_plus_u;
    delete[] INPUT.orbital_corr;
    INPUT.orbital_corr = new int[INPUT.ntype];
    for (int i = 0; i < INPUT.ntype; ++i)
    {
        INPUT.orbital_corr[i] = PARAM.inp.orbital_corr[i];
    }
    delete[] INPUT.hubbard_u;
    INPUT.hubbard_u = new double[INPUT.ntype];
    for (int i = 0; i < INPUT.ntype; ++i)
    {
        INPUT.hubbard_u[i] = PARAM.inp.sup.hubbard_u[i];
    }
    INPUT.omc = PARAM.inp.omc;
    INPUT.yukawa_potential = PARAM.inp.yukawa_potential;
    INPUT.yukawa_lambda = PARAM.inp.yukawa_lambda;
    INPUT.uramping = PARAM.inp.sup.uramping;
    INPUT.dft_plus_dmft = PARAM.inp.dft_plus_dmft;
    INPUT.rpa = PARAM.inp.rpa;
    INPUT.of_kinetic = PARAM.inp.of_kinetic;
    INPUT.of_method = PARAM.inp.of_method;
    INPUT.of_conv = PARAM.inp.of_conv;
    INPUT.of_tole = PARAM.inp.of_tole;
    INPUT.of_tolp = PARAM.inp.of_tolp;
    INPUT.of_tf_weight = PARAM.inp.of_tf_weight;
    INPUT.of_vw_weight = PARAM.inp.of_vw_weight;
    INPUT.of_wt_alpha = PARAM.inp.of_wt_alpha;
    INPUT.of_wt_beta = PARAM.inp.of_wt_beta;
    INPUT.of_wt_rho0 = PARAM.inp.of_wt_rho0;
    INPUT.of_hold_rho0 = PARAM.inp.of_hold_rho0;
    INPUT.of_lkt_a = PARAM.inp.of_lkt_a;
    INPUT.of_full_pw = PARAM.inp.of_full_pw;
    INPUT.of_full_pw_dim = PARAM.inp.of_full_pw_dim;
    INPUT.of_read_kernel = PARAM.inp.of_read_kernel;
    INPUT.of_kernel_file = PARAM.inp.of_kernel_file;
    INPUT.bessel_nao_smooth = PARAM.inp.bessel_nao_smooth;
    INPUT.bessel_nao_sigma = PARAM.inp.bessel_nao_sigma;
    INPUT.bessel_nao_ecut = PARAM.inp.bessel_nao_ecut;
    INPUT.bessel_nao_rcut = PARAM.inp.sup.bessel_nao_rcut;
    INPUT.bessel_nao_rcuts = PARAM.inp.bessel_nao_rcuts;
    INPUT.bessel_nao_tolerence = PARAM.inp.bessel_nao_tolerence;
    INPUT.bessel_descriptor_lmax = PARAM.inp.bessel_descriptor_lmax;
    INPUT.bessel_descriptor_smooth = PARAM.inp.bessel_descriptor_smooth;
    INPUT.bessel_descriptor_sigma = PARAM.inp.bessel_descriptor_sigma;
    INPUT.bessel_descriptor_ecut = PARAM.inp.bessel_descriptor_ecut;
    INPUT.bessel_descriptor_rcut = PARAM.inp.bessel_descriptor_rcut;
    INPUT.bessel_descriptor_tolerence = PARAM.inp.bessel_descriptor_tolerence;
    INPUT.device = PARAM.inp.device;
    INPUT.precision = PARAM.inp.precision;
    INPUT.test_skip_ewald = PARAM.inp.test_skip_ewald;
    INPUT.use_paw = PARAM.inp.use_paw;
}

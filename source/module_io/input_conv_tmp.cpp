#include <sstream>
#include <string>

#include "input.h"
#include "input_conv.h"
#include "module_parameter/parameter.h"

void Input_Conv::tmp_convert()
{
    INPUT.stru_file = PARAM.get().stru_file;
    INPUT.pseudo_dir = PARAM.get().pseudo_dir;
    INPUT.orbital_dir = PARAM.get().orbital_dir;
    INPUT.read_file_dir = PARAM.get().read_file_dir;
    INPUT.kpoint_file = PARAM.get().kpoint_file;
    INPUT.wannier_card = PARAM.get().wannier_card;
    INPUT.latname = PARAM.get().latname;
    INPUT.calculation = PARAM.get().calculation;
    INPUT.esolver_type = PARAM.get().esolver_type;
    INPUT.pseudo_rcut = PARAM.get().pseudo_rcut;
    INPUT.pseudo_mesh = PARAM.get().pseudo_mesh;
    INPUT.ntype = PARAM.get().ntype;
    INPUT.nbands = PARAM.get().nbands;
    INPUT.nbands_istate = PARAM.get().nbands_istate;
    INPUT.pw_seed = PARAM.get().pw_seed;
    INPUT.init_vel = PARAM.get().init_vel;
    INPUT.ref_cell_factor = PARAM.get().ref_cell_factor;

    INPUT.kpar = PARAM.get().kpar;
    INPUT.berry_phase = PARAM.get().berry_phase;
    INPUT.gdir = PARAM.get().gdir;
    INPUT.kspacing[0] = PARAM.get().kspacing[0];
    INPUT.kspacing[1] = PARAM.get().kspacing[1];
    INPUT.kspacing[2] = PARAM.get().kspacing[2];
    INPUT.min_dist_coef = PARAM.get().min_dist_coef;
    INPUT.towannier90 = PARAM.get().towannier90;
    INPUT.nnkpfile = PARAM.get().nnkpfile;
    INPUT.wannier_spin = PARAM.get().wannier_spin;
    INPUT.wannier_method = PARAM.get().wannier_method;
    INPUT.out_wannier_mmn = PARAM.get().out_wannier_mmn;
    INPUT.out_wannier_amn = PARAM.get().out_wannier_amn;
    INPUT.out_wannier_unk = PARAM.get().out_wannier_unk;
    INPUT.out_wannier_eig = PARAM.get().out_wannier_eig;
    INPUT.out_wannier_wvfn_formatted = PARAM.get().out_wannier_wvfn_formatted;

    INPUT.nche_sto = PARAM.get().nche_sto;
    INPUT.nbands_sto = PARAM.get().nbands_sto;
    INPUT.seed_sto = PARAM.get().seed_sto;
    INPUT.initsto_ecut = PARAM.get().initsto_ecut;
    INPUT.emax_sto = PARAM.get().emax_sto;
    INPUT.emin_sto = PARAM.get().emin_sto;
    INPUT.bndpar = PARAM.get().bndpar;
    INPUT.initsto_freq = PARAM.get().initsto_freq;
    INPUT.method_sto = PARAM.get().method_sto;
    INPUT.npart_sto = PARAM.get().npart_sto;
    INPUT.cal_cond = PARAM.get().cal_cond;
    INPUT.cond_che_thr = PARAM.get().cond_che_thr;
    INPUT.cond_smear = PARAM.get().cond_smear;
    INPUT.cond_dw = PARAM.get().cond_dw;
    INPUT.cond_wcut = PARAM.get().cond_wcut;
    INPUT.cond_dt = PARAM.get().cond_dt;
    INPUT.cond_dtbatch = PARAM.get().cond_dtbatch;
    INPUT.cond_fwhm = PARAM.get().cond_fwhm;
    INPUT.cond_nonlocal = PARAM.get().cond_nonlocal;

    INPUT.dft_functional = PARAM.get().dft_functional;
    INPUT.xc_temperature = PARAM.get().xc_temperature;
    INPUT.nspin = PARAM.get().nspin;
    INPUT.lmaxmax = PARAM.get().lmaxmax;

    INPUT.basis_type = PARAM.get().basis_type;
    INPUT.ks_solver = PARAM.get().ks_solver;
    INPUT.cal_force = PARAM.get().cal_force;
    INPUT.force_thr = PARAM.get().force_thr;
    INPUT.stress_thr = PARAM.get().stress_thr;
    INPUT.press1 = PARAM.get().press1;
    INPUT.press2 = PARAM.get().press2;
    INPUT.press3 = PARAM.get().press3;
    INPUT.cal_stress = PARAM.get().cal_stress;
    INPUT.nstream = PARAM.get().nstream;
    INPUT.fixed_axes = PARAM.get().fixed_axes;
    INPUT.fixed_ibrav = PARAM.get().fixed_ibrav;
    INPUT.fixed_atoms = PARAM.get().fixed_atoms;
    INPUT.relax_method = PARAM.get().relax_method;
    INPUT.relax_new = PARAM.get().relax_new;
    INPUT.relax_cg_thr = PARAM.get().relax_cg_thr;
    INPUT.relax_scale_force = PARAM.get().relax_scale_force;
    INPUT.gamma_only = PARAM.get().gamma_only;
    INPUT.gamma_only_local = PARAM.get().sup.gamma_only_local;
    INPUT.fft_mode = PARAM.get().fft_mode;
    INPUT.ecutwfc = PARAM.get().ecutwfc;
    INPUT.ecutrho = PARAM.get().ecutrho;
    INPUT.erf_ecut = PARAM.get().erf_ecut;
    INPUT.erf_height = PARAM.get().erf_height;
    INPUT.erf_sigma = PARAM.get().erf_sigma;
    INPUT.ncx = PARAM.get().sup.ncx;
    INPUT.ncy = PARAM.get().sup.ncy;
    INPUT.ncz = PARAM.get().sup.ncz;
    INPUT.nx = PARAM.get().nx;
    INPUT.ny = PARAM.get().ny;
    INPUT.nz = PARAM.get().nz;
    INPUT.bx = PARAM.get().bx;
    INPUT.by = PARAM.get().by;
    INPUT.bz = PARAM.get().bz;
    INPUT.ndx = PARAM.get().ndx;
    INPUT.ndy = PARAM.get().ndy;
    INPUT.ndz = PARAM.get().ndz;
    INPUT.diago_proc = PARAM.get().diago_proc;
    INPUT.pw_diag_nmax = PARAM.get().pw_diag_nmax;
    INPUT.diago_cg_prec = PARAM.get().diago_cg_prec;
    INPUT.pw_diag_ndim = PARAM.get().pw_diag_ndim;
    INPUT.pw_diag_thr = PARAM.get().pw_diag_thr;
    INPUT.nb2d = PARAM.get().nb2d;
    INPUT.nurse = PARAM.get().nurse;
    INPUT.nbspline = PARAM.get().nbspline;
    INPUT.colour = PARAM.get().colour;
    INPUT.t_in_h = PARAM.get().t_in_h;
    INPUT.vl_in_h = PARAM.get().vl_in_h;
    INPUT.vnl_in_h = PARAM.get().vnl_in_h;
    INPUT.vh_in_h = PARAM.get().vh_in_h;
    INPUT.vion_in_h = PARAM.get().vion_in_h;
    INPUT.test_force = PARAM.get().test_force;
    INPUT.test_stress = PARAM.get().test_stress;
    INPUT.scf_thr = PARAM.get().scf_thr;
    INPUT.scf_thr_type = PARAM.get().scf_thr_type;
    INPUT.scf_nmax = PARAM.get().scf_nmax;
    INPUT.relax_nmax = PARAM.get().relax_nmax;
    INPUT.out_level = PARAM.get().out_level;
    INPUT.out_md_control = PARAM.get().sup.out_md_control;
    INPUT.smearing_method = PARAM.get().smearing_method;
    INPUT.smearing_sigma = PARAM.get().smearing_sigma;

    INPUT.init_wfc = PARAM.get().init_wfc;
    INPUT.init_chg = PARAM.get().init_chg;
    INPUT.psi_initializer = PARAM.get().psi_initializer;
    INPUT.chg_extrap = PARAM.get().chg_extrap;
    INPUT.mem_saver = PARAM.get().mem_saver;
    INPUT.printe = PARAM.get().printe;
    INPUT.out_freq_elec = PARAM.get().out_freq_elec;
    INPUT.out_freq_ion = PARAM.get().out_freq_ion;
    INPUT.out_chg = PARAM.get().out_chg;
    INPUT.out_pot = PARAM.get().out_pot;
    INPUT.out_wfc_pw = PARAM.get().out_wfc_pw;
    INPUT.out_wfc_r = PARAM.get().out_wfc_r;
    INPUT.out_dos = PARAM.get().out_dos;
    INPUT.out_band = PARAM.get().out_band;
    INPUT.out_proj_band = PARAM.get().out_proj_band;
    INPUT.cal_syns = PARAM.get().cal_syns;
    INPUT.dmax = PARAM.get().dmax;
    INPUT.out_interval = PARAM.get().out_interval;
    INPUT.out_app_flag = PARAM.get().out_app_flag;
    INPUT.out_ndigits = PARAM.get().out_ndigits;
    INPUT.out_mat_r = PARAM.get().out_mat_r;
    INPUT.out_alllog = PARAM.get().out_alllog;
    INPUT.out_element_info = PARAM.get().out_element_info;
    INPUT.out_bandgap = PARAM.get().out_bandgap;
    INPUT.dos_emin_ev = PARAM.get().dos_emin_ev;
    INPUT.dos_emax_ev = PARAM.get().dos_emax_ev;
    INPUT.dos_edelta_ev = PARAM.get().dos_edelta_ev;
    INPUT.dos_scale = PARAM.get().dos_scale;
    INPUT.dos_nche = PARAM.get().dos_nche;
    INPUT.dos_setemin = PARAM.get().sup.dos_setemin;
    INPUT.dos_setemax = PARAM.get().sup.dos_setemax;
    INPUT.dos_sigma = PARAM.get().dos_sigma;
    INPUT.lcao_ecut = PARAM.get().lcao_ecut;
    INPUT.lcao_dk = PARAM.get().lcao_dk;
    INPUT.lcao_dr = PARAM.get().lcao_dr;
    INPUT.lcao_rmax = PARAM.get().lcao_rmax;
    INPUT.search_radius = PARAM.get().search_radius;
    INPUT.search_pbc = PARAM.get().search_pbc;
    INPUT.onsite_radius = PARAM.get().onsite_radius;

    INPUT.mdp = PARAM.get_mdp();

    INPUT.vdw_method = PARAM.get().vdw_method;
    INPUT.vdw_s6 = PARAM.get().vdw_s6;
    INPUT.vdw_s8 = PARAM.get().vdw_s8;
    INPUT.vdw_a1 = PARAM.get().vdw_a1;
    INPUT.vdw_a2 = PARAM.get().vdw_a2;
    INPUT.vdw_d = PARAM.get().vdw_d;
    INPUT.vdw_abc = PARAM.get().vdw_abc;
    INPUT.vdw_cutoff_radius = PARAM.get().vdw_cutoff_radius;
    INPUT.vdw_radius_unit = PARAM.get().vdw_radius_unit;
    INPUT.vdw_cn_thr = PARAM.get().vdw_cn_thr;
    INPUT.vdw_cn_thr_unit = PARAM.get().vdw_cn_thr_unit;
    INPUT.vdw_C6_file = PARAM.get().vdw_C6_file;
    INPUT.vdw_C6_unit = PARAM.get().vdw_C6_unit;
    INPUT.vdw_R0_file = PARAM.get().vdw_R0_file;
    INPUT.vdw_R0_unit = PARAM.get().vdw_R0_unit;
    INPUT.vdw_cutoff_type = PARAM.get().vdw_cutoff_type;
    INPUT.vdw_cutoff_period = PARAM.get().vdw_cutoff_period;

    INPUT.ocp = PARAM.get().ocp;
    INPUT.ocp_set = PARAM.get().ocp_set;
    INPUT.out_mul = PARAM.get().out_mul;
    INPUT.noncolin = PARAM.get().noncolin;
    INPUT.lspinorb = PARAM.get().lspinorb;
    INPUT.soc_lambda = PARAM.get().soc_lambda;

    INPUT.restart_save = PARAM.get().restart_save;
    INPUT.restart_load = PARAM.get().restart_load;
    INPUT.cell_factor = PARAM.get().cell_factor;
    INPUT.dft_plus_u = PARAM.get().dft_plus_u;
    delete[] INPUT.orbital_corr;
    INPUT.orbital_corr = new int[INPUT.ntype];
    for (int i = 0; i < INPUT.ntype; ++i)
    {
        INPUT.orbital_corr[i] = PARAM.get().orbital_corr[i];
    }
    delete[] INPUT.hubbard_u;
    INPUT.hubbard_u = new double[INPUT.ntype];
    for (int i = 0; i < INPUT.ntype; ++i)
    {
        INPUT.hubbard_u[i] = PARAM.get().sup.hubbard_u[i];
    }
    INPUT.omc = PARAM.get().omc;
    INPUT.yukawa_potential = PARAM.get().yukawa_potential;
    INPUT.yukawa_lambda = PARAM.get().yukawa_lambda;
    INPUT.uramping = PARAM.get().sup.uramping;
    INPUT.dft_plus_dmft = PARAM.get().dft_plus_dmft;
    INPUT.rpa = PARAM.get().rpa;
    INPUT.of_kinetic = PARAM.get().of_kinetic;
    INPUT.of_method = PARAM.get().of_method;
    INPUT.of_conv = PARAM.get().of_conv;
    INPUT.of_tole = PARAM.get().of_tole;
    INPUT.of_tolp = PARAM.get().of_tolp;
    INPUT.of_tf_weight = PARAM.get().of_tf_weight;
    INPUT.of_vw_weight = PARAM.get().of_vw_weight;
    INPUT.of_wt_alpha = PARAM.get().of_wt_alpha;
    INPUT.of_wt_beta = PARAM.get().of_wt_beta;
    INPUT.of_wt_rho0 = PARAM.get().of_wt_rho0;
    INPUT.of_hold_rho0 = PARAM.get().of_hold_rho0;
    INPUT.of_lkt_a = PARAM.get().of_lkt_a;
    INPUT.of_full_pw = PARAM.get().of_full_pw;
    INPUT.of_full_pw_dim = PARAM.get().of_full_pw_dim;
    INPUT.of_read_kernel = PARAM.get().of_read_kernel;
    INPUT.of_kernel_file = PARAM.get().of_kernel_file;
    INPUT.bessel_nao_smooth = PARAM.get().bessel_nao_smooth;
    INPUT.bessel_nao_sigma = PARAM.get().bessel_nao_sigma;
    INPUT.bessel_nao_ecut = PARAM.get().bessel_nao_ecut;
    INPUT.bessel_nao_rcut = PARAM.get().sup.bessel_nao_rcut;
    INPUT.bessel_nao_rcuts = PARAM.get().bessel_nao_rcuts;
    INPUT.bessel_nao_tolerence = PARAM.get().bessel_nao_tolerence;
    INPUT.bessel_descriptor_lmax = PARAM.get().bessel_descriptor_lmax;
    INPUT.bessel_descriptor_smooth = PARAM.get().bessel_descriptor_smooth;
    INPUT.bessel_descriptor_sigma = PARAM.get().bessel_descriptor_sigma;
    INPUT.bessel_descriptor_ecut = PARAM.get().bessel_descriptor_ecut;
    INPUT.bessel_descriptor_rcut = PARAM.get().bessel_descriptor_rcut;
    INPUT.bessel_descriptor_tolerence = PARAM.get().bessel_descriptor_tolerence;
    INPUT.device = PARAM.get().device;
    INPUT.precision = PARAM.get().precision;
    INPUT.test_skip_ewald = PARAM.get().test_skip_ewald;
    INPUT.use_paw = PARAM.get().use_paw;
}

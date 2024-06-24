#include "parameter.h"

Parameter PARAM;

Parameter::Parameter()
{
}
Parameter::~Parameter()
{
}

std::string Parameter::get_suffix() const
{
    return suffix;
}
std::string Parameter::get_latname() const
{
    return latname;
}
std::string Parameter::get_stru_file() const
{
    return stru_file;
}
std::string Parameter::get_kpoint_file() const
{
    return kpoint_file;
}
std::string Parameter::get_pseudo_dir() const
{
    return pseudo_dir;
}
std::string Parameter::get_orbital_dir() const
{
    return orbital_dir;
}
double Parameter::get_pseudo_rcut() const
{
    return pseudo_rcut;
}
bool Parameter::get_pseudo_mesh() const
{
    return pseudo_mesh;
}
int Parameter::get_lmaxmax() const
{
    return lmaxmax;
}
std::string Parameter::get_dft_functional() const
{
    return dft_functional;
}
double Parameter::get_xc_temperature() const
{
    return xc_temperature;
}
std::string Parameter::get_calculation() const
{
    return calculation;
}
std::string Parameter::get_esolver_type() const
{
    return esolver_type;
}
int Parameter::get_ntype() const
{
    return ntype;
}
int Parameter::get_nspin() const
{
    return nspin;
}
std::vector<double> Parameter::get_kspacing() const
{
    return kspacing;
}
double Parameter::get_min_dist_coef() const
{
    return min_dist_coef;
}
int Parameter::get_nbands() const
{
    return nbands;
}
int Parameter::get_nbands_istate() const
{
    return nbands_istate;
}
std::string Parameter::get_bands_to_print() const
{
    return bands_to_print;
}
std::string Parameter::get_symmetry() const
{
    return symmetry;
}
bool Parameter::get_init_vel() const
{
    return init_vel;
}
double Parameter::get_symmetry_prec() const
{
    return symmetry_prec;
}
bool Parameter::get_symmetry_autoclose() const
{
    return symmetry_autoclose;
}
double Parameter::get_nelec() const
{
    return nelec;
}
double Parameter::get_nelec_delta() const
{
    return nelec_delta;
}
double Parameter::get_nupdown() const
{
    return nupdown;
}
bool Parameter::get_two_fermi() const
{
    return two_fermi;
}
bool Parameter::get_out_mul() const
{
    return out_mul;
}
bool Parameter::get_noncolin() const
{
    return noncolin;
}
bool Parameter::get_lspinorb() const
{
    return lspinorb;
}
int Parameter::get_kpar() const
{
    return kpar;
}
int Parameter::get_bndpar() const
{
    return bndpar;
}
int Parameter::get_out_freq_elec() const
{
    return out_freq_elec;
}
bool Parameter::get_dft_plus_dmft() const
{
    return dft_plus_dmft;
}
bool Parameter::get_rpa() const
{
    return rpa;
}
int Parameter::get_printe() const
{
    return printe;
}
int Parameter::get_mem_saver() const
{
    return mem_saver;
}
int Parameter::get_diago_proc() const
{
    return diago_proc;
}
int Parameter::get_nbspline() const
{
    return nbspline;
}
std::string Parameter::get_wannier_card() const
{
    return wannier_card;
}
double Parameter::get_soc_lambda() const
{
    return soc_lambda;
}
bool Parameter::get_cal_force() const
{
    return cal_force;
}
int Parameter::get_out_freq_ion() const
{
    return out_freq_ion;
}
int Parameter::get_elpa_num_thread() const
{
    return elpa_num_thread;
}
std::string Parameter::get_device() const
{
    return device;
}
std::string Parameter::get_precision() const
{
    return precision;
}

double Parameter::get_ecutwfc() const
{
    return ecutwfc;
}
double Parameter::get_ecutrho() const
{
    return ecutrho;
}
double Parameter::get_erf_ecut() const
{
    return erf_ecut;
}
double Parameter::get_erf_height() const
{
    return erf_height;
}
double Parameter::get_erf_sigma() const
{
    return erf_sigma;
}
int Parameter::get_fft_mode() const
{
    return fft_mode;
}
int Parameter::get_pw_diag_nmax() const
{
    return pw_diag_nmax;
}
int Parameter::get_diago_cg_prec() const
{
    return diago_cg_prec;
}
int Parameter::get_pw_diag_ndim() const
{
    return pw_diag_ndim;
}
bool Parameter::get_diago_full_acc() const
{
    return diago_full_acc;
}
double Parameter::get_pw_diag_thr() const
{
    return pw_diag_thr;
}
int Parameter::get_nb2d() const
{
    return nb2d;
}
double Parameter::get_scf_thr() const
{
    return scf_thr;
}
int Parameter::get_scf_thr_type() const
{
    return scf_thr_type;
}
std::string Parameter::get_init_wfc() const
{
    return init_wfc;
}
bool Parameter::get_psi_initializer() const
{
    return psi_initializer;
}
std::string Parameter::get_init_chg() const
{
    return init_chg;
}
std::string Parameter::get_chg_extrap() const
{
    return chg_extrap;
}
int Parameter::get_out_chg() const
{
    return out_chg;
}
int Parameter::get_out_pot() const
{
    return out_pot;
}
int Parameter::get_out_wfc_pw() const
{
    return out_wfc_pw;
}
bool Parameter::get_out_wfc_r() const
{
    return out_wfc_r;
}
int Parameter::get_out_dos() const
{
    return out_dos;
}
std::vector<int> Parameter::get_out_band() const
{
    return out_band;
}
bool Parameter::get_out_proj_band() const
{
    return out_proj_band;
}
bool Parameter::get_restart_save() const
{
    return restart_save;
}
bool Parameter::get_restart_load() const
{
    return restart_load;
}
std::string Parameter::get_read_file_dir() const
{
    return read_file_dir;
}
int Parameter::get_nx() const
{
    return nx;
}
int Parameter::get_ny() const
{
    return ny;
}
int Parameter::get_nz() const
{
    return nz;
}
int Parameter::get_ncx() const
{
    return ncx;
}
int Parameter::get_ncy() const
{
    return ncy;
}
int Parameter::get_ncz() const
{
    return ncz;
}
int Parameter::get_ndx() const
{
    return ndx;
}
int Parameter::get_ndy() const
{
    return ndy;
}
int Parameter::get_ndz() const
{
    return ndz;
}
double Parameter::get_cell_factor() const
{
    return cell_factor;
}
int Parameter::get_pw_seed() const
{
    return pw_seed;
}

int Parameter::get_method_sto() const
{
    return method_sto;
}
int Parameter::get_npart_sto() const
{
    return npart_sto;
}
int Parameter::get_nbands_sto() const
{
    return nbands_sto;
}
std::string Parameter::get_nbndsto_str() const
{
    return nbndsto_str;
}
int Parameter::get_nche_sto() const
{
    return nche_sto;
}
double Parameter::get_emin_sto() const
{
    return emin_sto;
}
double Parameter::get_emax_sto() const
{
    return emax_sto;
}
int Parameter::get_seed_sto() const
{
    return seed_sto;
}
double Parameter::get_initsto_ecut() const
{
    return initsto_ecut;
}
int Parameter::get_initsto_freq() const
{
    return initsto_freq;
}
bool Parameter::get_cal_cond() const
{
    return cal_cond;
}
double Parameter::get_cond_che_thr() const
{
    return cond_che_thr;
}
double Parameter::get_cond_dw() const
{
    return cond_dw;
}
double Parameter::get_cond_wcut() const
{
    return cond_wcut;
}
double Parameter::get_cond_dt() const
{
    return cond_dt;
}
int Parameter::get_cond_dtbatch() const
{
    return cond_dtbatch;
}
int Parameter::get_cond_smear() const
{
    return cond_smear;
}
double Parameter::get_cond_fwhm() const
{
    return cond_fwhm;
}
bool Parameter::get_cond_nonlocal() const
{
    return cond_nonlocal;
}

// ==============   #Parameters (4.Relaxation) ===========================
std::string Parameter::get_ks_solver() const
{
    return ks_solver;
}
int Parameter::get_scf_nmax() const
{
    return scf_nmax;
}
int Parameter::get_relax_nmax() const
{
    return relax_nmax;
}
bool Parameter::get_out_stru() const
{
    return out_stru;
}
double Parameter::get_force_thr() const
{
    return force_thr;
}
double Parameter::get_force_thr_ev() const
{
    return force_thr_ev;
}
double Parameter::get_force_thr_ev2() const
{
    return force_thr_ev2;
}
double Parameter::get_relax_cg_thr() const
{
    return relax_cg_thr;
}
double Parameter::get_stress_thr() const
{
    return stress_thr;
}
double Parameter::get_press1() const
{
    return press1;
}
double Parameter::get_press2() const
{
    return press2;
}
double Parameter::get_press3() const
{
    return press3;
}
double Parameter::get_relax_bfgs_w1() const
{
    return relax_bfgs_w1;
}
double Parameter::get_relax_bfgs_w2() const
{
    return relax_bfgs_w2;
}
double Parameter::get_relax_bfgs_rmax() const
{
    return relax_bfgs_rmax;
}
double Parameter::get_relax_bfgs_rmin() const
{
    return relax_bfgs_rmin;
}
double Parameter::get_relax_bfgs_init() const
{
    return relax_bfgs_init;
}
bool Parameter::get_cal_stress() const
{
    return cal_stress;
}
std::string Parameter::get_fixed_axes() const
{
    return fixed_axes;
}
bool Parameter::get_fixed_ibrav() const
{
    return fixed_ibrav;
}
bool Parameter::get_fixed_atoms() const
{
    return fixed_atoms;
}
std::string Parameter::get_relax_method() const
{
    return relax_method;
}
bool Parameter::get_relax_new() const
{
    return relax_new;
}
double Parameter::get_relax_scale_force() const
{
    return relax_scale_force;
}
std::string Parameter::get_out_level() const
{
    return out_level;
}
bool Parameter::get_out_md_control() const
{
    return out_md_control;
}
bool Parameter::get_out_dm() const
{
    return out_dm;
}
bool Parameter::get_out_dm1() const
{
    return out_dm1;
}
bool Parameter::get_out_bandgap() const
{
    return out_bandgap;
}
bool Parameter::get_use_paw() const
{
    return use_paw;
}

bool Parameter::get_deepks_out_labels() const
{
    return deepks_out_labels;
}
bool Parameter::get_deepks_scf() const
{
    return deepks_scf;
}
bool Parameter::get_deepks_bandgap() const
{
    return deepks_bandgap;
}
bool Parameter::get_deepks_equiv() const
{
    return deepks_equiv;
}
bool Parameter::get_deepks_out_unittest() const
{
    return deepks_out_unittest;
}
std::string Parameter::get_deepks_model() const
{
    return deepks_model;
}

// ==============   #Parameters (5.LCAO) ===========================
std::string Parameter::get_basis_type() const
{
    return basis_type;
}
bool Parameter::get_gamma_only() const
{
    return gamma_only;
}
bool Parameter::get_gamma_only_local() const
{
    return gamma_only_local;
}

double Parameter::get_search_radius() const
{
    return search_radius;
}
bool Parameter::get_search_pbc() const
{
    return search_pbc;
}
double Parameter::get_lcao_ecut() const
{
    return lcao_ecut;
}
double Parameter::get_lcao_dk() const
{
    return lcao_dk;
}
double Parameter::get_lcao_dr() const
{
    return lcao_dr;
}
double Parameter::get_lcao_rmax() const
{
    return lcao_rmax;
}

std::vector<int> Parameter::get_out_mat_hs() const
{
    return out_mat_hs;
}
bool Parameter::get_out_mat_hs2() const
{
    return out_mat_hs2;
}
bool Parameter::get_out_mat_dh() const
{
    return out_mat_dh;
}
bool Parameter::get_out_mat_xc() const
{
    return out_mat_xc;
}
bool Parameter::get_out_hr_npz() const
{
    return out_hr_npz;
}
bool Parameter::get_out_dm_npz() const
{
    return out_dm_npz;
}
bool Parameter::get_dm_to_rho() const
{
    return dm_to_rho;
}
int Parameter::get_out_interval() const
{
    return out_interval;
}
bool Parameter::get_out_app_flag() const
{
    return out_app_flag;
}
int Parameter::get_out_ndigits() const
{
    return out_ndigits;
}
bool Parameter::get_out_mat_t() const
{
    return out_mat_t;
}
bool Parameter::get_out_element_info() const
{
    return out_element_info;
}
bool Parameter::get_out_mat_r() const
{
    return out_mat_r;
}
int Parameter::get_out_wfc_lcao() const
{
    return out_wfc_lcao;
}
int Parameter::get_bx() const
{
    return bx;
}
int Parameter::get_by() const
{
    return by;
}
int Parameter::get_bz() const
{
    return bz;
}
int Parameter::get_nstream() const
{
    return nstream;
}

// ==============   #Parameters (6.Smearing) ===========================

std::string Parameter::get_smearing_method() const
{
    return smearing_method;
}
double Parameter::get_smearing_sigma() const
{
    return smearing_sigma;
}
// ==============   #Parameters (7.Charge Mixing) ======================
std::string Parameter::get_mixing_mode() const
{
    return mixing_mode;
}
double Parameter::get_mixing_beta() const
{
    return mixing_beta;
}
int Parameter::get_mixing_ndim() const
{
    return mixing_ndim;
}
double Parameter::get_mixing_restart() const
{
    return mixing_restart;
}
double Parameter::get_mixing_gg0() const
{
    return mixing_gg0;
}
double Parameter::get_mixing_beta_mag() const
{
    return mixing_beta_mag;
}
double Parameter::get_mixing_gg0_mag() const
{
    return mixing_gg0_mag;
}
double Parameter::get_mixing_gg0_min() const
{
    return mixing_gg0_min;
}
double Parameter::get_mixing_angle() const
{
    return mixing_angle;
}
bool is_mixing_tau() const
{
    return mixing_tau;
}
bool is_mixing_dftu() const
{
    return mixing_dftu;
}
bool is_mixing_dmr() const
{
    return mixing_dmr;
}
// ==============   #Parameters (8.DOS) ===============================
double Parameter::get_dos_emin_ev() const
{
    return dos_emin_ev;
}
double Parameter::get_dos_emax_ev() const
{
    return dos_emax_ev;
}
bool Parameter::get_dos_setemin() const
{
    return dos_setemin;
}
bool Parameter::get_dos_setemax() const
{
    return dos_setemax;
}
double Parameter::get_dos_edelta_ev() const
{
    return dos_edelta_ev;
}
double Parameter::get_dos_scale() const
{
    return dos_scale;
}
double Parameter::get_dos_sigma() const
{
    return dos_sigma;
}
int Parameter::get_dos_nche() const
{
    return dos_nche;
}

// ==============   #Parameters (9.Molecular dynamics) ================
MD_para Parameter::get_mdp() const
{
    return mdp;
}
double Parameter::get_ref_cell_factor() const
{
    return ref_cell_factor;
}
bool Parameter::get_cal_syns() const
{
    return cal_syns;
}
double Parameter::get_dmax() const
{
    return dmax;
}

// =======   #Parameters (10.Electric field and dipole correction) ====
bool Parameter::get_efield_flag() const
{
    return efield_flag;
}
bool Parameter::get_dip_cor_flag() const
{
    return dip_cor_flag;
}
int Parameter::get_efield_dir() const
{
    return efield_dir;
}
double Parameter::get_efield_pos_max() const
{
    return efield_pos_max;
}
double Parameter::get_efield_pos_dec() const
{
    return efield_pos_dec;
}
double Parameter::get_efield_amp() const
{
    return efield_amp;
}

bool Parameter::get_gate_flag() const
{
    return gate_flag;
}
double Parameter::get_zgate() const
{
    return zgate;
}
bool Parameter::get_relax() const
{
    return relax;
}
bool Parameter::get_block() const
{
    return block;
}
double Parameter::get_block_down() const
{
    return block_down;
}
double Parameter::get_block_up() const
{
    return block_up;
}
double Parameter::get_block_height() const
{
    return block_height;
}
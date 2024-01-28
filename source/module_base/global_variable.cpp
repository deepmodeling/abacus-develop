//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#include "global_variable.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
namespace GlobalV
{

//----------------------------------------------------------
// EXPLAIN : Basic Global Variables
// In practice calculation, these values are set in
// input.cpp.
//----------------------------------------------------------
int NBANDS = 0;
int NBANDS_ISTATE = 0; // default number.
int NLOCAL = 0; // total number of local basis.

double  KSPACING[3] = {0.0,0.0,0.0};
double MIN_DIST_COEF = 0.2;

double PSEUDORCUT;
bool PSEUDO_MESH;

std::string CALCULATION = "scf";
std::string ESOLVER_TYPE = "ksdft";
int EFIELD_FLAG = 0; // 5: add electric field
int DIP_COR_FLAG = 0; // 7: add dipole field
bool GATE_FLAG = false;    // add gate field
bool out_app_flag = true;  // whether output r(R), H(R), S(R), T(R), and dH(R) matrices in an append manner during MD  liuyu 2023-03-20
int out_ndigits = 8;
std::string DFT_FUNCTIONAL = "default";
double XC_TEMPERATURE = 0.0;
int NSPIN = 1; // LDA
bool TWO_EFERMI = 0; // two fermi energy, exist only magnetization is fixed.
double nupdown = 0.0;
int CURRENT_SPIN = 0;
int CURRENT_K = 0;
int CAL_FORCE = 0; // if cal_force >1, means do the grid integration 'cal_force' times.
double FORCE_THR = 1.0e-3;
bool CAL_STRESS = false;
double PRESS1 = 0.0;
double PRESS2 = 0.0;
double PRESS3 = 0.0;
double PRESSURE = 0.0;
std::string RELAX_METHOD = "bfgs";
std::string OUT_LEVEL = "ie";
double relax_scale_force = 0.5;
bool relax_new = true;
bool fixed_atoms = false;
int OUT_FREQ_ELEC = 0;
int OUT_FREQ_ION = 0;
int RELAX_NMAX = 20;
int md_prec_level = 0;
int SCF_NMAX = 100;

bool use_paw = false;
bool use_uspp = false;
bool double_grid = false;

std::string BASIS_TYPE = "pw"; // xiaohui add 2013-09-01
std::string KS_SOLVER = "cg"; // xiaohui add 2013-09-01
double SEARCH_RADIUS = -1.0;
bool SEARCH_PBC = true;
bool SPARSE_MATRIX = false;

int DIAGO_PROC = 0;
int PW_DIAG_NMAX = 30;
int DIAGO_CG_PREC = 1; // mohan add 2012-03-31
int PW_DIAG_NDIM = 4;
double PW_DIAG_THR = 1.0e-2;
int NB2D = 1;

double SCF_THR = 1.0e-9;
int SCF_THR_TYPE = 1;

double DQ = 0.010; // space between Q points of the reciprocal radial tab
int NQX = 10000; // number of points describing reciprocal radial tab
int NQXQ = 10000; // number of points describing reciprocal radial tab for Q

int NURSE = 0; // used for debug.
bool COLOUR = 0;
bool GAMMA_ONLY_LOCAL = 0; // mohan add 2010-10-20
bool GAMMA_ONLY_PW = 0; // mohan add 2012-06-05

int T_IN_H = 1; // mohan add 2010-11-28
int VL_IN_H = 1;
int VNL_IN_H = 1;
int VH_IN_H = 1;
int VION_IN_H = 1;
int ZEEMAN_IN_H = 1;
double STRESS_THR = 0.5; // LiuXh add 20180515 liuyu update 2023-05-10

int ocp = 0;
std::string ocp_set = "none";
std::vector<double> ocp_kb;
// int ocp_n=0;
// double ocp_kb[10000];
bool out_mul = false; // qifeng add 2019/9/10
//----------------------------------------------------------
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : NPROC( global number of process )
// NAME : KPAR( global number of pools )
// NAME : MY_RANK( global index of process )
// NAME : MY_POOL( global index of pool (count in pool))
// NAME : NPROC_IN_POOL( local number of process in a pool.)
// NAME : RANK_IN_POOL( global index of pool (count in process),
//  		my_rank in each pool)
//----------------------------------------------------------
int NPROC = 1;
int KPAR = 1;
int NSTOGROUP = 1;
int MY_RANK = 0;
int MY_POOL = 0;
int MY_STOGROUP = 0;
int NPROC_IN_POOL = 1;
int NPROC_IN_STOGROUP = 1;
int RANK_IN_POOL = 0;
int RANK_IN_STOGROUP = 0;
int DRANK = -1; // mohan add 2012-01-13, must be -1, so we can recognize who didn't in DIAG_WORLD
int DSIZE = KPAR;
int DCOLOR = -1;
int GRANK = MY_RANK;
int GSIZE = DSIZE;

//----------------------------------------------------------
// EXPLAIN : The input file name and directory
//----------------------------------------------------------
std::string global_in_card = "INPUT";
std::string stru_file = "STRU";
std::string global_kpoint_card = "KPT";
std::string global_wannier_card;

std::string global_pseudo_dir = "";
std::string global_orbital_dir = ""; // liuyu add 2021-08-14

// std::string global_pseudo_type = "auto";
std::string global_epm_pseudo_card;
std::string global_out_dir;
std::string global_readin_dir; // zhengdy modified
std::string global_stru_dir;
std::string global_matrix_dir;

std::ofstream ofs_running;
std::ofstream ofs_warning;
std::ofstream ofs_info; // output math lib info
std::ofstream ofs_device; // output device info

//----------------------------------------------------------
// EXPLAIN : test level for each class
//----------------------------------------------------------
int test_input = 0;
int test_winput = 0;
int test_kpoint = 0;
int test_atom = 0;
int test_unitcell = 0;
int test_symmetry = 0;

int test_fft = 0;
int test_pw = 0;
int test_elec = 0;

int test_wf = 0;
int test_charge = 0;
int test_potential = 0;
int test_energy = 0;
// for test purpose, skip ewald calculation
bool test_skip_ewald = false;
//----------------------------------------------------------
// module_hamilt_lcao/hamilt_lcaodft
//----------------------------------------------------------
int test_atom_arrange = 0;
int test_atom_input = 0;
int test_grid = 0; // 4 now
int test_grid_driver = 0; // 4 now
int test_overlap = 0;
int TEST_FORCE = 0; // mohan add 2011-03-18
int TEST_STRESS = 0; // zhengdy add 2018-05-16
int test_gridt = 0; // mohan add 2011-03-17
//----------------------------------------------------------
// src_pseudo
//----------------------------------------------------------
int test_pseudo_cell = 0; // 2 : output readin data
int test_pp = 0; // pp: pseudopotential
int test_kmesh = 0;
int test_relax_method = 0;
//----------------------------------------------------------
// src_tools
//----------------------------------------------------------
int test_deconstructor = 0;

// added by zhengdy-soc
bool NONCOLIN = false;
bool LSPINORB = false;
bool DOMAG = false;
bool DOMAG_Z = false;
int NPOL = 1;
double soc_lambda = 1.0;

bool FINAL_SCF = false; // LiuXh add 20180619

bool deepks_out_labels = false; // caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_scf = false; // caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_bandgap = false; // for bandgap label. QO added 2021-12-15
bool deepks_out_unittest = false;

bool deepks_setorb = false;

bool out_element_info = false; // added by zhengdy 2021-11-26

bool imp_sol = false; // implicit solvation.  sunml added 2022-04-04
double eb_k = 80.0;
double tau = 1.0798 * 1e-5;
double sigma_k = 0.6;
double nc_k = 0.00037;

bool dft_plus_u = false; //DFTU control

//Xinyang Dong added for rpa
bool rpa_setorb = false;
std::vector<std::string> rpa_orbitals;

std::string of_kinetic = "wt";
std::string of_method = "tn";
std::string of_conv = "energy";
double of_tole = 2e-6;
double of_tolp = 1e-5;
double of_tf_weight = 1.;
double of_vw_weight = 1.;
double of_wt_alpha = 5./6.;  
double of_wt_beta = 5./6.;
double of_wt_rho0 = 0.;
bool of_hold_rho0 = false;
double of_lkt_a = 1.3;
bool of_full_pw = true;
int of_full_pw_dim = 0;
bool of_read_kernel = false;
std::string of_kernel_file = "WTkernel.txt";

// mixing parameters
std::string MIXING_MODE = "broyden";
double MIXING_BETA = 0.7;
int MIXING_NDIM = 8;
int MIXING_RESTART = 0;
double MIXING_GG0 = 1.00;
double MIXING_BETA_MAG = 1.6;
double MIXING_GG0_MAG = 1.00;
double MIXING_GG0_MIN = 0.1;
double MIXING_ANGLE = 0.0;
bool MIXING_TAU = 0;

//==========================================================
// device flags added by denghui
//==========================================================
std::string device_flag = "unknown";
//==========================================================
// precision flags added by denghui
//==========================================================
std::string precision_flag = "unknown";

std::string chg_extrap = "";
int out_pot = 0;

std::string init_chg = "";

std::string init_wfc = "atomic";
bool psi_initializer = false;

int out_chg = 0;
double nelec = 0;
bool out_bandgap = false; // QO added for bandgap printing
int out_interval = 1;    // convert from out_hsR_interval liuyu 2023-04-18

bool out_mat_xc = false; // output Vxc in KS-wfc representation for GW calculation

//==========================================================
// Deltaspin related
//==========================================================
bool sc_mag_switch = 0;
bool decay_grad_switch = 0;
double sc_thr = 1.0e-6;
int nsc = 100;
int nsc_min = 2;
int sc_scf_nmin = 2;
double alpha_trial = 0.01; // eV/uB^2
double sccut = 3;          // eV/uB
std::string sc_file = "none";

//==========================================================
// Quasiatomic orbital related
//==========================================================
bool qo_switch = false;
std::string qo_basis = "hydrogen";
std::vector<std::string> qo_strategy = {};
double qo_thr = 1.0e-6;
std::vector<double> qo_screening_coeff = {};

//==========================================================
// PEXSI related
//==========================================================
int pexsi_npole = 54;
int pexsi_inertia = 1;
int pexsi_nmax = 80;
// int pexsi_symbolic = 1;
int pexsi_comm = 1;
int pexsi_storage = 1;
int pexsi_ordering = 0;
int pexsi_row_ordering = 1;
int pexsi_nproc = 1;
int pexsi_symm = 1;
int pexsi_trans = 0;
int pexsi_method = 1;
int pexsi_nproc_pole = 1;
// double pexsi_spin = 2;
double pexsi_temp = 0.0001;
double pexsi_gap = 0;
double pexsi_delta_e = 20.0;
double pexsi_mu_lower = -10;
double pexsi_mu_upper = 10;
double pexsi_mu = 0.0;
double pexsi_mu_thr = 0.05;
double pexsi_mu_expand = 0.3;
double pexsi_mu_guard = 0.2;
double pexsi_elec_thr = 0.001;
double pexsi_zero_thr = 1e-10;
} // namespace GlobalV

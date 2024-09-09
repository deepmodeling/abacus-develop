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
int NLOCAL = 0;        // total number of local basis.

int NSPIN = 1;       // LDA
bool TWO_EFERMI = false; // two fermi energy, exist only magnetization is fixed.
double nupdown = 0.0;
double FORCE_THR = 1.0e-3;
bool CAL_STRESS = false;
double PRESSURE = 0.0;
std::string RELAX_METHOD = "bfgs";
std::string OUT_LEVEL = "ie";
bool relax_new = true;
bool fixed_atoms = false;
int SCF_NMAX = 100;

bool use_paw = false;
bool use_uspp = false;
bool double_grid = false;
std::string KS_SOLVER = "cg";  // xiaohui add 2013-09-01
double SEARCH_RADIUS = -1.0;

int PW_DIAG_NMAX = 30;
int PW_DIAG_NDIM = 4;
double PW_DIAG_THR = 1.0e-2;
int NB2D = 1;

int SCF_THR_TYPE = 1;

double DQ = 0.010; // space between Q points of the reciprocal radial tab
int NQX = 10000;   // number of points describing reciprocal radial tab
int NQXQ = 10000;  // number of points describing reciprocal radial tab for Q

bool COLOUR = false;
bool GAMMA_ONLY_LOCAL = false; // mohan add 2010-10-20
bool GAMMA_ONLY_PW = false;    // mohan add 2012-06-05

int ZEEMAN_IN_H = 1;

// int ocp_n=0;
bool out_mul = false; // qifeng add 2019/9/10
//----------------------------------------------------------
// EXPLAIN : Parallel information
//----------------------------------------------------------

int NPROC = 1; ///< global number of process
int KPAR = 1;  ///< global number of pools
int KPAR_LCAO = 1; ///< global number of pools for LCAO diagonalization only
int NSTOGROUP = 1;
int MY_RANK = 0; ///< global index of process
int MY_POOL = 0; ///< global index of pool (count in pool)
int MY_STOGROUP = 0;
int NPROC_IN_POOL = 1; ///< local number of process in a pool
int NPROC_IN_STOGROUP = 1;
int RANK_IN_POOL = 0; ///< global index of pool (count in process), my_rank in each pool
int RANK_IN_STOGROUP = 0;
int DRANK = -1; ///< mohan add 2012-01-13, must be -1, so we can recognize who
                ///< didn't in DIAG_WORLD
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


// std::string global_pseudo_type = "auto";
std::string global_epm_pseudo_card;
std::string global_out_dir;
std::string global_readin_dir; // zhengdy modified
std::string global_stru_dir;
std::string global_matrix_dir;

std::ofstream ofs_running;
std::ofstream ofs_warning;
std::ofstream ofs_info;   // output math lib info
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

int test_pw = 0;

int test_wf = 0;
int test_charge = 0;
int test_potential = 0;
int test_energy = 0;
//----------------------------------------------------------
// module_hamilt_lcao/hamilt_lcaodft
//----------------------------------------------------------
int test_atom_input = 0;
int test_grid = 0;        // 4 now
int test_grid_driver = 0; // 4 now
int test_overlap = 0;
int TEST_FORCE = 0;  // mohan add 2011-03-18
int TEST_STRESS = 0; // zhengdy add 2018-05-16
int test_gridt = 0;  // mohan add 2011-03-17
//----------------------------------------------------------
// src_pseudo
//----------------------------------------------------------
int test_pseudo_cell = 0; // 2 : output readin data
int test_pp = 0;          // pp: pseudopotential
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


bool deepks_out_labels = false; // caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_scf = false; // caoyu add 2021-10-16 for DeePKS, wenfei 2022-1-16
bool deepks_bandgap = false; // for bandgap label. QO added 2021-12-15
int  deepks_v_delta = 0; // for v_delta label. xinyuan added 2023-2-15

bool deepks_equiv = false;

bool deepks_setorb = false;




// Xinyang Dong added for rpa
std::vector<std::string> rpa_orbitals;

std::string of_kinetic = "wt";
std::string of_method = "tn";
std::string of_conv = "energy";
double of_tole = 2e-6;
double of_tolp = 1e-5;
double of_tf_weight = 1.;
double of_vw_weight = 1.;
double of_wt_alpha = 5. / 6.;
double of_wt_beta = 5. / 6.;
double of_wt_rho0 = 0.;
bool of_hold_rho0 = false;
double of_lkt_a = 1.3;
bool of_full_pw = true;
int of_full_pw_dim = 0;
bool of_read_kernel = false;
std::string of_kernel_file = "WTkernel.txt";

// mixing parameters
double MIXING_GG0 = 1.00;
double MIXING_GG0_MAG = 1.00;
bool MIXING_DMR = false;

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

double nelec = 0;
bool out_bandgap = false; // QO added for bandgap printing

//==========================================================
// Deltaspin related
//==========================================================
double sc_thr = 1.0e-6;

//==========================================================
// Quasiatomic orbital related
//==========================================================
double qo_thr = 1.0e-6;
std::vector<double> qo_screening_coeff = {};

// on-site orbitals
double onsite_radius = 0.0;
} // namespace GlobalV

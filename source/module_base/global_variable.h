//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace GlobalV
{
//==========================================================
// EXPLAIN : Basic Global Variables
//==========================================================

extern int NBANDS;
extern int NBANDS_ISTATE; // 1.05 // mohan add 2011-03-22
extern int NLOCAL; // 1.1 // mohan add 2009-05-29

extern double KSPACING[3];
extern double MIN_DIST_COEF;

extern double PSEUDORCUT;
extern bool PSEUDO_MESH;

extern std::string CALCULATION; // 2 "scf";"nscf" ;"symmetry"
extern std::string ESOLVER_TYPE;
extern int EFIELD_FLAG; // 5 add electric field
extern int DIP_COR_FLAG; // 7 add dipole correction
extern bool GATE_FLAG;     // add gate field
extern bool out_app_flag;  // whether output r(R), H(R), S(R), T(R), and dH(R) matrices in an append manner during MD  liuyu 2023-03-20

extern std::string DFT_FUNCTIONAL; // 6.5 change the DFT functional from input file.
extern double XC_TEMPERATURE;
extern int NSPIN; // 7
extern bool TWO_EFERMI; // 7.5 two fermi energy, exist if nupdown isn't zero.
extern double nupdown;
extern int CURRENT_SPIN; // 8
extern int CURRENT_K; // 8

extern int CAL_FORCE; // 8.1
extern double FORCE_THR; // 8.2
extern bool CAL_STRESS; // 8.25 calcualte the stress
extern double PRESS1;
extern double PRESS2;
extern double PRESS3;
extern double PRESSURE;
extern std::string RELAX_METHOD;
extern std::string OUT_LEVEL;
extern int OUT_FREQ_ELEC;
extern int OUT_FREQ_ION;

extern double relax_scale_force;
extern bool relax_new;

extern bool use_paw;
extern bool use_uspp;
extern bool double_grid;

extern bool fixed_atoms;

extern int RELAX_NMAX; // 8.3
extern int SCF_NMAX;   // 8.4
extern int md_prec_level;    // liuyu 2023-03-13

extern std::string BASIS_TYPE; // xiaohui add 2013-09-01
extern std::string KS_SOLVER; // xiaohui add 2013-09-01
extern double SEARCH_RADIUS; // 11.1 // mohan add 2011-03-10
extern bool SEARCH_PBC; // 11.2 // mohan add 2011-03-10
extern bool SPARSE_MATRIX; // 11.3 // mohan add 2009-03-13

// added by zhengdy-soc
extern bool NONCOLIN; // 0 : collinear ; 1 : non-collinear
extern bool LSPINORB; // 0 : no soc ; 1 : has soc
extern bool DOMAG; // 1 : calculate the magnetism with x, y, z component
extern bool DOMAG_Z; // 1 : constrain the magnetism to z axis
extern int NPOL; // 1 : no soc; 2 : has soc
extern double soc_lambda; // soc modulator factor, from 0 to 1

extern int DIAGO_PROC; // 12.1 number of processors used to diag.
extern int PW_DIAG_NMAX; // 13
extern int DIAGO_CG_PREC; // 13.1
extern int PW_DIAG_NDIM; // 14
extern double PW_DIAG_THR; // 15 pw_diag_thr
extern int NB2D; // 16.5 dividsion of 2D_matrix.

extern double SCF_THR; // 17
extern int SCF_THR_TYPE; // type of the criterion of scf_thr, 1: reci drho for pw, 2: real drho for lcao

extern double DQ; // 19 mohan add 2009-09-10
extern int NQX; // 20 mohan add 2009-09-10
extern int NQXQ; // liuyu add 2023-10-03

extern int NURSE; // 21 mohan add 2010-09-10
extern bool COLOUR; // mohan add 2011-04-26
extern bool GAMMA_ONLY_LOCAL; // 22 : mohan add 2010-10-20
extern bool GAMMA_ONLY_PW; // mohan add 2012-06-05

extern int T_IN_H; // 23, calculate T in H or not.
extern int VL_IN_H; // 24, calculate Vl in H or not.
extern int VNL_IN_H; // 25, calculate Vnl in H or not.
extern int VH_IN_H; // 26, calculate Vh in H or not.
extern int VION_IN_H; // 28, calculate Vion_loc in H or not.
extern double STRESS_THR; // LiuXh add 20180515

extern int ocp;
// extern int ocp_n;
extern std::string ocp_set;
extern std::vector<double> ocp_kb;
// extern double ocp_kb[10000];
extern bool out_mul; // qifeng add 2019/9/10
//========================================================================
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : NPROC( global number of process )
// NAME : KPAR( global number of pools )
// NAME : MY_RANK( global index of process )
// NAME : MY_POOL( global index of pool (count in pool))
// NAME : NPROC_IN_POOL( local number of process in a pool.)
// NAME : RANK_IN_POOL( global index of pool (count in process),
//  	  MY_RANK in each pool)
// NAME : DRANK( index of diag world)
// NAME : DSIZE( number of processors in diag world, only 1 DWORLD exist)
// NAME : DCOLOR( color of each group)
// NAME : GRANK( index of grid world)
// NAME : GSIZE( number of processors in each grid world)
//========================================================================
extern int NPROC;
extern int KPAR;
extern int NSTOGROUP;
extern int MY_RANK;
extern int MY_POOL;
extern int MY_STOGROUP;
extern int NPROC_IN_POOL;
extern int NPROC_IN_STOGROUP;
extern int RANK_IN_POOL;
extern int RANK_IN_STOGROUP;
extern int DRANK;
extern int DSIZE;
extern int DCOLOR;
extern int GRANK;
extern int GSIZE;

//==========================================================
// EXPLAIN : readin file dir, output file std::ofstream
// GLOBAL VARIABLES :
// NAME : global_in_card
// NAME : stru_file
// NAME : global_kpoint_card
// NAME : global_wannier_card
// NAME : global_pseudo_dir
// NAME : global_pseudo_type // mohan add 2013-05-20 (xiaohui add 2013-06-23)
// NAME : global_out_dir
// NAME : ofs_running( contain information during runnnig)
// NAME : ofs_warning( contain warning information, including error)
//==========================================================
extern std::string global_in_card;
extern std::string stru_file;
extern std::string global_kpoint_card;
extern std::string global_wannier_card;

extern std::string global_pseudo_dir;
// extern std::string global_pseudo_type; // mohan add 2013-05-20 (xiaohui add 2013-06-23)
extern std::string global_out_dir;
extern std::string global_orbital_dir; // liuyu add 2021-08-14
extern std::string global_readin_dir; // zhengdy modified
extern std::string global_stru_dir;   // liuyu add 2022-05-24 for MD STRU
extern std::string global_matrix_dir; // liuyu add 2022-09-19 for HS matrix outpu, jiyy modified 2023-01-23 for R matrix output

extern std::ofstream ofs_running;
extern std::ofstream ofs_warning;
extern std::ofstream ofs_info;
extern std::ofstream ofs_device;

//==========================================================
// EXPLAIN : test level for each class
//==========================================================
extern int test_input;
extern int test_winput;
extern int test_kpoint;
extern int test_atom;
extern int test_unitcell;
extern int test_symmetry;

extern int test_fft;
extern int test_pw;
extern int test_elec;

extern int test_wf;
extern int test_charge;
extern int test_potential;
extern int test_energy;
extern bool test_skip_ewald;
//==========================================================
// src_onscaling
//==========================================================
extern int test_atom_arrange;
extern int test_atom_input;
extern int test_grid;
extern int test_grid_driver;
extern int test_overlap;
extern int TEST_FORCE; // mohan add 2011-03-18
extern int TEST_STRESS; // zhengdy add 2018-05-16
extern int test_gridt; // mohan add 2011-03-17
//==========================================================
// src_pseudo
//==========================================================
extern int test_pseudo_cell;
extern int test_pp;
extern int test_kmesh;
extern int test_relax_method;
//==========================================================
// src_tools
//==========================================================
extern int test_deconstructor;

extern bool FINAL_SCF; // LiuXh add 20180619

extern bool
    deepks_out_labels; // (need libnpy) prints energy and force labels and descriptors for training, wenfei 2022-1-12
extern bool
    deepks_scf; //(need libnpy and libtorch) if set 1, a trained model would be needed to cal V_delta and F_delta
extern bool deepks_bandgap; // for bandgap label. QO added 2021-12-15

extern bool deepks_setorb;

extern bool deepks_out_unittest; // if set 1, prints intermediate quantities that shall be used for making unit test

extern std::string deepks_model; // needed when deepks_scf=1

// the following 3 are used when generating jle.orb
extern int bessel_lmax; // lmax used in descriptor, mohan added 2021-01-03
extern double bessel_rcut;
extern double bessel_tol;

// whether or not output information for each element
extern bool out_element_info;

// implicit solvation
extern bool imp_sol; // sunml added 2022-04-04
extern double eb_k;
extern double tau;
extern double sigma_k;
extern double nc_k;

// DFTU control
extern bool dft_plus_u;
// rpa related
extern bool rpa_setorb;
extern std::vector<std::string> rpa_orbitals;

// ofdft sunliang add on 2022-05-11
extern std::string of_kinetic; // Kinetic energy functional, such as TF, VW, WT
extern std::string of_method;  // optimization method, include cg1, cg2, tn (default), bfgs
extern std::string of_conv;    // select the convergence criterion, potential, energy (default), or both
extern double of_tole;    // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
extern double of_tolp;    // tolerance of potential for determining the convergence, default=1e-5 in a.u.
extern double of_tf_weight;  // weight of TF KEDF
extern double of_vw_weight;  // weight of vW KEDF
extern double of_wt_alpha;   // parameter alpha of WT KEDF
extern double of_wt_beta;    // parameter beta of WT KEDF
extern double of_wt_rho0;    // set the average density of system, in Bohr^-3
extern bool of_hold_rho0;   // If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be set to 1 automaticly if of_wt_rho0 is not zero.
extern double of_lkt_a;    // parameter a of LKT KEDF
extern bool of_full_pw;     // If set to 1, ecut will be ignored while collecting planewaves, so that all planewaves will be used.
extern int of_full_pw_dim;  // If of_full_pw = 1, the dimention of FFT will be testricted to be (0) either odd or even; (1) odd only; (2) even only.
extern bool of_read_kernel; // If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from formula. Only usable for WT KEDF.
extern std::string of_kernel_file; // The name of WT kernel file.

// mixing parameters
extern std::string MIXING_MODE;
extern double MIXING_BETA;
extern int MIXING_NDIM;
extern double MIXING_GG0;
extern bool MIXING_TAU;
extern double MIXING_BETA_MAG;
extern double MIXING_GG0_MAG;
extern double MIXING_GG0_MIN;

//==========================================================
// device flags added by denghui
//==========================================================
extern std::string device_flag;
//==========================================================
// precision flags added by denghui
//==========================================================
extern std::string precision_flag;

extern std::string chg_extrap;
extern int out_pot;

extern std::string init_chg; //  output charge if out_chg > 0, and output every "out_chg" elec step.
/// @brief method to initialize wavefunction
/// @author kirk0830, 20230920
extern std::string init_wfc; 
/// @brief whether use the new psi initializer to initialize psi
/// @author ykhuang, 20230920
extern bool psi_initializer;
extern int out_chg;

extern double nelec;
extern bool out_bandgap;
extern int out_interval;

// Deltaspin related
extern bool sc_mag_switch; // 0: no deltaspin; 1: constrain atomic magnetic moments;
extern bool decay_grad_switch; // 0: decay grad will be set to zero; 1: with decay grad set for some elements
extern double sc_thr;
extern int nsc;
extern int nsc_min;
extern int sc_scf_nmin;
extern double alpha_trial;
extern double sccut;
extern std::string sc_file;

// Quasiatomic orbitals related
extern bool out_kpt_reduction;
} // namespace GlobalV
#endif

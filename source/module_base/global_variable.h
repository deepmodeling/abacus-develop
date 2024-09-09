//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace GlobalV
{
//==========================================================
// EXPLAIN : Basic Global Variables
//==========================================================

extern int NBANDS;
extern int NLOCAL;        // 1.1 // mohan add 2009-05-29


extern int NSPIN;       // 7
extern bool TWO_EFERMI; // 7.5 two fermi energy, exist if nupdown isn't zero.
extern double nupdown;

extern double FORCE_THR; // 8.2
extern bool CAL_STRESS;  // 8.25 calcualte the stress

extern double PRESSURE;
extern std::string RELAX_METHOD;
extern std::string OUT_LEVEL;

extern bool relax_new;

extern bool use_paw;
extern bool use_uspp;
extern bool double_grid;

extern bool fixed_atoms;

extern int SCF_NMAX;      // 8.4

extern std::string KS_SOLVER;  // xiaohui add 2013-09-01
extern double SEARCH_RADIUS;   // 11.1 // mohan add 2011-03-10

// added by zhengdy-soc
extern bool NONCOLIN;     // 0 : collinear ; 1 : non-collinear
extern bool LSPINORB;     // 0 : no soc ; 1 : has soc
extern bool DOMAG;        // 1 : calculate the magnetism with x, y, z component
extern bool DOMAG_Z;      // 1 : constrain the magnetism to z axis
extern int NPOL;          // 1 : no soc; 2 : has soc

extern int PW_DIAG_NMAX;   // 13
extern int PW_DIAG_NDIM;   // 14
extern double PW_DIAG_THR; // 15 pw_diag_thr
extern int NB2D;           // 16.5 dividsion of 2D_matrix.

extern int SCF_THR_TYPE; // type of the criterion of scf_thr, 1: reci drho for
                         // pw, 2: real drho for lcao

extern double DQ; // 19 mohan add 2009-09-10
extern int NQX;   // 20 mohan add 2009-09-10
extern int NQXQ;  // liuyu add 2023-10-03

extern bool COLOUR;           // mohan add 2011-04-26
extern bool GAMMA_ONLY_LOCAL; // 22 : mohan add 2010-10-20
extern bool GAMMA_ONLY_PW;    // mohan add 2012-06-05


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

//========================================================================
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : KPAR_LCAO ( global number of pools for LCAO diagonalization only)
//========================================================================
extern int KPAR_LCAO;

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

// extern std::string global_pseudo_type; // mohan add 2013-05-20 (xiaohui add
// 2013-06-23)
extern std::string global_out_dir;
extern std::string global_readin_dir;  // zhengdy modified
extern std::string global_stru_dir;    // liuyu add 2022-05-24 for MD STRU
extern std::string global_matrix_dir;  // liuyu add 2022-09-19 for HS matrix outpu, jiyy
                                       // modified 2023-01-23 for R matrix output

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

extern int test_pw;

extern int test_wf;
extern int test_charge;
extern int test_potential;
extern int test_energy;
//==========================================================
// src_onscaling
//==========================================================
extern int test_atom_input;
extern int test_grid;
extern int test_grid_driver;
extern int test_overlap;
extern int TEST_FORCE;  // mohan add 2011-03-18
extern int TEST_STRESS; // zhengdy add 2018-05-16
extern int test_gridt;  // mohan add 2011-03-17
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


extern bool deepks_out_labels; // (need libnpy) prints energy and force labels
                               // and descriptors for training, wenfei 2022-1-12
extern bool deepks_scf;        //(need libnpy and libtorch) if set 1, a trained model
                               // would be needed to cal V_delta and F_delta
extern bool deepks_bandgap;    // for bandgap label. QO added 2021-12-15

extern int deepks_v_delta; // for v_delta label. xinyuan added 2023-2-15

extern bool deepks_equiv; //whether to use equviariant version of DeePKS

extern bool deepks_setorb;



// implicit solvation

// DFTU control
// rpa related
extern std::vector<std::string> rpa_orbitals;

// mixing parameters
extern double MIXING_GG0;
extern double MIXING_GG0_MAG;

//==========================================================
// device flags added by denghui
//==========================================================
extern std::string device_flag;
//==========================================================
// precision flags added by denghui
//==========================================================

extern std::string chg_extrap;
extern int out_pot;

extern std::string init_chg; //  output charge if out_chg > 0, and output every
                             //  "out_chg" elec step.
/// @brief method to initialize wavefunction
/// @author kirk0830, 20230920
extern std::string init_wfc;
/// @brief whether use the new psi initializer to initialize psi
/// @author ykhuang, 20230920
extern bool psi_initializer;

extern double nelec;
extern bool out_bandgap;

// Deltaspin related
extern double sc_thr;

// Quasiatomic orbital related
extern double qo_thr;
extern std::vector<double> qo_screening_coeff;

// radius of on-site orbitals
extern double onsite_radius;
} // namespace GlobalV
#endif

//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-07
//==========================================================
#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
//==========================================================
// EXPLAIN : Basic Global Variables
//==========================================================

extern int 		NBANDS;	
extern int		NBANDS_ISTATE;		// 1.05 // mohan add 2011-03-22
extern int 		NLOCAL;				// 1.1 // mohan add 2009-05-29

extern string	CALCULATION;		// 2 "scf";"nscf" ;"symmetry"
extern int		EFIELD;				// 5 add electric field
extern int 		DIPOLE;				// 7 add dipole correction


extern string   DFT_FUNCTIONAL;		// 6.5 change the DFT functional from input file.
extern int 		NSPIN;				// 7
extern bool		TWO_EFERMI; 		// 7.5 mohan add 2011-04-03, two fermi energy, exist if magnetization is fixed.
extern int 		CURRENT_SPIN;		// 8

extern int 		FORCE;				// 8.1
extern double	FORCE_THR;			// 8.2
extern bool 	STRESS;				// 8.25 calcualte the stress
extern double  PRESS1;
extern double  PRESS2;
extern double  PRESS3;
extern double PRESSURE;
extern string	MOVE_IONS;
extern string	OUT_LEVEL;

extern int		NSTEP;				// 8.3
extern int 		NITER;				// 8.4

extern string	BASIS_TYPE; //xiaohui add 2013-09-01
extern string	KS_SOLVER; //xiaohui add 2013-09-01
extern double	SEARCH_RADIUS;		// 11.1 // mohan add 2011-03-10
extern bool		SEARCH_PBC;			// 11.2 // mohan add 2011-03-10
extern bool		SPARSE_MATRIX;		// 11.3 // mohan add 2009-03-13

//added by zhengdy-soc
extern bool     NONCOLIN; //0 : collinear ; 1 : non-collinear
extern bool     LSPINORB; //0 : no soc ; 1 : has soc
extern bool     DOMAG;    //1 : calculate the magnetism with x, y, z component
extern bool     DOMAG_Z;  //1 : constrain the magnetism to z axis
extern int      NPOL;     //1 : no soc; 2 : has soc
extern int      PRENSPIN; //NSPIN used before, for restart with soc


extern int		DIAGO_PROC;			// 12.1 number of processors used to diag.
extern int 		DIAGO_CG_MAXITER;	// 13
extern int		DIAGO_CG_PREC;		// 13.1
extern int 		DIAGO_DAVID_NDIM;	// 14
extern double 	ETHR;				// 15 ethr
extern double	FS_REF_ENERGY;		// 16
extern int		NB2D;				// 16.5 dividsion of 2D_matrix.

extern double 	DRHO2;				// 17

extern string	RESTART_MODE;		// 18

extern double DQ;	// 19 mohan add 2009-09-10
extern int NQX;	// 20 mohan add 2009-09-10

extern int NURSE; // 21 mohan add 2010-09-10
extern bool COLOUR; // mohan add 2011-04-26
extern bool GAMMA_ONLY_LOCAL; // 22 : mohan add 2010-10-20
extern bool GAMMA_ONLY_PW; // mohan add 2012-06-05

extern int T_IN_H; // 23, calculate T in H or not.
extern int VL_IN_H; // 24, calculate Vl in H or not.
extern int VNL_IN_H; // 25, calculate Vnl in H or not.
extern double   STRESS_THR; //LiuXh add 20180515

extern int tddft;	//fuxiang add 2016-10-27
extern double td_dr2;
extern double td_dt;
extern double td_force_dt;
extern int val_elec_01;
extern int val_elec_02;
extern int val_elec_03;
extern int vext;
extern int vext_dire;

extern int ocp;
//extern int ocp_n;
extern string ocp_set;
extern vector<double> ocp_kb;
//extern double ocp_kb[10000];
extern int  mulliken;//qifeng add 2019/9/10
//========================================================================
// EXPLAIN : Parallel information
// GLOBAL VARIABLES :
// NAME : NPROC( global number of process )
// NAME : NPOOL( global number of pools )
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
extern int NPOOL;
extern int MY_RANK;
extern int MY_POOL;
extern int NPROC_IN_POOL;
extern int RANK_IN_POOL;
extern int DRANK;
extern int DSIZE;
extern int DCOLOR;
extern int GRANK;
extern int GSIZE;

//==========================================================
// EXPLAIN : readin file dir, output file ofstream
// GLOBAL VARIABLES :
// NAME : global_in_card
// NAME : global_atom_card
// NAME : global_kpoint_card
// NAME : global_wannier_card
// NAME : global_pseudo_dir
// NAME : global_pseudo_type // mohan add 2013-05-20 (xiaohui add 2013-06-23)
// NAME : global_out_dir
// NAME : ofs_running( contain information during runnnig)
// NAME : ofs_warning( contain warning information, including error)
//==========================================================
extern string	global_in_card;
extern string	global_atom_card;
extern string 	global_kpoint_card;
extern string	global_wannier_card;

extern string	global_pseudo_dir;
extern string   global_pseudo_type; // mohan add 2013-05-20 (xiaohui add 2013-06-23)
extern string 	global_out_dir;

extern ofstream ofs_running;
extern ofstream ofs_warning;

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
//==========================================================
// src_onscaling
//==========================================================
extern int test_atom_arrange;
extern int test_atom_input;
extern int test_grid;
extern int test_grid_driver;
extern int test_overlap;
extern int TEST_FORCE;// mohan add 2011-03-18
extern int TEST_STRESS;//zhengdy add 2018-05-16
extern int test_gridt; // mohan add 2011-03-17
//==========================================================
// src_pseudo
//==========================================================
extern int test_pseudo_cell;
extern int test_pp;
extern int test_kmesh;
extern int test_ion_dynamics;
//==========================================================
// src_tools
//==========================================================
extern int test_deconstructor;

extern bool FINAL_SCF; //LiuXh add 20180619
extern int NEW_DM;  // Shen Yu add 2019/5/9

#endif

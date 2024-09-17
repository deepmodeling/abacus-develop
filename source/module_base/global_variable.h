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

extern bool CAL_STRESS;  // 8.25 calcualte the stress

extern std::string RELAX_METHOD;


extern bool use_uspp;

extern std::string KS_SOLVER;  // xiaohui add 2013-09-01
extern double SEARCH_RADIUS;   // 11.1 // mohan add 2011-03-10


extern int PW_DIAG_NDIM;   // 14
extern double PW_DIAG_THR; // 15 pw_diag_thr
extern int NB2D;           // 16.5 dividsion of 2D_matrix.

                         // pw, 2: real drho for lcao

extern double DQ; // 19 mohan add 2009-09-10
extern int NQX;   // 20 mohan add 2009-09-10
extern int NQXQ;  // liuyu add 2023-10-03


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
// NAME : KPAR_LCAO ( global number of pools for LCAO diagonalization only)
//========================================================================
extern int NPROC;
extern int KPAR;
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
extern std::string stru_file;
// extern std::string global_pseudo_type; // mohan add 2013-05-20 (xiaohui add
// 2013-06-23)
extern std::ofstream ofs_running;
extern std::ofstream ofs_warning;
extern std::ofstream ofs_info;
extern std::ofstream ofs_device;

// mixing parameters

//==========================================================
// device flags added by denghui
//==========================================================
extern std::string device_flag;
//==========================================================
// precision flags added by denghui
//==========================================================

                             //  "out_chg" elec step.
/// @brief method to initialize wavefunction
/// @author kirk0830, 20230920
/// @brief whether use the new psi initializer to initialize psi
/// @author ykhuang, 20230920

extern double nelec;

// Deltaspin related

// Quasiatomic orbital related

// radius of on-site orbitals
} // namespace GlobalV
#endif

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
// NAME : ofs_running( contain information during runnnig)
// NAME : ofs_warning( contain warning information, including error)
//==========================================================
extern std::ofstream ofs_running;
extern std::ofstream ofs_warning;
extern std::ofstream ofs_info;
extern std::ofstream ofs_device;
} // namespace GlobalV
#endif

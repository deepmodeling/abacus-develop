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
double nupdown = 0.0;

bool use_uspp = false;
std::string KS_SOLVER = "cg";  // xiaohui add 2013-09-01
double SEARCH_RADIUS = -1.0;

int NB2D = 1;


//----------------------------------------------------------
// EXPLAIN : Parallel information
//----------------------------------------------------------

int NPROC = 1; ///< global number of process
int KPAR = 1;  ///< global number of pools
int KPAR_LCAO = 1; ///< global number of pools for LCAO diagonalization only
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


std::ofstream ofs_running;
std::ofstream ofs_warning;
std::ofstream ofs_info;   // output math lib info
std::ofstream ofs_device; // output device info

//----------------------------------------------------------
// src_tools
//----------------------------------------------------------

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

std::vector<std::string> rpa_orbitals;

//==========================================================
// device flags added by denghui
//==========================================================
std::string device_flag = "unknown";

double nelec = 0;


// on-site orbitals
} // namespace GlobalV

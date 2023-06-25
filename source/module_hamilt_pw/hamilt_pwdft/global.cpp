#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------
namespace GlobalC
{
#ifdef __EXX
Exx_Info exx_info;
Exx_Lip exx_lip(exx_info.info_lip);
#endif
pseudopot_cell_vnl ppcell;
UnitCell ucell;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07
Restart restart; // Peize Lin add 2020.04.04
}

//Magnetism mag;															

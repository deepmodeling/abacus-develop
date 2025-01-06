#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------
namespace GlobalC
{
#ifdef __EXX
    Exx_Info exx_info;
#endif
UnitCell ucell;
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07
Restart restart; // Peize Lin add 2020.04.04
}

//Magnetism mag;															

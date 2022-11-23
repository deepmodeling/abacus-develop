#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------
namespace GlobalC
{
K_Vectors kv; // mem check in in here.
Structure_Factor sf;
ModulePW::PW_Basis* rhopw;
ModulePW::PW_Basis_Big* bigpw;
ModulePW::PW_Basis_K* wfcpw;
energy en;
wavefunc wf;
#ifdef __MPI
#ifdef __EXX
Exx_Global exx_global;
Exx_Lip exx_lip(exx_global.info);
#endif
#endif
pseudopot_cell_vnl ppcell;
UnitCell ucell;
Charge_Mixing CHR_MIX;
ModuleSymmetry::Symmetry symm;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07
Restart restart; // Peize Lin add 2020.04.04
}

//Magnetism mag;															

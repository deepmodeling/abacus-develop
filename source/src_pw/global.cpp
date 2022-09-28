#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------
namespace GlobalC
{
K_Vectors kv; // mem check in in here.
Use_FFT UFFT; // mohan add 2010-07-22
Structure_Factor sf;
ModulePW::PW_Basis* rhopw;
ModulePW::PW_Basis_Big* bigpw;
ModulePW::PW_Basis_K* wfcpw;
energy en;
wavefunc wf;
Hamilt hm;
#ifdef __MPI
#ifdef __EXX
Exx_Global exx_global;
Exx_Lip exx_lip(exx_global.info);
#endif
#endif
pseudopot_cell_vnl ppcell;
UnitCell_pseudo ucell;
Charge_Broyden CHR;
Potential pot;
ModuleSymmetry::Symmetry symm;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07
Vdwd2_Parameters vdwd2_para;	// Peize Lin add 2021-03-09
Vdwd3_Parameters vdwd3_para;	// jiyy add 2021-05-02
Restart restart; // Peize Lin add 2020.04.04
}

//Magnetism mag;															

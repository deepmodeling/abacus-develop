#include "global.h"
//----------------------------------------------------------
// init "GLOBAL CLASS" object
//----------------------------------------------------------

kvect kv; // mem check in in here.
Use_FFT UFFT; // mohan add 2010-07-22
FFT fftwan;
output out;

PW_Basis pw;
Stochastic_WF STO_WF;
energy en;
wavefunc wf;
Hamilt hm;
Exx_Global exx_global;
Exx_Lip exx_lip(exx_global.info);

UnitCell_pseudo ucell;
pseudopot_cell_vnl ppcell;
xcfunc xcf;
Charge_Broyden CHR;
Magnetism mag;

Potential pot;
Symmetry symm;
Parallel_Grid Pgrid; //mohan add 2010-06-06 
Parallel_Kpoints Pkpoints; // mohan add 2010-06-07

Vdwd2_Parameters vdwd2_para;	// Peize Lin add 2021-03-09
Vdwd3_Parameters vdwd3_para;	// jiyy add 2021-05-02			

Soc soc; //added by zhengdy-soc														

Restart restart; // Peize Lin add 2020.04.04

#ifndef GLOBAL_H
#define GLOBAL_H

#include "../run_pw.h"
#include "../src_global/global_variable.h"
#include "../src_global/global_function.h"
#include "pw_basis.h"
#include "energy.h"
#include "VNL_in_pw.h"
#include "charge_broyden.h"
#include "potential.h"
#include "xc_type.h"
#include "hamilt.h"
#include "../src_ions/ions.h"
#include "wavefunc.h"
#include "use_fft.h"
#include "klist.h"
#include "../src_io/output.h"
#include "magnetism.h"
#include "vdwd2.h"
#include "vdwd2_parameters.h"
#include "vdwd3.h"
#include "vdwd3_parameters.h"	
#include "../src_global/restart.h" 
#include "exx_global.h"
#include "exx_lip.h"
#include "../src_parallel/ft.h"

//==========================================================
// EXPLAIN : define "GLOBAL CLASS"
//==========================================================
extern Use_FFT UFFT;
extern FFT fftwan;
extern K_Vectors kv;
extern output out;

extern PW_Basis pw;
extern Stochastic_WF STO_WF; //qianrui add 2021-2-5
extern energy en;
extern wavefunc wf;
extern Hamilt hm;
extern Exx_Global exx_global;
extern Exx_Lip exx_lip;

#include "symmetry.h"
#include "module_cell/unitcell_pseudo.h"
#include "../src_parallel/parallel_grid.h"
#include "../src_parallel/parallel_kpoints.h"
extern pseudopot_cell_vnl ppcell;
extern xcfunc xcf;
extern Charge_Broyden CHR;
extern Potential pot;
extern Symmetry symm;
extern Magnetism mag;
extern UnitCell_pseudo ucell;
extern Parallel_Grid Pgrid; 
extern Parallel_Kpoints Pkpoints;
extern Vdwd2_Parameters vdwd2_para;		// Peize Lin add 2021.03.09
extern Vdwd3_Parameters vdwd3_para;		// jiyy add 2021-05-02	
extern Restart restart;	// Peize Lin add 2020.04.04

#include "../src_lcao/global_fp.h"


#endif

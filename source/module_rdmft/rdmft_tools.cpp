//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
// #include "rdmft.h"
#include "module_rdmft/rdmft_tools.h"

#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_psi/psi.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>

// used by class Veff_rdmft
//#include "veff_lcao.h"
//#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/pot_local.h"
#include "module_elecstate/potentials/pot_xc.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

// for test use dgemm_
#include "module_base/matrix.h"
#include "module_base/blas_connector.h"

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"  //test





namespace rdmft
{

















}




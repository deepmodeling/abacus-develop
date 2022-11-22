#ifndef GLOBAL_FP_H
#define GLOBAL_FP_H

#include "../module_neighbor/sltk_grid_driver.h"
#include "../module_gint/grid_technique.h"
#include "src_pdiag/pdiag_double.h"
#include "local_orbital_wfc.h"
#include "local_orbital_charge.h"
#include "LCAO_matrix.h"
#include "LCAO_gen_fixedH.h"
#include "LCAO_hamilt.h" 
#include "../module_orbital/ORB_read.h"
#include "../module_orbital/ORB_gen_tables.h"
#ifdef __EXX
#include "../src_ri/exx_lcao.h"
#include "module_ri/Exx_LRI.h"
#endif

namespace GlobalC
{
extern Grid_Driver GridD;
extern Pdiag_Double ParaO;
#ifdef __EXX
extern Exx_Lcao exx_lcao; // Peize Lin add 2016-12-03
extern Exx_LRI<double> exx_lri_double; // Peize Lin add 2022-08-06
extern Exx_LRI<std::complex<double>> exx_lri_complex; // Peize Lin add 2022-08-06
#endif
}
#endif

#ifndef CTRL_OUTPUT_FP_H 
#define CTRL_OUTPUT_FP_H 

#include "source_estate/elecstate_lcao.h"

namespace ModuleIO
{

	void ctrl_output_fp(UnitCell& ucell, 
			elecstate::ElecState* pelec,	
			ModulePW::PW_Basis_Big* pw_big,
			ModulePW::PW_Basis* pw_rhod,
			Charge &chr,
			surchem &solvent,
			Parallel_Grid &para_grid,
			const int istep);

}
#endif

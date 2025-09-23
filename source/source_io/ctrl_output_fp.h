#ifndef CTRL_OUTPUT_FP_H 
#define CTRL_OUTPUT_FP_H 

#include "source_estate/elecstate_lcao.h"

namespace ModuleIO
{
	void ctrl_output_fp(UnitCell& ucell, 
		    elecstate::ElecState* pelec,	
			const int istep);
}
#endif

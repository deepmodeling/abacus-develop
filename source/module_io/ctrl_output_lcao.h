#ifndef CTRL_OUTPUT_LCAO_H 
#define CTRL_OUTPUT_LCAO_H 

#include <complex>

#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate_lcao.h" // use elecstate::ElecState

namespace ModuleIO
{
	template <typename TK, typename TR>
		void ctrl_output_lcao(const UnitCell& ucell, 
				const elecstate::ElecStateLCAO<TK>* pelec, 
				const Parallel_Orbitals& pv,
				const int istep);
}
#endif

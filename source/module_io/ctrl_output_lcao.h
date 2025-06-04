#ifndef CTRL_OUTPUT_LCAO_H 
#define CTRL_OUTPUT_LCAO_H 

#include <complex>

#include "module_cell/unitcell.h" // use UnitCell
#include "module_cell/klist.h" // use K_Vectors
#include "module_elecstate/elecstate_lcao.h" // use elecstate::ElecStateLCAO<TK> 
#include "module_psi/psi.h" // use Psi<TK>
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>

namespace ModuleIO
{
	template <typename TK, typename TR>
		void ctrl_output_lcao(const UnitCell& ucell, 
				const K_Vectors& kv,
				const elecstate::ElecStateLCAO<TK>* pelec, 
				const Parallel_Orbitals& pv,
				const psi::Psi<TK>* psi,
				hamilt::HamiltLCAO<TK, TR>* p_hamilt,
				const int istep);
}
#endif

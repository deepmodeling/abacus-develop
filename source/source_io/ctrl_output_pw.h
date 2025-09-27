#ifndef CTRL_OUTPUT_PW_H 
#define CTRL_OUTPUT_PW_H 

#include "source_estate/elecstate_lcao.h"

namespace ModuleIO
{

template <typename Device>
void ctrl_iter_pw(const int istep, 
		const int iter, 
		const double &conv_esolver,
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
        const Input_para& inp);

template <typename Device>
void ctrl_scf_pw(elecstate::ElecState* pelec);

template <typename Device>
void ctrl_runner_pw(UnitCell& ucell, 
		elecstate::ElecState* pelec,	
		ModulePW::PW_Basis_Big* pw_big,
		ModulePW::PW_Basis* pw_rhod,
		Charge &chr,
		surchem &solvent,
		Parallel_Grid &para_grid,
		const int istep);

}
#endif

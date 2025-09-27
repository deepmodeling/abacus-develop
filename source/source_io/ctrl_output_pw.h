#ifndef CTRL_OUTPUT_PW_H 
#define CTRL_OUTPUT_PW_H 

#include "source_psi/psi_init.h" // about psi
#include "source_estate/elecstate_lcao.h" // use pelec

namespace ModuleIO
{

void ctrl_iter_pw(const int istep, 
		const int iter, 
		const double &conv_esolver,
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
        const Input_para& inp);

template <typename T, typename Device>
void ctrl_scf_pw(elecstate::ElecState* pelec,
        const Charge &chr,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
		const ModulePW::PW_Basis *pw_rho,
		const ModulePW::PW_Basis *pw_rhod,
		const ModulePW::PW_Basis_Big *pw_big,
        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
        psi::Psi<T, Device>* kspw_psi,
        psi::Psi<std::complex<double>, Device>* __kspw_psi,
        const Device* ctx,
        const Parallel_Grid &para_grid,
        const Input_para& inp);

template <typename T, typename Device>
void ctrl_runner_pw(UnitCell& ucell, 
		elecstate::ElecState* pelec,	
        ModulePW::PW_Basis_K* pw_wfc,
        ModulePW::PW_Basis* pw_rho,
        ModulePW::PW_Basis* pw_rhod,
		Charge &chr,
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
		psi::Psi<T, Device>* kspw_psi,
		psi::Psi<std::complex<double>, Device>* __kspw_psi,
        Structure_Factor &sf,
        pseudopot_cell_vnl &ppcell,
		surchem &solvent,
        const Device* ctx,
        Parallel_Grid &para_grid,
        const Input_para& inp);

}
#endif

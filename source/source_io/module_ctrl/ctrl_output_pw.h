#ifndef CTRL_OUTPUT_PW_H 
#define CTRL_OUTPUT_PW_H 

#include "source_base/module_device/device.h" // use Device
#include "source_psi/psi.h"                   // define psi
#include "source_estate/elecstate_lcao.h"     // use pelec
#include "source_psi/setup_psi_pw.h" // use Setup_Psi class

namespace ModuleIO
{

// print out information in 'iter_finish' in ESolver_KS_PW
void ctrl_iter_pw(const int istep, 
		const int iter, 
		const double &conv_esolver,
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
        const Input_para& inp);

// print out information in 'after_scf' in ESolver_KS_PW
template <typename T, typename Device>
void ctrl_scf_pw(const int istep,
        UnitCell& ucell,
        elecstate::ElecState* pelec,
        const Charge &chr,
		const K_Vectors &kv,
		const ModulePW::PW_Basis_K *pw_wfc,
		const ModulePW::PW_Basis *pw_rho,
		const ModulePW::PW_Basis *pw_rhod,
		const ModulePW::PW_Basis_Big *pw_big,
        Setup_Psi_pw<T, Device> &stp,
        const Device* ctx,
        const base_device::AbacusDevice_t &device, // mohan add 2025-10-15
        const Parallel_Grid &para_grid,
        const Input_para& inp);

// print out information in 'after_all_runners' in ESolver_KS_PW
template <typename T, typename Device>
void ctrl_runner_pw(UnitCell& ucell, 
		elecstate::ElecState* pelec,	
        ModulePW::PW_Basis_K* pw_wfc,
        ModulePW::PW_Basis* pw_rho,
        ModulePW::PW_Basis* pw_rhod,
		Charge &chr,
        K_Vectors &kv,
        Setup_Psi_pw<T, Device> &stp,
        Structure_Factor &sf,
        pseudopot_cell_vnl &ppcell,
		surchem &solvent,
        const Device* ctx,
        Parallel_Grid &para_grid,
        const Input_para& inp);

}
#endif

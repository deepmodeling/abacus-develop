#ifndef SETUP_ESTATE_PW_H
#define SETUP_ESTATE_PW_H

#include "source_base/module_device/device.h" // use Device
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_estate/elecstate.h"
#include "source_pw/module_pwdft/VL_in_pw.h"
#include "source_pw/module_pwdft/VSep_in_pw.h"

namespace elecstate
{

template <typename T, typename Device>
void setup_estate_pw(UnitCell& ucell, // unitcell 
		K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState* &pelec, // pointer of electrons
		Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* &vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K* pw_wfc,  // pw for wfc
		ModulePW::PW_Basis* pw_rho, // pw for rho
		ModulePW::PW_Basis* pw_rhod, // pw for rhod
        ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp); // input parameters

template <typename T, typename Device>
void teardown_estate_pw(elecstate::ElecState* &pelec, VSep* &vsep_cell); 

}


#endif

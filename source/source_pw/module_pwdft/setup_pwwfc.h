#ifndef SETUP_PWWFC_H
#define SETUP_PWWFC_H

#include "source_io/module_parameter/parameter.h" // input parameters
#include "source_cell/unitcell.h" // cell information
#include "source_cell/klist.h" // k-points
#include "source_basis/module_pw/pw_basis.h" // pw_rho
#include "source_basis/module_pw/pw_basis_k.h" // pw_wfc 

namespace pw
{

void teardown_pwwfc(ModulePW::PW_Basis_K* &pw_wfc);

void setup_pwwfc(const Input_para& inp,
		const UnitCell& ucell, 
		const ModulePW::PW_Basis& pw_rho,
		K_Vectors& kv,
		ModulePW::PW_Basis_K* &pw_wfc);

}

#endif

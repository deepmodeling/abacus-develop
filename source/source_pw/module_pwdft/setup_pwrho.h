#ifndef SETUP_PWRHO_H
#define SETUP_PWRHO_H

#include "source_cell/unitcell.h" // use UnitCell
#include "source_pw/module_pwdft/structure_factor.h" // use Structure_Factor
#include "source_basis/module_pw/pw_basis.h" // use PW_Basis
#include "source_io/module_parameter/input_parameter.h" // use Input_para

namespace pw
{

void setup_pwrho(
		UnitCell& ucell, // unitcell 
        const bool double_grid, // for USPP
        bool &pw_rho_flag, // flag for allocation of pw_rho
		ModulePW::PW_Basis* &pw_rho, // pw for rhod
		ModulePW::PW_Basis* &pw_rhod, // pw for rhod
		ModulePW::PW_Basis_Big* &pw_big, // pw for rhod
		const std::string &classname,
		const Input_para& inp); // input parameters *


void teardown_pwrho(bool &pw_rho_flag,
		const bool double_grid,
		ModulePW::PW_Basis* &pw_rho, // pw for rhod
		ModulePW::PW_Basis* &pw_rhod); // pw for rhod

}



#endif

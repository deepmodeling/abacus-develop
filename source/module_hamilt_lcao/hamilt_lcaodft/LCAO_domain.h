#ifndef LCAO_DOMAIN_H
#define LCAO_DOMAIN_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"


namespace LCAO_domain
{

//! prepare grid integration 
void grid_prepare(
		const Grid_Technique& gt, 
        Gint_Gamma &gint_gamma,
        Gint_k &gint_k,
		const ModulePW::PW_Basis& rhopw, 
		const ModulePW::PW_Basis_Big& bigpw);


}

#endif

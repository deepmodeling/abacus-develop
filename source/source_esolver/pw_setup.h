#ifndef PW_SETUP_H
#define PW_SETUP_H

//! Input parameters
#include "source_io/module_parameter/parameter.h"
//! Unit cell information
#include "source_cell/unitcell.h"
//! Plane wave basis
#include "source_basis/module_pw/pw_basis.h"
//! K points in Brillouin zone
#include "source_cell/klist.h"
//! Plane wave basis set for k points
#include "source_basis/module_pw/pw_basis_k.h"

namespace ModuleESolver
{

void pw_setup(const Input_para& inp,
		const UnitCell& ucell, 
		const ModulePW::PW_Basis& pw_rho,
		K_Vectors& kv,
		ModulePW::PW_Basis_K& pw_wfc);

}

#endif

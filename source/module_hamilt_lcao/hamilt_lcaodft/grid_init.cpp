#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

namespace LCAO_domain
{

//--------------------------------------------
// prepare grid network for Gint(grid integral)
//--------------------------------------------
void grid_prepare(
		const Grid_Technique& gt, 
        Gint_Gamma &gint_gamma,
        Gint_k &gint_k,
		const UnitCell& ucell,
        const LCAO_Orbitals& orb,
		const ModulePW::PW_Basis& rhopw, 
		const ModulePW::PW_Basis_Big& bigpw)
{
    ModuleBase::TITLE("LCAO_domain","grid_prepare");
    ModuleBase::timer::tick("LCAO_domain","grid_prepare");
    if(PARAM.globalv.gamma_only_local)
    {
		gint_gamma.prep_grid(
				gt, 
				bigpw.nbx, 
				bigpw.nby, 
				bigpw.nbzp, 
				bigpw.nbzp_start,
				rhopw.nxyz, 
				bigpw.bx, 
				bigpw.by, 
				bigpw.bz, 
				bigpw.bxyz, 
				bigpw.nbxx,
				rhopw.ny, 
				rhopw.nplane, 
				rhopw.startz_current,
				&ucell,
				&orb);
	}
    else // multiple k-points
    {
        // cal the grid integration of 'Vl' matrix for l-points algorithms.
		gint_k.prep_grid(
				gt, 
				bigpw.nbx, 
				bigpw.nby, 
				bigpw.nbzp, 
				bigpw.nbzp_start,
				rhopw.nxyz, 
				bigpw.bx, 
				bigpw.by, 
				bigpw.bz, 
				bigpw.bxyz, 
				bigpw.nbxx,
				rhopw.ny, 
				rhopw.nplane, 
				rhopw.startz_current,
				&ucell,
				&orb);
	}

    ModuleBase::timer::tick("LCAO_domain","grid_prepare");
    return;
}

}

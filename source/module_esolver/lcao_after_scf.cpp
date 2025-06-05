#include "esolver_ks_lcao.h"

#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_io/cube_io.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/output_sk.h"
#include "module_io/write_elecstat_pot.h"
#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_io/print_info.h"

#include <memory>

//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"

// 2025-06-04
#include "module_io/ctrl_output_lcao.h"
#include <iostream>

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");

    //------------------------------------------------------------------
    //! 1) call after_scf() of ESolver_KS
    //------------------------------------------------------------------
    ESolver_KS<TK>::after_scf(ucell, istep, conv_esolver);

    //------------------------------------------------------------------
    //! 2) output of lcao every few ionic steps 
    //------------------------------------------------------------------
	const auto* estate = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec);
	if(!estate)
	{
		ModuleBase::WARNING_QUIT("ModuleIO::ctrl_output_lcao","pelec does not exist");
	}

    if(istep % PARAM.inp.out_interval == 0)
    {
/*
        ModuleIO::ctrl_output_lcao<TK, TR>(ucell, 
				this->kv,
				this->pelec, 
				this->pv, 
				this->gd,
				this->psi,
				this->p_hamilt,
				this->two_center_bundle_,
				this->GK,
				this->orb_,
				this->pw_wfc,
				this->pw_rho,
				this->GridT,
				this->pw_big,
				this->sf,
				this->rdmft_solver,
#ifdef __DEEPKS
				this->ld,
#endif
#ifdef __EXX
				this->exd,
				this->exc,
#endif
				istep);
*/
    }

    //------------------------------------------------------------------
    //! 3) Clean up RA, which is used to serach for adjacent atoms
    //------------------------------------------------------------------
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        RA.delete_grid();
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver

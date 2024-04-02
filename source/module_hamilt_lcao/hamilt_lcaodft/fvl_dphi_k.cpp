#include "FORCE_k.h"

#include <map>
#include <unordered_map>

#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


// calculate the force due to < phi | Vlocal | dphi >
void Force_LCAO_k::cal_fvl_dphi_k(const bool isforce,
		const bool isstress,
		LCAO_Matrix &lm,
        Gint_k &gint_k,
		const elecstate::Potential* pot_in,
		ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& svl_dphi,
		double** DM_R)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvl_dphi_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
        return;
    }

    assert(lm.DHloc_fixedR_x != NULL);
    assert(lm.DHloc_fixedR_y != NULL);
    assert(lm.DHloc_fixedR_z != NULL);

    int istep = 1;

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        const double* vr_eff1 = pot_in->get_effective_v(GlobalV::CURRENT_SPIN);
        const double* vofk_eff1 = nullptr;
        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            vofk_eff1 = pot_in->get_effective_vofk(GlobalV::CURRENT_SPIN);
        }

        //--------------------------------
        // Grid integration here.
        //--------------------------------
        // fvl_dphi can not be set to zero here if Vna is used
        if (isstress || isforce)
        {
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                Gint_inout inout(DM_R,
                                 is,
                                 vr_eff1,
                                 vofk_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force_meta);
                gint_k.cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(DM_R, is, vr_eff1, isforce, isstress, &fvl_dphi, &svl_dphi, Gint_Tools::job_type::force);
                gint_k.cal_gint(&inout);
            }
        }
    }

    if (isstress)
    {
        StressTools::stress_fill(-1.0, GlobalC::ucell.omega, svl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
    return;
}

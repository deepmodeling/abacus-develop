#include "FORCE.h"

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
template<>
void Force_LCAO<std::complex<double>>::cal_fvl_dphi(
    const bool isforce,
    const bool isstress,
    const elecstate::Potential* pot_in,
    TGint<std::complex<double>>::type& gint,
    ModuleBase::matrix& fvl_dphi,
    ModuleBase::matrix& svl_dphi)
{
    ModuleBase::TITLE("Force_LCAO", "cal_fvl_dphi");
    if (!isforce && !isstress)
    {
        return;
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_fvl_dphi");

    const int nspin = GlobalV::NSPIN;

    for (int is = 0; is < nspin; ++is)
    {
        const double* vr_eff1 = pot_in->get_effective_v(is);
        const double* vofk_eff1 = nullptr;

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            vofk_eff1 = pot_in->get_effective_vofk(is);
        }

        //--------------------------------
        // Grid integration here.
        //--------------------------------
        // fvl_dphi can not be set to zero here if Vna is used
        if (isstress || isforce)
        {
            // meta GGA functionals
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                Gint_inout inout(is,
                                 vr_eff1,
                                 vofk_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force_meta);

                gint.cal_gint(&inout);
            }
            else
            {
				Gint_inout inout(is, 
						vr_eff1, 
						isforce, 
						isstress, 
						&fvl_dphi, 
						&svl_dphi, 
						Gint_Tools::job_type::force);

                gint.cal_gint(&inout);
            }
        }
    }

    if (isstress)
    {
        StressTools::stress_fill(-1.0, GlobalC::ucell.omega, svl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_fvl_dphi");
    return;
}

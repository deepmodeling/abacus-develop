#pragma once
#include "pulay_force_stress.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_parameter/parameter.h"
namespace PulayForceStress
{
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const elecstate::Potential* pot, ///< [in] potential on grid
        typename TGint<TK>::type& gint,
        const bool& isforce,
        const bool& isstress,
        const bool& set_dmr_gint)
    {
        if (set_dmr_gint) { gint.transfer_DM2DtoGrid(dm.get_DMR_vector()); }    // 2d block to grid
        const int nspin = PARAM.inp.nspin;
        for (int is = 0; is < nspin; ++is)
        {
            const double* vr_eff1 = pot->get_effective_v(is);
            const double* vofk_eff1 = nullptr;
            if (XC_Functional::get_ked_flag())
            {
                vofk_eff1 = pot->get_effective_vofk(is);
                Gint_inout inout(is, vr_eff1, vofk_eff1, isforce, isstress, &f, &s, Gint_Tools::job_type::force_meta);
                gint.cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(is, vr_eff1, isforce, isstress, &f, &s, Gint_Tools::job_type::force);
                gint.cal_gint(&inout);
            }
        }
        if (isstress) { StressTools::stress_fill(-1.0, ucell.omega, s); }
    }
}
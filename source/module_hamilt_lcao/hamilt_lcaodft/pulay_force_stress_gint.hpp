#pragma once
#include "pulay_force_stress.h"
#include "module_hamilt_lcao/hamilt_lcaodft/stress_tools.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_interface.h"
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
        const int nspin = PARAM.inp.nspin;

#ifndef __NEW_GINT
        if (set_dmr_gint) { gint.transfer_DM2DtoGrid(dm.get_DMR_vector()); }    // 2d block to grid
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
#else
        std::vector<const double*> vr_eff(nspin, nullptr);
        std::vector<const double*> vofk_eff(nspin, nullptr);
        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            for (int is = 0; is < nspin; ++is)
            {
                vr_eff[is] = pot->get_effective_v(is);
                vofk_eff[is] = pot->get_effective_vofk(is);
            }
            ModuleGint::cal_gint_fvl_meta(nspin, vr_eff, vofk_eff, dm.get_DMR_vector(), isforce, isstress, &f, &s);
        }
        else
        {
            for(int is = 0; is < nspin; ++is)
            {
                vr_eff[is] = pot->get_effective_v(is);
            }
            ModuleGint::cal_gint_fvl(nspin, vr_eff, dm.get_DMR_vector(), isforce, isstress, &f, &s);
        }
#endif

        if (isstress) { StressTools::stress_fill(-1.0, ucell.omega, s); }
    }
}
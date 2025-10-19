#pragma once
#include "pulay_fs.h"
#include "source_lcao/stress_tools.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/temp_gint/gint_interface.h"
namespace PulayForceStress
{
    template<typename TK, typename TR>
    void cal_pulay_fs(
        ModuleBase::matrix& f,  ///< [out] force
        ModuleBase::matrix& s,  ///< [out] stress
        const elecstate::DensityMatrix<TK, TR>& dm,  ///< [in] density matrix
        const UnitCell& ucell,  ///< [in] unit cell
        const elecstate::Potential* pot, ///< [in] potential on grid
        const bool& isforce,
        const bool& isstress,
        const bool& set_dmr_gint)
    {
        const int nspin = PARAM.inp.nspin;

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

        if (isstress) { StressTools::stress_fill(-1.0, ucell.omega, s); }
    }
}

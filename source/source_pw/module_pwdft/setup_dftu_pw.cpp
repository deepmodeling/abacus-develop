#include "source_pw/module_pwdft/setup_dftu_pw.h"
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 2025-11-06
#include "source_pw/module_pwdft/dftu_base_io.h" // mohan add 2025-11-08
#include "source_io/module_parameter/parameter.h"

namespace DFTU_BASE
{

void iter_init_dftu_pw(const int iter,
                        const int istep,
                        Plus_U_Base& dftu, // mohan add 2025-11-06
                        const void* psi,
                        const ModuleBase::matrix& wg,
                        const UnitCell& ucell,
                        Charge_Mixing* p_chgmix,
                        const int* isk)
{
    if (!p_chgmix || !PARAM.inp.dft_plus_u)
    {
        return;
    }

    if (iter == 1 && istep == 0)
    {
        return;
    }

    if (dftu.get_occ_mat_ctrl() != 2)
    {
        dftu.cal_occ_pw(psi, wg, ucell, p_chgmix, isk);
    }
    DFTU_BASE::output(dftu, ucell, PARAM.inp.out_chg[0], PARAM.globalv.global_out_dir, PARAM.inp.nspin, PARAM.globalv.npol);
}

}

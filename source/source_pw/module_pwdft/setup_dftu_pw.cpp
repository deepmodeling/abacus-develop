#include "source_pw/module_pwdft/setup_dftu_pw.h"
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 2025-11-06
#include "source_pw/module_pwdft/dftu_output.h" // mohan add 2025-11-08
#include "source_io/module_parameter/parameter.h"

namespace pw
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

    // With occupation-matrix control 2 the matrix is intentionally frozen.
    // It still defines a DFT+U potential, so rebuild that potential before
    // the first Hamiltonian construction of every SCF step.  The ordinary
    // occupation-update path keeps its historical first-step guard.
    if (dftu.get_occ_mat_ctrl() == 2)
    {
        dftu.update_eff_pot_pw(ucell);
    }
    else if (!(iter == 1 && istep == 0))
    {
        dftu.cal_occ_pw(psi, wg, ucell, p_chgmix, isk);
    }
    dftu_io::output(dftu, ucell, PARAM.inp.out_chg[0], PARAM.globalv.global_out_dir, PARAM.inp.nspin, PARAM.globalv.npol);
}

}

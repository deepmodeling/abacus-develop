#include "source_pw/module_pwdft/dftu_pw.h"
#include "source_lcao/module_dftu/dftu.h"

namespace pw
{

void iter_init_dftu_pw(const int iter,
                       const int istep,
                       Plus_U& dftu,
                       const void* psi,
                       const ModuleBase::matrix& wg,
                       const UnitCell& ucell,
                       const Input_para& inp)
{
    if (!inp.dft_plus_u)
    {
        return;
    }

    if (iter == 1 && istep == 0)
    {
        return;
    }

    if (dftu.omc != 2)
    {
        dftu.cal_occ_pw(iter, psi, wg, ucell, inp.mixing_beta);
    }
    dftu.output(ucell);
}

}

#include <sstream>
#include <string>

#include "input.h"
#include "input_conv.h"
#include "module_parameter/parameter.h"

void Input_Conv::tmp_convert()
{
    INPUT.stru_file = PARAM.inp.stru_file;
    INPUT.pw_seed = PARAM.inp.pw_seed;
    INPUT.init_vel = PARAM.inp.init_vel;
    INPUT.ref_cell_factor = PARAM.inp.ref_cell_factor;
    INPUT.bessel_nao_rcut = PARAM.globalv.bessel_nao_rcut;
    INPUT.cond_dtbatch = PARAM.inp.cond_dtbatch;
    INPUT.nche_sto = PARAM.inp.nche_sto;
    INPUT.mdp = PARAM.mdp;

    const int ntype = PARAM.inp.ntype;
    delete[] INPUT.orbital_corr;
    INPUT.orbital_corr = new int[ntype];
    for (int i = 0; i < ntype; ++i)
    {
        INPUT.orbital_corr[i] = PARAM.inp.orbital_corr[i];
    }
    delete[] INPUT.hubbard_u;
    INPUT.hubbard_u = new double[ntype];
    for (int i = 0; i < ntype; ++i)
    {
        INPUT.hubbard_u[i] = PARAM.globalv.hubbard_u[i];
    }
    INPUT.bands_to_print_ = PARAM.inp.bands_to_print;
}

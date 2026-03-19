#include "source_pw/module_pwdft/update_cell_pw.h"
#include "source_base/global_variable.h"
#include "source_base/global_function.h"

namespace pw
{

void update_cell_pw(const UnitCell& ucell,
                    pseudopot_cell_vnl& ppcell,
                    const K_Vectors& kv,
                    ModulePW::PW_Basis_K* pw_wfc,
                    const Input_para& inp)
{
    ModuleBase::TITLE("pw", "update_cell_pw");

    if (!ucell.cell_parameter_updated)
    {
        return;
    }

    ppcell.rescale_vnl(ucell.omega);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    pw_wfc->initgrids(ucell.lat0, ucell.latvec, pw_wfc->nx, pw_wfc->ny, pw_wfc->nz);
    pw_wfc->initparameters(false, inp.ecutwfc, kv.get_nks(), kv.kvec_d.data());
    pw_wfc->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
}

}

#include "param_update.h"
#include "source_io/module_parameter/parameter.h"

namespace elecstate {

void ParamUpdater::update_from_atoms_info(const AtomsInfoResult& atoms_info)
{
    PARAM.input.nelec = atoms_info.nelec;
    PARAM.input.nbands = atoms_info.nbands;
    PARAM.input.nupdown = atoms_info.nupdown;
    PARAM.sys.use_uspp = atoms_info.use_uspp;
    PARAM.sys.nbands_l = atoms_info.nbands_l;
    PARAM.sys.ks_run = atoms_info.ks_run;
}

}
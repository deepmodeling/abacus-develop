#include "elecstate.h"
#include "source_base/parallel_reduce.h"
#include "source_estate/module_pot/h_hartree_pw.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"

namespace elecstate
{

double ElecState::get_hartree_energy()
{
    return H_Hartree_pw::hartree_energy;
}

double ElecState::get_etot_efield()
{
    return Efield::etotefield;
}

double ElecState::get_etot_gatefield()
{
    return Gatefield::etotgatefield;
}

double ElecState::get_solvent_model_Ael()
{
    return surchem::Ael;
}

double ElecState::get_solvent_model_Acav()
{
    return surchem::Acav;
}

double ElecState::get_dftu_energy()
{
    return dftu_energy_;
}

double ElecState::get_local_pp_energy()
{
    double local_pseudopot_energy = 0.; // electron-ion interaction energy from local pseudopotential
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        local_pseudopot_energy
            += BlasConnector::dot(this->charge->rhopw->nrxx, 
                                  this->pot->get_fixed_v(), 
                                  1, 
                                  this->charge->rho[is], 1)
                                  * this->charge->rhopw->omega / this->charge->rhopw->nxyz;
    }
    Parallel_Reduce::reduce_pool(local_pseudopot_energy);
    return local_pseudopot_energy;
}

} // namespace elecstate

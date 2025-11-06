#include "elecstate.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_lcao/module_dftu/dftu.h" // mohan add 2025-11-06

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
    return Plus_U::get_energy();
}

double ElecState::get_local_pp_energy()
{
    double local_pseudopot_energy = 0.; // electron-ion interaction energy from local pseudopotential
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        local_pseudopot_energy
            += BlasConnector::dot(this->charge->rhopw->nrxx, this->pot->get_fixed_v(), 1, this->charge->rho[is], 1)
               * this->charge->rhopw->omega / this->charge->rhopw->nxyz;
    }
    Parallel_Reduce::reduce_all(local_pseudopot_energy);
    return local_pseudopot_energy;
}

} // namespace elecstate

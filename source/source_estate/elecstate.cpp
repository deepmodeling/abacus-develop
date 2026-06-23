#include "elecstate.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "source_estate/module_charge/charge.h"
#include "module_pot/potential_new.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "occupy.h"

namespace elecstate
{

ElecState::ElecState(Charge* chr_in, ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis_Big* bigpw_in)
{
    this->charge = chr_in;
    this->charge->set_rhopw(rhopw_in);
    this->bigpw = bigpw_in;
    this->eferm.two_efermi = PARAM.globalv.two_fermi;
}

ElecState::~ElecState()
{
    if (this->pot != nullptr)
    {
        delete this->pot;
        this->pot = nullptr;
    }
}

const double* ElecState::getRho(int spin) const
{
    return &(this->charge->rho[spin][0]);
}



void ElecState::init_nelec_spin()
{
    this->nelec_spin.resize(PARAM.inp.nspin);
    if (PARAM.inp.nspin == 2)
    {
        this->nelec_spin[0] = (PARAM.inp.nelec + PARAM.inp.nupdown) / 2.0;
        this->nelec_spin[1] = (PARAM.inp.nelec - PARAM.inp.nupdown) / 2.0;
    }
}


void ElecState::init_ks(Charge* chr_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        const ModulePW::PW_Basis_Big* bigpw_in)
{
    this->charge = chr_in;
    this->klist = klist_in;
    this->bigpw = bigpw_in;
    // init nelec_spin with nelec and nupdown
    this->init_nelec_spin();
    // initialize ekb and wg
    this->ekb.create(nk_in, PARAM.globalv.nbands_l);
    this->wg.create(nk_in, PARAM.globalv.nbands_l);
}

} // namespace elecstate

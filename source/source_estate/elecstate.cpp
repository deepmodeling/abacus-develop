#include "elecstate.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/memory.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_title.h"
#include "occupy.h"

namespace elecstate
{

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

void ElecState::init_scf(const int istep, 
                         const UnitCell& ucell,
                         const Parallel_Grid& pgrid,
                         const ModuleBase::ComplexMatrix& strucfac, 
                         const bool* numeric,
                         ModuleSymmetry::Symmetry& symm, 
                         const void* wfcpw)
{
    //! core correction potential.
    this->charge->set_rho_core(ucell,strucfac, numeric);

    //! other effective potentials need charge density,
    // choose charge density from ionic step 0.
    if (istep == 0)
    {
        this->charge->init_rho(this->eferm,ucell, pgrid, strucfac, symm, (const void*)this->klist, wfcpw);
        this->charge->check_rho(); // check the rho
    }

    //! renormalize the charge density
    this->charge->renormalize_rho();

    //! initialize the potential
    this->pot->init_pot(istep, this->charge);
}


void ElecState::init_ks(Charge* chg_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        ModulePW::PW_Basis* rhopw_in,
                        const ModulePW::PW_Basis_Big* bigpw_in)
{
    this->charge = chg_in;
    this->klist = klist_in;
    this->charge->set_rhopw(rhopw_in);
    this->bigpw = bigpw_in;
    // init nelec_spin with nelec and nupdown
    this->init_nelec_spin();
    // initialize ekb and wg
    this->ekb.create(nk_in, PARAM.globalv.nbands_l);
    this->wg.create(nk_in, PARAM.globalv.nbands_l);
}

} // namespace elecstate

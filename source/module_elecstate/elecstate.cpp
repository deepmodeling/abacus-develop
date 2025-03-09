#include "elecstate.h"
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/tool_title.h"
#include "occupy.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    return &(this->charge->rho[spin][0]);
}

void ElecState::fixed_weights(const std::vector<double>& ocp_kb, const int& nbands, const double& nelec)
{
    assert(nbands > 0);
    assert(nelec > 0.0);

    const double ne_thr = 1.0e-5;

    const int num = this->klist->get_nks() * nbands;
    if (num != ocp_kb.size())
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                 "size of occupation array is wrong , please check ocp_set");
    }

    double num_elec = 0.0;
    for (int i = 0; i < ocp_kb.size(); ++i)
    {
        num_elec += ocp_kb[i];
    }

    if (std::abs(num_elec - nelec) > ne_thr)
    {
        ModuleBase::WARNING_QUIT("ElecState::fixed_weights",
                                 "total number of occupations is wrong , please check ocp_set");
    }

    for (int ik = 0; ik < this->wg.nr; ++ik)
    {
        for (int ib = 0; ib < this->wg.nc; ++ib)
        {
            this->wg(ik, ib) = ocp_kb[ik * this->wg.nc + ib];
        }
    }
    this->skip_weights = true;

    return;
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
    if (!PARAM.inp.use_paw)
    {
        this->charge->set_rho_core(ucell,strucfac, numeric);
    }
    else
    {
        this->charge->set_rho_core_paw();
    }

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

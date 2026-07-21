#include "elecstate.h"
#include "source_io/module_chgpot/write_init.h"

namespace elecstate
{

void init_scf(const UnitCell& ucell,
              const Parallel_Grid& pgrid,
              const ModuleBase::ComplexMatrix& strucfac, 
              const bool* numeric,
              const int istep,
              const std::string& out_dir,
              const Input_para& inp,
              ElecState* pelec)
{
    //! core correction potential.
    pelec->charge->set_rho_core(ucell, strucfac, numeric);

    //! renormalize the charge density
    if(PARAM.inp.init_chg != "dm_no_renormalize")
    {
        pelec->charge->renormalize_rho();
    }

    //! initialize the potential
    pelec->pot->init_pot(pelec->charge);

    //! output the initial potential
    ModuleIO::write_pot_init(ucell, pgrid, pelec, istep, out_dir, inp, PARAM.globalv.two_fermi);
}

} // namespace elecstate

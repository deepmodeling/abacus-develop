#include "Cell_PW.h"
#include "ions.h"
#include "../module_base/timer.h"

void Cell_PW::opt_cells_pw(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Cell_PW", "opt_cells_pw");
    ModuleBase::timer::tick("Cell_PW", "opt_cells_pw");


    // ion optimization begins
    // electron density optimization is included in ion optimization

    Ions ions;
    ions.opt_ions_pw(p_esolver);
    
    ModuleBase::timer::tick("Cell_PW", "opt_cells_pw");
}

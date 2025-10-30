#ifndef INIT_DM_H
#define INIT_DM_H

#include "source_cell/unitcell.h" // use unitcell
#include "source_estate/elecstate_lcao.h"// use ElecStateLCAO
#include "source_psi/psi.h" // use electronic wave functions
#include "source_estate/module_charge/charge.h" // use charge

namespace elecstate
{

template <typename TK> 
void init_dm(UnitCell& ucell,
		ElecStateLCAO<TK>* pelec,
        psi::Psi<TK>* psi,
		Charge &chr,
        const int iter,
        const int exx_two_level_step);

}

#endif

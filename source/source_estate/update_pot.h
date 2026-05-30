#ifndef UPDATE_POT_H
#define UPDATE_POT_H

#include "source_cell/unitcell.h"
#include "source_estate/elecstate.h"

namespace elecstate
{

void update_pot(UnitCell& ucell, // unitcell 
		elecstate::ElecState* &pelec, // pointer of electrons
		const Charge &chr,
        const bool conv_esolver); // charge density
}


#endif

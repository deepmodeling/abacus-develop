#include "source_estate/update_pot.h"
#include "source_estate/cal_ux.h"

void elecstate::update_pot(UnitCell& ucell, // unitcell 
		elecstate::ElecState* &pelec, // pointer of electrons
		const Charge &chr,
        const bool conv_esolver
         ) // charge density
{
    if (!conv_esolver)
    {
        elecstate::cal_ux(ucell);
        pelec->pot->update_from_charge(&chr, &ucell);
        pelec->f_en.descf = pelec->cal_delta_escf();
    }
    else
    {
        pelec->cal_converged();
    }
}

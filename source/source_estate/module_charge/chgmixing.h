#ifndef CHGMIXING_H
#define CHGMIXING_H

#include "source_estate/module_charge/charge_mixing.h"
#include "source_io/module_parameter/input_parameter.h" // use Input_para

namespace module_charge
{

void chgmixing(const int iter,
        Charge_Mixing* p_chgmix,
		const Input_para& inp); // input parameters

}


#endif

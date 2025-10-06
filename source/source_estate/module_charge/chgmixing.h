#ifndef CHGMIXING_H
#define CHGMIXING_H

#include "source_estate/module_charge/charge_mixing.h"
#include "source_io/module_parameter/input_parameter.h" // use Input_para

namespace module_charge
{

void chgmixing_ks_pw(const int iter,
        Charge_Mixing* p_chgmix,
		const Input_para& inp); // input parameters

void chgmixing_ks_lcao(const int iter, // scf iteration number
        Charge_Mixing* p_chgmix, // charge mixing class
        const int nnr, // dimension of density matrix
		const Input_para& inp); // input parameters

}


#endif

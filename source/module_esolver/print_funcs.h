#ifndef PRINT_FUNCTIONS_H
#define PRINT_FUNCTIONS_H

#include "module_parameter/input_parameter.h"
#include "module_basis/module_pw/pw_basis_k.h"

namespace Print_functions
{
    void print_wfcfft(
		const Input_para& inp, 
        ModulePW::PW_Basis_K& pw_wfc,
		std::ofstream& ofs);
}	

#endif

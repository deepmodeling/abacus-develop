#ifndef WRITE_WFC_PW_H
#define WRITE_WFC_PW_H
#include "module_psi/psi.h"
#include "module_pw/pw_basis_k.h"
#include "module_cell/klist.h"

namespace ModuleIO
{
	void write_wfc_pw( const std::string &fn, const psi::Psi<std::complex<double>> &psi, const K_Vectors* p_kv, 
                    const ModulePW::PW_Basis_K *wfc_basis);
}

#endif

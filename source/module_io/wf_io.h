#ifndef WF_IO_H
#define WF_IO_H
#include "module_psi/psi.h"
#include "module_pw/pw_basis_k.h"
#include "src_pw/klist.h"

namespace WF_io
{
	void write_wfc( const std::string &fn, const psi::Psi<std::complex<double>> &psi, const K_Vectors* p_kv, 
                    const ModulePW::PW_Basis_K *wfc_basis);
    void read_wfc(  const std::string &fn, const psi::Psi<std::complex<double>> &psi, const K_Vectors* p_kv, 
                    const ModulePW::PW_Basis_K *wfc_basis);
}

#endif

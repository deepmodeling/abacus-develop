#ifndef WRITE_WFC_R_H
#define WRITE_WFC_R_H

#ifdef __MPI
#include "mpi.h"
#endif

#include <complex>
#include <string>
#include <vector>

#include "source_base/complexmatrix.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "module_psi/psi.h"

namespace ModuleIO
{
	// write ||wfc_r|| for all k-points and all bands
	// Input: wfc_g[ik](ib,ig)
	// loop order is for(z){for(y){for(x)}}
void write_psi_r_1(const UnitCell& ucell,
                   const psi::Psi<std::complex<double>>& wfc_g,
                   const ModulePW::PW_Basis_K* wfcpw,
                   const std::string& folder_name,
                   const bool& square,
                   const K_Vectors& kv);

// Input: wfc_g(ib,ig)
// Output: wfc_r[ir]
std::vector<std::complex<double>> cal_wfc_r(const ModulePW::PW_Basis_K* wfcpw,
                                            const psi::Psi<std::complex<double>>& wfc_g,
                                            const int ik,
                                            const int ib);

// Input: chg_r[ir]
void write_chg_r_1(const UnitCell& ucell,
                   const ModulePW::PW_Basis_K* wfcpw,
                   const std::vector<double>& chg_r,
                   const std::string& file_name
                   #ifdef __MPI
                   ,MPI_Request& mpi_request
                   #endif
                   );
}

#endif

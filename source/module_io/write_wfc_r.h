//======================
// AUTHOR : Peize Lin
// DATE :   2021-11-21
//======================

#ifndef WRITE_WFC_R_H
#define WRITE_WFC_R_H

#ifdef __MPI
#include "mpi.h"
#endif

#include <complex>
#include <string>
#include <vector>

#include "module_base/complexmatrix.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_psi/psi.h"

namespace ModuleIO
{
	// write ||wfc_r|| for all k-points and all bands
	// Input: wfc_g[ik](ib,ig)
	// loop order is for(z){for(y){for(x)}}
void write_psi_r_1(const psi::Psi<std::complex<double>>& wfc_g,
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
#ifdef __MPI
void write_chg_r_1(const ModulePW::PW_Basis_K* wfcpw,
                   const std::vector<double>& chg_r,
                   const std::string& file_name,
                   MPI_Request& mpi_request);
#else
void write_chg_r_1(const ModulePW::PW_Basis_K* wfcpw, const std::vector<double>& chg_r, const std::string& file_name);
#endif
}

#endif

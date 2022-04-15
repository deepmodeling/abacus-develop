//======================
// AUTHOR : Peize Lin
// DATE :   2021-11-21
//======================

#ifndef WRITE_WFC_REALSPACE_H
#define WRITE_WFC_REALSPACE_H

#include "module_base/complexmatrix.h"
#include "module_base/vector3.h"
#include <vector>
#include <complex>
#include <string>

#ifdef __MPI
#include "mpi.h"
#endif

namespace Write_Wfc_Realspace
{
	// write ||wfc_r|| for all k-points and all bands
	// Input: wfc_g[ik](ib,ig)
	// loop order is for(z){for(y){for(x)}}
    void write_wfc_realspace_1(const ModuleBase::ComplexMatrix*const wfc_g, const std::string &folder_name, const bool& square);

	// Input: wfc_g(ib,ig)
	// Output: wfc_r[ir]
    std::vector<std::complex<double>> cal_wfc_r(const ModuleBase::ComplexMatrix &wfc_g, const int ik, const int ib);

	// Input: chg_r[ir]
#ifdef __MPI
	void write_charge_realspace_1(const std::vector<double> &chg_r, const std::string &file_name, MPI_Request &mpi_request);
#else
	void write_charge_realspace_1(const std::vector<double> &chg_r, const std::string &file_name);
#endif
}

#endif
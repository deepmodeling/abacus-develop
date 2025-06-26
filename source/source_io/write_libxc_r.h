//======================
// AUTHOR : Peize Lin
// DATE :   2024-09-12
//======================

#ifndef WRITE_LIBXC_R_H
#define WRITE_LIBXC_R_H

#ifdef USE_LIBXC

#include <vector>
#include <fstream>

	class Charge;
	namespace ModulePW{ class PW_Basis_Big; }
	namespace ModulePW{ class PW_Basis; }

namespace ModuleIO
{
	extern void write_libxc_r(
		const int order,
		const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge &chr,
		const ModulePW::PW_Basis_Big &pw_big,
		const ModulePW::PW_Basis &pw_rhod);
	
  #ifdef __MPI
	extern void write_cube_core(
		std::ofstream &ofs_cube,
		const int bz,
		const int nbz,
		const int nplane,
		const int startz_current,
		const double*const data,
		const int nxy,
		const int nz,
		const int nld,
		const int n_data_newline);
  #else
	extern void write_cube_core(
		std::ofstream &ofs_cube,
		const double*const data,
		const int nxy,
		const int nz,
		const int n_data_newline);
  #endif
}

#endif // USE_LIBXC

#endif // WRITE_LIBXC_R_H
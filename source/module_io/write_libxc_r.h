//======================
// AUTHOR : Peize Lin
// DATE :   2024-09-12
//======================

#ifndef WRITE_LIBXC_R_H
#define WRITE_LIBXC_R_H

#ifdef USE_LIBXC

#include <vector>

	class Charge;
	namespace ModuleESolver{ class ESolver_FP; }

namespace ModuleIO
{
	extern void write_libxc_r(
		const int order,
		const std::vector<int> &func_id,
		const int &nrxx, // number of real-space grid
		const double &omega, // volume of cell
		const double tpiba,
		const Charge* const chr,
		const ModuleESolver::ESolver_FP*const esolver);
}

#endif // USE_LIBXC

#endif // WRITE_LIBXC_R_H
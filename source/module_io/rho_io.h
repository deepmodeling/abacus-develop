#ifndef RHO_IO_H
#define RHO_IO_H
#include <string>
#include "module_cell/unitcell.h"
#ifdef __MPI
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#endif

namespace ModuleIO
{
	bool read_rho(
#ifdef __MPI
		Parallel_Grid* Pgrid,
#endif
		const int &is,
		const int &nspin,
		const std::string &fn,
		double* rho,
		int& nx,
		int& ny,
		int& nz,
		double& ef,
		UnitCell& ucell,
		int &prenspin);

	void write_rho(const double* rho_save,
		    const int &is,
		    const int &iter,
		    const std::string &fn,
		    const int &precision = 11,
		    const bool for_plot = false);//mohan add 2007-10-17
}

#endif

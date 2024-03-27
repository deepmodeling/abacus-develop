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
		int my_rank,
		std::string esolver_type,
		int rank_in_stogroup,
		const int& is,
		std::ofstream& ofs_running,
		const int &nspin,
		const std::string &fn,
		double* rho,
		int& nx,
		int& ny,
		int& nz,
		double& ef,
		const UnitCell* ucell,
		int &prenspin);

	void write_rho(
#ifdef __MPI
		const int& bz,
		const int& nbz,
		const int& nplane,
		const int& startz_current,
#endif
		const double* rho_save,
		const int& is,
		const int& nspin,
		const int& iter,
		const std::string& fn,
		const int& nx,
		const int& ny,
		const int& nz,
		const double& ef,
		const UnitCell* ucell,
		const int &precision = 11);//mohan add 2007-10-17
}

#endif
